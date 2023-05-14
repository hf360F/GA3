"""
EXTERNAL CONSTRAINTS:
Total mass < 1 kg
Total length < 0.35 m
Total copper tube length < 3.5 m

Compressor characteristics fixed (see csv)
Cold stream inlet fixed as 20 C
Hot stream inlet fixed as 60 C
"""

import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import least_squares, root


def K_2d_interp(sigma, Ret):
    dfKc = pd.read_csv("data/Turb_Kc.csv").to_numpy()
    dfKe = pd.read_csv("data/Turb_Ke.csv").to_numpy()

    Ret = np.clip(Ret, 3000, 10000)
    sigma = np.clip(sigma, 0, 1)

    Kc = scipy.interpolate.griddata(dfKc[:, :-1], dfKc[:, -1], np.column_stack((Ret, sigma)), 'cubic')
    Ke = scipy.interpolate.griddata(dfKe[:, :-1], dfKe[:, -1], np.column_stack((Ret, sigma)), 'cubic')

    return Kc + Ke


def K(sigma, Ret):
    """Entrance and exit loss coefficient for a tube.  Assumes turbulent flow in tube.
       Will warn for out of range Reynolds number, but match to the closest curve.

    Args:
        sigma (float): Ratio of total tube area to shell area.
        Ret (float): Tube Reynolds number.

    Returns:
        float: Sum of Kc and Ke loss coefficients.
    """

    polyOrder = 3

    # Reynolds range covered by curves
    Remin = 3000
    Re2 = 5000
    Remax = 10000

    interpolationFlag = False

    if Ret > Remax:
        dfKc = pd.read_csv("data/Turb_10000_Kc.csv")
        dfKe = pd.read_csv("data/Turb_10000_Ke.csv")
        warnings.warn(f"Tube Re = {Ret:.2f} out of range (Remax = {Remax:.2f}), matching to closest curve.")

    elif Ret > Re2:
        dfKc = pd.read_csv("data/Turb_5000_Kc.csv")
        dfKc2 = pd.read_csv("data/Turb_10000_Kc.csv")
        dfKe = pd.read_csv("data/Turb_5000_Ke.csv")
        dfKe2 = pd.read_csv("data/Turb_10000_Ke.csv")

        interpolationFlag = False  # Can interpolate between two curves of different Re

    elif Ret >= Remin:
        dfKc = pd.read_csv("data/Turb_3000_Kc.csv")
        dfKc2 = pd.read_csv("data/Turb_5000_Kc.csv")
        dfKe = pd.read_csv("data/Turb_3000_Ke.csv")
        dfKe2 = pd.read_csv("data/Turb_5000_Ke.csv")

        interpolationFlag = True

    else:
        dfKc = pd.read_csv("data/Turb_3000_Kc.csv")
        dfKe = pd.read_csv("data/Turb_3000_Ke.csv")
        warnings.warn(f"Tube Re = {Ret:.2f} out of range (Remin = {Remin:.2f}), matching to closest curve.")

    # Fit polynomials to digitised data
    polyCoeffsKc = np.polyfit(dfKc["sigma"], dfKc["Kc"], polyOrder)
    polyCoeffsKe = np.polyfit(dfKe["sigma"], dfKe["Ke"], polyOrder)

    if interpolationFlag:  # If we need second set of polynomials for interpolating between two curves
        polyCoeffsKc2 = np.polyfit(dfKc2["sigma"], dfKc2["Kc"], polyOrder)
        polyCoeffsKe2 = np.polyfit(dfKe2["sigma"], dfKe2["Ke"], polyOrder)

        if Ret >= Re2:
            weightingFactor = (Ret - Re2) / (Remax - Re2)
        else:
            weightingFactor = (Ret - Remin) / (Re2 - Remin)

        Kc = (1 - weightingFactor) * np.poly1d(polyCoeffsKc)(sigma) + weightingFactor * np.poly1d(polyCoeffsKc2)(sigma)
        Ke = (1 - weightingFactor) * np.poly1d(polyCoeffsKe)(sigma) + weightingFactor * np.poly1d(polyCoeffsKe2)(sigma)

    else:
        Kc = np.poly1d(polyCoeffsKc)(sigma)
        Ke = np.poly1d(polyCoeffsKe)(sigma)

    return Kc + Ke


def chicSolver(hx, pump):
    """Find intersection of pump and HX characteristic to set operating point.
    The HX flow path to use is inferred from pump.pump_type.
    
    Args:
        hx (HX object): Heat exchanger design.
        pump (Pump object): Pump design.

    Returns:
        (float, float): Mass flow rate, kg/s, and HX pressure drop (= pump pressure rise), Pa
    """

    if pump.pump_type is Pump.HOT:
        rho = hx.hotStream["rho"]
        hx_dp = hx.hydraulicAnalysisTube

        # We have found the intersection of the curves when pressure change over pump and HX sum to zero:
    else:
        rho = hx.coldStream["rho"]
        hx_dp = hx.hydraulicAnalysisShell

    def f(mdot):
        return hx_dp(mdot) - pump.dp(mdot / rho)

    solution = root(f, x0=0.2 * np.ones(len(hx.Nt)))

    if not solution.success:
        raise ValueError("Unable to intersect pump and heat exchanger characteristics!")
    else:
        return solution.x, pump.dp(solution.x / rho)


class HX:
    def __init__(self, coldStream, hotStream, kt, epst, lt, do, di, Nt, Y, isSquare, Np, Nb, B, G, ds, dn):
        """Heat exchanger design class.

        Args:
            coldStream (dict): Dictionary of cold inlet temperature 'Ti' C, specific heat capacity 'cp' J/kgK,
                fluid density 'rho' kg/m^3, thermal conductivity 'k' W/mK, viscocity 'mu' kg/ms.
            hotStream (dict): Dictionary of hot stream properties, as for cold stream.
            kt (float): Thermal conducitivty of tube, W/mK.
            epst (float): Effective roughness height of tube, m.
            lt (float|ndarray): Length of tube section, m.
            do (float): Tube outer diameter, m.
            di (float): Tube inner diameter, m.
            Nt (int|ndarray): Number of tubes.
            Y (float|ndarray): Tube pitch, m.
            isSquare (bool): Are tubes arranged in square pattern, if False tubes are asssumed to be in triangular pattern.
            Np (int): Number of passes.
            Nb (int): Number of baffles.
            B (float): Baffle pich, m.
            G (float): Baffle cut, m.
            ds (float): Shell diameter, m.
            dn (float): Nozzle diameter for both streams, m.
        """

        self.coldStream = coldStream
        self.hotStream = hotStream
        self.kt = kt
        self.epst = epst
        self.lt = lt
        self.do = do
        self.di = di
        self.Nt = Nt
        self.Y = Y
        self.isSquare = isSquare
        self.Np = Np
        self.Nb = Nb
        self.B = B
        self.G = G
        self.ds = ds
        self.dn = dn

        self.Attot = self.Nt * self.di ** 2 * np.pi / 4  # Tube total flowpath area, m^2
        self.Apipe = self.ds ** 2 * np.pi / 4
        self.sigma = self.Attot * self.Np / self.Apipe  # Note scaling with number of passes Np
        self.An = self.dn ** 2 * np.pi / 4
        self.F = 1
        if Np != 1:
            raise (NotImplementedError("F factor for multi-pass setups not implemented yet"))

    def hydraulicAnalysisTube(self, mdot, verbose=False):
        """Perform pressure drop analysis on tube flow path for given mdot.

        Args:
            mdot (float): Tube mass flow rate, kg/s. 
            verbose (bool): Whether or not to print intermediate values (velocities, areas, etc).       
        """

        Vt = mdot / (self.hotStream["rho"] * self.Attot)  # Bulk tube velocity, m/s
        Ret = self.hotStream["rho"] * Vt * self.di / self.hotStream["mu"]  # Tube Reynolds

        # Haaland's approximation of Colebrook-White for Darcy friction factor
        fTube = (-1.8 * np.log10((self.epst / (self.di * 3.7)) ** 1.11 + (6.9 / Ret))) ** (-2)
        dpFric = fTube * (self.lt / self.di) * 0.5 * self.hotStream["rho"] * (Vt ** 2)

        # Look up total loss factor, and calculate friction
        Ktot = K_2d_interp(self.sigma, Ret)
        dpEnds = Ktot * 0.5 * self.hotStream["rho"] * Vt ** 2

        # Nozzle loss
        Vn = mdot / (self.hotStream["rho"] * self.An)
        dpNozzles = 2 * 0.5 * self.hotStream["rho"] * (Vn ** 2)

        tubeTotaldP = dpFric + dpEnds + dpNozzles

        if verbose:
            print(f"\nTUBE HYDRAULIC ANALYSIS SUMMARY FOR mdot = {mdot:.2f} kg/s\n")
            print(f"Tube flow area: {self.Attot:.6f} m^2")
            print(f"Tube bulk velocity: {Vt:.2f} m/s")
            print(f"Tube Reynolds number: {Ret:.0f}")
            print(f"Tube friction factor: {fTube:.6f}")
            print(f"Total inlet/exit loss factor: {Ktot:.3f}")
            print(f"Total pressure drop: {tubeTotaldP:.0f} Pa (friction {dpFric:.0f},"
                  f" ends {dpEnds:.0f}, nozzles {dpNozzles:.0f})\n")

        return tubeTotaldP

    def hydraulicAnalysisShell(self, mdot, verbose=False):
        """Perform pressure drop analysis on shell flow path for given mdot.

        Args:
            mdot (float): Shell mass flow rate, kg/s. 
            verbose (bool): Whether or not to print intermediate values (velocities, areas, etc).       
        """

        # Shell loss
        self.As = self.ds * (self.Y - self.do) * self.B / self.Y  # Approximation of shell flow area
        Vs = mdot / (self.coldStream["rho"] * self.As)

        self.dseff = self.ds * (self.As / self.Apipe)
        Res = self.coldStream["rho"] * Vs * self.dseff / self.coldStream["mu"]

        if self.isSquare:
            a = 0.34
        else:
            a = 0.2

        shelldp1 = 4 * a * (Res ** (-0.15)) * self.Nt * self.coldStream["rho"] * (Vs ** 2)

        # Nozzle loss
        self.An = self.dn ** 2 * np.pi / 4
        Vn = mdot / (self.hotStream["rho"] * self.An)
        dpNozzles = 2 * 0.5 * self.coldStream["rho"] * (Vn ** 2)

        shellTotaldp = shelldp1 + dpNozzles

        if verbose:
            print(f"\nSHELL HYDRAULIC ANALYSIS SUMMARY FOR mdot = {mdot:.2f} kg/s\n")
            print(f"Shell flow area: {self.As:.6f} m^2")
            print(f"Shell characteristic velocity: {Vs:.2f} m/s")
            print(f"Shell Reynolds number: {Res:.0f}")
            print(f"Total pressure drop: {shellTotaldp:.0f} Pa (shell {shelldp1:.0f},"
                  f"  nozzles {dpNozzles:.0f})\n")

        return shellTotaldp

    def plotHXChics(self, mdotMin=0, mdotMax=1, n=100):
        """Plot HX flow characteristics.

        Args:
            mdotMin (float): Minimum mass flow to evaluate pressure drop at, kg/s.
            mdotMax (float): Maximum mass flow to evaluate pressure drop at, kg/s.
            n (int): Number of mass flow rates to evaluate between mdotMin and mdotMax.
        """

        mdots = np.linspace(mdotMin, mdotMax, n + 1)
        dpsTube, dpsShell = np.zeros_like(mdots), np.zeros_like(mdots)

        for i in range(n + 1):
            dpsTube[i] = self.hydraulicAnalysisTube(mdot=mdots[i])
            dpsShell[i] = self.hydraulicAnalysisShell(mdot=mdots[i])

        plt.plot(mdots, dpsTube / 1000, label="Tube flow path", color="red")
        plt.plot(mdots, dpsShell / 1000, label="Shell flow path", color="blue")
        plt.xlim(mdotMin, mdotMax)
        plt.ylim(0)
        plt.xlabel("Mass flow rate, $kg/s$")
        plt.ylabel("Pressure drop, $kPa$")
        plt.title("Heat exchanger characteristics")
        plt.legend()
        plt.grid()
        plt.show()

    def thermalAnalysis(self, mdot_t, mdot_s):
        """
        Perform thermal analysis on the HX given the 2 mass flowrates.

        Args:
            mdot_t: Tube mass flow, kg/s.
            mdot_s: Shell mass flow rate, kg/s

        Returns:
            Q: Resultant heat transfer rate, W.
        """
        # Reynolds numbers
        Ret = mdot_t * self.di / self.Attot / self.hotStream["mu"]
        Res = mdot_s * self.ds * (self.As / self.Apipe) / self.As / self.coldStream["mu"]

        # Pr and geometry constant are const
        Pr = 4.31
        c = 0.15 if self.isSquare else 0.2

        # Nusselt number correlations
        Nut = 0.023 * Ret ** 0.8 * Pr ** 0.3
        Nus = c * Res ** 0.6 * Pr ** 0.3

        # convection ht coeffs
        hi = Nut * self.hotStream["k"] / self.di
        ho = Nus * self.coldStream["k"] / self.do

        # total ht coeff
        H = 1 / (1 / hi + self.di * np.log(self.do / self.di) / 2 / self.kt + self.di / self.do / ho)
        # H = 9878

        # inlet temps
        Thi = self.hotStream["Ti"] * np.ones_like(mdot_t)
        Tci = self.coldStream["Ti"] * np.ones_like(mdot_s)

        # SFEE
        HA = H * np.pi * self.di * self.lt * self.Nt

        def LMTD(Tho_, Tco_):
            T1 = Thi - Tco_
            T2 = Tho_ - Tci
            mask = np.isclose(T1, T2)
            return mask * T1 + (1 - mask) * (T1 - T2) / np.log(T1 / T2)

        def f(To):
            Tho, Tco = np.split(To, 2)
            Q = HA * LMTD(Tho, Tco) * self.F
            return np.concatenate(
                (mdot_t * self.hotStream["cp"] * (Thi - Tho) - Q, mdot_s * self.coldStream["cp"] * (Tco - Tci) - Q))

        x0 = np.concatenate((40 * np.ones_like(mdot_s), 30 * np.ones_like(mdot_s)))

        res = least_squares(f, x0, bounds=(np.tile(Tci, 2), np.tile(Thi,2)))
        Tho, Tco = np.split(res.x,2)
        print(f'Tho: {Tho}\n'
              f'Tco: {Tco}')

        print(res)

        Q = HA * LMTD(Tho, Tco)

        return Q


class Pump:
    COLD, HOT = range(2)

    def __init__(self, pump_type: int) -> None:
        """
        Pump characteristics class.
        :param pump_type: The pump variety, one of Pump.COLD and Pump.HOT.
        """
        self.pump_type = pump_type
        match self.pump_type:
            case Pump.COLD:
                type_str = 'cold'
            case Pump.HOT:
                type_str = 'hot'
            case _:
                raise "invalid pump type"

        data = np.genfromtxt(f"data/{type_str}.csv", delimiter=',')
        self.flowrate_data = data[:, 0]
        self.flowMin = np.min(data[:, 0])
        self.flowMax = np.max(data[:, 0])
        self.dp_data = data[:, 1]
        self.poly = np.poly1d(np.polyfit(self.flowrate_data, self.dp_data, 3))

    def dp(self, flowrate: np.ndarray | float) -> np.ndarray | float:
        """
        Return the pressure rise over the pump as a function of the required flowrate.

        :param flowrate: The required flowrate, m^3/s
        :return: The total pressure rise, Pa.
        """

        # if (flowrate > self.flowMax) or (flowrate < self.flowMin):
        #     warnings.warn(f"Flowrate {flowrate:.5f} m^3/s lies outside of "
        #                   f"pump curve domain ({self.flowMin:.5f} to {self.flowMax:.5f})")

        return np.clip(self.poly(flowrate), 0, None)
