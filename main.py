"""
EXTERNAL CONSTRAINTS:
Total mass < 1 kg
Total length < 0.35 m
Total copper tube length < 3.5 m

Compressor characteristics fixed (see csv)
Cold stream inlet fixed as 20 C
Hot stream inlet fixed as 60 C
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import root, root_scalar


class HX:
    def __init__(self, coldStream, hotStream, kt, epst, lt, do, di, Nt, Y, isSquare, Np, Nb, G, ds, dn):
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
            Nt (int): Number of tubes.
            Y (float): Tube pitch, m.
            isSquare (bool): Are tubes arranged in square pattern, if False tubes are asssumed to be in triangular pattern.
            Np (int): Number of passes.
            Nb (int): Number of baffles.
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
        self.G = G
        self.ds = ds
        self.dn = dn
        self.F = 0.9

        self.dfKc = pd.read_csv("data/Turb_Kc.csv").to_numpy()
        self.dfKe = pd.read_csv("data/Turb_Ke.csv").to_numpy()

        if Np != 1:
            raise (NotImplementedError("F factor for multi-pass setups not implemented yet"))

    @property
    def Attot(self):
        return self.Nt * self.di ** 2 * np.pi / 4  # Tube total flowpath area, m^2

    @property
    def Apipe(self):
        return self.ds ** 2 * np.pi / 4

    @property
    def sigma(self):
        return self.Attot * self.Np / self.Apipe  # Note scaling with number of passes Np

    @property
    def An(self):
        return self.dn ** 2 * np.pi / 4

    @property
    def B(self):
        return self.lt / (self.Nb + 1)

    @property
    def As(self):
        return self.ds * (self.Y - self.do) * self.B / self.Y  # Approximation of shell flow area

    @property
    def dseff(self):
        return self.ds * (self.As / self.Apipe)

    def K(self, sigma, Ret):
        Ret = np.clip(Ret, 3000, 10000)
        sigma = np.clip(sigma, 0, 1)

        Kc = scipy.interpolate.griddata(self.dfKc[:, :-1], self.dfKc[:, -1], (Ret, sigma), 'cubic')
        Ke = scipy.interpolate.griddata(self.dfKe[:, :-1], self.dfKe[:, -1], (Ret, sigma), 'cubic')

        return Kc + Ke

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
        Ktot = self.K(self.sigma, Ret)
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

    def chicSolver(self, pump):
        """Find intersection of pump and HX characteristic to set operating point.
        The HX flow path to use is inferred from pump.pump_type.

        Args:
            pump (Pump object): Pump design.

        Returns:
            (float, float): Mass flow rate, kg/s, and HX pressure drop (= pump pressure rise), Pa
        """

        if pump.pump_type is Pump.HOT:
            rho = self.hotStream["rho"]
            hx_dp = self.hydraulicAnalysisTube

            # We have found the intersection of the curves when pressure change over pump and HX sum to zero:
        else:
            rho = self.coldStream["rho"]
            hx_dp = self.hydraulicAnalysisShell

        def f(mdot):
            return hx_dp(mdot) - pump.dp(mdot / rho)

        solution = root_scalar(f, bracket=[0.001, 1])

        if not solution.converged:
            raise ValueError("Unable to intersect pump and heat exchanger characteristics!")
        else:
            return solution.root, pump.dp(solution.root / rho)

    def hydraulicAnalysisShell(self, mdot, verbose=False):
        """Perform pressure drop analysis on shell flow path for given mdot.

        Args:
            mdot (float): Shell mass flow rate, kg/s. 
            verbose (bool): Whether or not to print intermediate values (velocities, areas, etc).       
        """

        # Shell loss
        Vs = mdot / (self.coldStream["rho"] * self.As)

        Res = self.coldStream["rho"] * Vs * self.dseff / self.coldStream["mu"]

        a = 0.34 if self.isSquare else 0.2

        shelldp1 = 4 * a * (Res ** (-0.15)) * self.Nt * self.coldStream["rho"] * (Vs ** 2)

        # Nozzle loss
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

    def thermalAnalysis(self, mdot_t: float, mdot_s: float) -> float:
        """
        Perform thermal analysis on the HX given the 2 mass flowrates.

        Args:
            mdot_t (float): Tube mass flow, kg/s.
            mdot_s (float): Shell mass flow, kg/s

        Returns:
            Q (float): Resultant heat transfer rate, W.
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
        Thi = self.hotStream["Ti"]
        Tci = self.coldStream["Ti"]

        # SFEE
        HA = H * np.pi * self.di * self.lt * self.Nt

        def LMTD(Tho_, Tco_):
            T1 = Thi - Tco_
            T2 = Tho_ - Tci
            return T1 if np.isclose(T1, T2) else (T1 - T2) / np.log(T1 / T2)

        # error function for LMTD solver
        def f(To: np.ndarray):
            """
            Returns error between Q calculated by LMTD.H.A.F and m.cp.dT

            Args:
                To (np.ndarray): [Tho, Tco]

            Returns:
                dQ (ndarray): [dQ1,dQ2]

            """
            Tho, Tco = To
            Q = HA * LMTD(Tho, Tco) * self.F
            return np.array(
                (mdot_t * self.hotStream["cp"] * (Thi - Tho) - Q, mdot_s * self.coldStream["cp"] * (Tco - Tci) - Q))

        # initial guess of outlet temperatures
        x0 = np.array([40, 30])  # least squares solver sets errors to zero, optimising for outlet temperatures
        res = root(f, x0)
        if not res.success:
            raise AssertionError(f"Could not solve for outlet temperatures and heat transfer")
        # recover vectors of outlet temperatures
        Tho, Tco = res.x

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
