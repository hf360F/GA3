"""
EXTERNAL CONSTRAINTS:
Total mass < 1 kg
Total length < 0.35 m
Total copper tube length < 3.5 m

Compressor characteristics fixed (see csv)
Cold stream inlet fixed as 20 C
Hot stream inlet fixed as 60 C
"""

import fluids.fittings as ft
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import root, root_scalar

from GA3_CONSTS import *


def Flookup(P, R, Nps):
    """Temperature delta correction factor, for Q = HAF * LMTD.

    Args:
        P (float): Ratio of cold stream inlet-outlet dT to hot to cold stream inlet dT.
        R (float): Heat capacity ratio, cold stream to hot stream.
        Nps (int): Number of shell passes, data only covers Nps = 1 or Nps = 2.

    Returns:
        (float): Correction factor F.
    """
    if Nps == 1:
        df = pd.read_csv("data/F1pass.csv").to_numpy()
    elif Nps == 2:
        df = pd.read_csv("data/F2pass.csv").to_numpy()
    else:
        raise (ValueError("Data for temperature delta correction only supports 1 or 2 shell passes"))

    R = np.clip(R, 0.4, 2)
    P = np.clip(P, 0, 1)

    F = scipy.interpolate.griddata(df[:, :-1], df[:, -1], np.column_stack((R, P)), "cubic")

    return np.clip(F, 0, 1)[0]


class HX:
    def __init__(self, coldStream, hotStream, kt, epst, lt, do, di, Nt, Y, isSquare, Nps, Npt, Nb, B, G, ds, dn):
        """Heat exchanger design class.

        Args:
            coldStream (dict): Dictionary of cold inlet temperature 'Ti' C, specific heat capacity 'cp' J/kgK,
                fluid density 'rho' kg/m^3, thermal conductivity 'k' W/mK, viscocity 'mu' kg/ms.
            hotStream (dict): Dictionary of hot stream properties, as for cold stream.
            kt (float): Thermal conducitivty of tube, W/mK.
            epst (float): Effective roughness height of tube, m.
            lt (float): Length of tube section, m.
            do (float): Tube outer diameter, m.
            di (float): Tube inner diameter, m.
            Nt (int): Number of tubes.
            Y (float): Tube pitch, m.
            isSquare (bool): Are tubes arranged in square pattern, if False tubes are asssumed to be in triangular pattern.
            Nps (int): Number of shell passes.
            Npt (int): Number of tube passes.
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
        self.Nps = Nps
        self.Npt = Npt
        self.Nb = Nb
        self.B = B
        self.G = G
        self.ds = ds
        self.dn = dn

        self.dfKc = pd.read_csv("data/Turb_Kc.csv").to_numpy()
        self.dfKe = pd.read_csv("data/Turb_Ke.csv").to_numpy()

        if Nps != 1:
            raise (NotImplementedError("F factor for multi-pass setups not implemented yet"))

    @property
    def Attot(self):
        return self.Nt / self.Npt * self.di ** 2 * np.pi / 4  # Tube total flowpath area, m^2

    @property
    def Apipe(self):
        return self.ds ** 2 * np.pi / 4

    @property
    def sigma(self):
        return self.Npt * self.Attot / self.Apipe

    @property
    def An(self):
        return self.dn ** 2 * np.pi / 4

    @property
    def As(self):
        return self.ds * self.B * (1 - self.do / self.Y) / self.Nps  # Approximation of shell flow area

    @property
    def dseff(self):
        # hydraulic diameter for rectangular channel
        a = self.B
        b = (self.Y - self.do)
        # return self.Nt ** 0.5 * (2 * a * b / (a + b))
        return self.ds * (self.As / self.Apipe)

    @property
    def mass(self):
        tube = M_TUBE * self.Nt * self.lt
        nozzles = M_NOZZLE * 4
        baffles = M_BAFFLE * self.Nb * (self.Apipe - 0.25 * self.Nt * np.pi * self.do ** 2)
        shell = M_SHELL * ((self.lt + 2 * END_WASTAGE) + 0.04 * 2)
        ends = M_ENDS * ((2 * 0.00635 * self.Apipe + 2 * 0.008 * (1 - self.sigma)) * self.Apipe)
        orings = M_ORING_L * 4 + M_ORING_S * self.Nt

        return tube + nozzles + baffles + shell + ends + orings

    def ft_darcy(self, Ret):
        # Haaland's approximation of Colebrook-White for Darcy friction factor
        return (-1.8 * np.log10((self.epst / (self.di * 3.7)) ** 1.11 + (6.9 / Ret))) ** (-2)

    def K(self, sigma, Ret):
        Ret = np.clip(Ret, 3000, 10000)
        sigma = np.clip(sigma, 0, 1)

        Kc = scipy.interpolate.griddata(self.dfKc[:, :-1], self.dfKc[:, -1], (Ret, sigma), 'cubic')
        Ke = scipy.interpolate.griddata(self.dfKe[:, :-1], self.dfKe[:, -1], (Ret, sigma), 'cubic')

        return Kc + Ke

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

    def hydraulicAnalysisTube(self, mdot, verbose=False):
        """Perform pressure drop analysis on tube flow path for given mdot.

        Args:
            mdot (float): Tube mass flow rate, kg/s. 
            verbose (bool): Whether or not to print intermediate values (velocities, areas, etc).       
        """

        Vt = mdot / (self.hotStream["rho"] * self.Attot)  # Bulk tube velocity, m/s
        Ret = self.hotStream["rho"] * Vt * self.di / self.hotStream["mu"]  # Tube Reynolds

        # Haaland's approximation of Colebrook-White for Darcy friction factor
        dpFric = self.ft_darcy(Ret) * (self.lt * self.Npt / self.di) * 0.5 * self.hotStream["rho"] * (Vt ** 2)

        # Sum entrance/exit loss factor, 180 degree bend loss factor (for Npt > 1), to calculate total minor loss
        Kreturn = ft.bend_rounded(Di=self.di, angle=180, fd=self.ft_darcy(Ret), rc=self.Y / 2, Re=Ret, method="Rennels")
        Kbaffle = 0.7
        Ktot = self.K(self.sigma, Ret) + (self.Npt - 1) * Kreturn + self.Nb * Kbaffle

        dpMinor = Ktot * 0.5 * self.hotStream["rho"] * Vt ** 2

        # Nozzle loss
        Vn = mdot / (self.hotStream["rho"] * self.An)
        dpNozzles = 0.5 * self.hotStream["rho"] * (Vn ** 2)

        tubeTotaldP = dpFric + dpMinor + dpNozzles

        if verbose:
            print(f"\nTUBE HYDRAULIC ANALYSIS SUMMARY FOR mdot = {mdot:.2f} kg/s\n")
            print(f"Tube flow area: {self.Attot:.6f} m^2")
            print(f"Tube bulk velocity: {Vt:.2f} m/s")
            print(f"Tube Reynolds number: {Ret:.0f}")
            print(f"Tube friction factor: {self.ft_darcy(Ret):.6f}")
            print(f"Sum of minor loss coefficients: {Ktot:.3f}")
            print(f"Total pressure drop: {tubeTotaldP:.0f} Pa (friction {dpFric:.0f},"
                  f" minor losses {dpMinor:.0f}, nozzles {dpNozzles:.0f})\n")

        return tubeTotaldP

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
        dpNozzles = 0.5 * self.coldStream["rho"] * (Vn ** 2)

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

    def HA(self, mdot_t, mdot_s):
        # Reynolds numbers
        Ret = mdot_t * self.di / self.Attot / self.hotStream["mu"]
        Res = mdot_s * self.dseff / self.As / self.coldStream["mu"]

        # Pr and geometry constant are const
        Pr = 4.31
        c = 0.15 if self.isSquare else 0.2

        # Nusselt number correlations

        # Tube-side
        # Dittus-Boelter
        # Nut = 0.026 * Ret ** 0.8 * Pr ** 0.3

        # Petukhov-Popov
        a = 1.07 + 900 / Ret + 0.63 / (1 + 10 * Pr)
        Cf = self.ft_darcy(Ret) / 4
        Nut = (Cf / 2) * Ret * Pr / (a + 12.7 * np.sqrt(Cf / 2) * (Pr ** (2 / 3) - 1))

        # Gnielinski
        # Nut = (Cf/2)*(Ret-1000)*Pr/(1+12.7*np.sqrt(Cf/2)*(Pr**(2/3)-1))

        # Shell-side needs work
        Nus = c * Res ** 0.6 * Pr ** 0.3

        # convection ht coeffs
        hi = Nut * self.hotStream["k"] / self.di
        ho = Nus * self.coldStream["k"] / self.do

        # total ht coeff
        H = 1 / (1 / hi + self.di * np.log(self.do / self.di) / 2 / self.kt + self.di / self.do / ho)

        # SFEE
        HA = H * np.pi * self.di * self.lt * self.Nt
        return HA

    def thermalAnalysis_NTU(self, mdot_t: float, mdot_s: float, verbose: bool = False) -> float:
        C1 = mdot_s * COLDSTREAM["cp"]
        C2 = mdot_t * HOTSTREAM["cp"]
        Cmin = np.min((C1, C2))
        Cmax = np.max((C1, C2))
        NTU = self.HA(mdot_t, mdot_s) / Cmin
        Cstar = Cmin / Cmax

        if self.Nb < 6:
            return None
            # raise ValueError(f"No correlation for Nb<6: Nb={self.Nb}")

        match self.Npt:
            case 1:
                # pure counterflow
                epsilon = (1 - np.exp(-NTU * (1 - Cstar))) / (1 - Cstar * np.exp(-NTU * (1 - Cstar)))
            case 2:
                # TEMA-E_1,2 HX
                Gamma = np.sqrt(1 + Cstar ** 2)
                epsilon = 2 / ((1 + Cstar) + Gamma + 1 / np.tanh(NTU * Gamma / 2))
            case x:
                return None
                # raise ValueError(f"No correlation defined for Npt>2: Npt={x}")

        Q = epsilon * Cmin * (HOTSTREAM["Ti"] - COLDSTREAM["Ti"])

        Tco = COLDSTREAM["Ti"] + Q / C1
        Tho = HOTSTREAM["Ti"] - Q / C2

        if verbose:
            print(f"\nHX THERMAL ANALYSIS SUMMARY - NTU\n\n"

                  f"  Heat transfer rate Q = {Q / 1000:.2f} kW.\n"
                  f"  NTU = {NTU:.3f}, effectiveness = {epsilon:.3f}.\n"
                  f"  mdoth = {mdot_t:.3f} kg/s, Tco = {Tco:.2f} C, deltaTc = {Q / C1:.2f} C.\n"
                  f"  mdotc = {mdot_s:.3f} kg/s, Tho = {Tho:.2f} C, deltaTh = {Q / C2:.2f} C.\n")
        return Q

    def thermalAnalysis_LMTD(self, mdot_t: float, mdot_s: float, verbose: bool = False, return_HAF=False) -> float:
        """
        Perform thermal analysis on the HX given the 2 mass flowrates.

        Args:
            mdot_t (float): Tube mass flow, kg/s.
            mdot_s (float): Shell mass flow, kg/s.
            verbose (bool): Indicate verbosity of output.

        Returns:
            Q (float): Resultant heat transfer rate, W.
        """
        HA = self.HA(mdot_t, mdot_s)

        # inlet temps
        Thi = self.hotStream["Ti"]
        Tci = self.coldStream["Ti"]

        def LMTD(Tho_, Tco_):
            T1 = Thi - Tco_
            T2 = Tho_ - Tci
            return T1 if np.isclose(T1, T2) else (T1 - T2) / np.log(T1 / T2)

        def f(x):
            """
            Returns error between Q calculated by LMTD.H.A.F and m.cp.dT.

            Args
                To (float, float): Hot and cold stream outlet temperatures, C.

            Returns:
                (float, float): Vector of dQs when found by LMTD H A F = (mdotcpdT)_hot and LMTD H A F = (mdotcpdT)_cold.

            """

            Tho, Tco = x
            F = Flookup(P=(Tco - Tci) / (Thi - Tci), R=(mdot_t / mdot_s), Nps=self.Nps)  # Update next F based on Q
            Q = HA * LMTD(Tho, Tco) * F
            return mdot_t * self.hotStream["cp"] * (Thi - Tho) - Q, mdot_s * self.coldStream["cp"] * (Tco - Tci) - Q

        # initial guess of outlet temperatures
        x0 = np.array((55, 25))
        # least squares solver sets errors to zero, optimising for outlet temperatures
        res = root(f, x0)
        if not res.success:
            raise AssertionError(f"Could not solve for outlet temperatures and heat transfer.")
        # recover vectors of outlet temperatures
        Tho, Tco = res.x
        F = Flookup(P=(Tco - Tci) / (Thi - Tci), R=(mdot_t / mdot_s), Nps=self.Nps)
        Q = HA * LMTD(Tho, Tco) * F

        if verbose:
            print(f"\nHX THERMAL ANALYSIS SUMMARY - LMTD\n\n"

                  f"  Heat transfer rate Q = {Q / 1000:.2f} kW.\n"
                  f"  temperature delta correction F = {F:.3f}.\n"
                  f"  mdoth = {mdot_t:.3f} kg/s, Tco = {Tco:.2f} C, deltaTc = {(Tco - self.coldStream['Ti']):.2f} C.\n"
                  f"  mdotc = {mdot_s:.3f} kg/s, Tho = {Tho:.2f} C, deltaTh = {(self.hotStream['Ti'] - Tho):.2f} C.\n")

        if not return_HAF:
            return Q
        else:
            return HA * F


class Pump:
    COLD, HOT = range(2)

    def __init__(self, pump_type: int, year: int) -> None:
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

        year = 2020 if year == 2019 else year
        assert year in (2017, 2018, 2020, 2022, 2023)

        data = np.genfromtxt(f"data/{type_str}-{year}.csv", delimiter=',')
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
        c = 0.5

        return np.clip(self.poly(flowrate), 0, None) + c * 0.5 * 990 * (flowrate / (np.pi * 0.01 ** 2)) ** 2
