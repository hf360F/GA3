"""
EXTERNAL CONSTRAINTS:
Total mass < 1 kg
Total length < 0.35 m
Total copper tube length < 3.5 m

Compressor characteristics fixed (see csv)
Cold stream inlet fixed as 20 C
Hot stream inlet fixed as 60 C
"""

from typing import Union
from scipy.optimize import minimize_scalar

import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import warnings

def K(sigma, Ret):
    """Entrance and exit loss coefficient for a tube.  Assumes turbulent flow in tube.
       Will warn for out of range Reynolds number, but match to closest curve.

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

    if Ret > Remax:
        dfKc = pd.read_csv("data/Turb_10000_Kc.csv")
        dfKe = dfKe3 = pd.read_csv("data/Turb_10000_Ke.csv")
        warnings.warn(f"Tube Re = {Ret:.2f} out of range (Remax = {Remax:.2f}), matching to closest curve.")
    elif Ret > (Re2 + Remax)/2:
        dfKc = pd.read_csv("data/Turb_10000_Kc.csv")
        dfKe = dfKe3 = pd.read_csv("data/Turb_10000_Ke.csv")
    elif Ret > (Remin + Re2)/2:
        dfKc = pd.read_csv("data/Turb_5000_Kc.csv")
        dfKe = dfKe3 = pd.read_csv("data/Turb_5000_Ke.csv")
    elif Ret >= Remin:
        dfKc = pd.read_csv("data/Turb_3000_Kc.csv")
        dfKe = dfKe3 = pd.read_csv("data/Turb_3000_Ke.csv")
    else:
        dfKc = pd.read_csv("data/Turb_3000_Kc.csv")
        dfKe = dfKe3 = pd.read_csv("data/Turb_3000_Ke.csv")
        warnings.warn(f"Tube Re = {Ret:.2f} out of range (Remin = {Remin:.2f}), matching to closest curve.")
  
    # Fit polynomials to digitised data
    polyCoeffsKc = np.polyfit(dfKc["sigma"], dfKc["Kc"], polyOrder)
    polyCoeffsKe = np.polyfit(dfKe["sigma"], dfKe["Ke"], polyOrder)

    # Evaluate polynomials
    Kc, Ke = 0, 0
    for i in range(polyOrder+1):
        Kc += polyCoeffsKc[i] * sigma**(polyOrder-i)
        Ke += polyCoeffsKe[i] * sigma**(polyOrder-i)

    return Kc+Ke

def chicSolver(HX, pump):
    """Find intersection of pump and HX characteristic to set operating point.
    The HX flow path to use is inferred from pump.pump_type.
    
    Args:
        HX (HX object): Heat exchanger design.
        pump (Pump object): Pump design.

    Returns:
        (float, float): Mass flow rate, kg/s, and HX pressure drop (= pump pressure rise), Pa
    """  

    if pump.pump_type is pump.HOT:
        rho = HX.hotStream["rho"]
        # We have found the intersection of the curves when pressure change over pump and HX sum to zero:
        def f(mdot):
            return abs(HX.hydraulicAnalysisTube(mdot) - pump.dp(mdot/rho)) # Note Pump.dp takes volumetric flow rate
    else:
        rho = HX.coldStream["rho"]
        def f(mdot):
            return abs(HX.hydraulicAnalysisShell(mdot) - pump.dp(mdot/rho))
        
    solution = minimize_scalar(f, bounds=[0, 1])

    if solution["success"] is False:
        raise ValueError("Unable to intersect pump and heat exchanger characteristics!")
    else:
        return (solution["x"], pump.dp(solution["x"]/rho))

class HX:
    def __init__(self, coldStream, hotStream, kt, epst, lt, do, di, Nt, Y, isSquare, Np, Nb, B, G, ds, dn):
        """Heat exchanger design class.

        Args:
            coldStream (dict): Dictionary of cold inlet temperature 'Tci' C, specific heat capacity 'cp' J/kgK, fluid density 'rho' kg/m^3, thermal conductivity 'k' W/mK, viscocity 'mu' kg/ms.
            hotStream (dict): Dictionary of hot stream properties, as for cold stream.
            kt (float): Thermal conducitivty of tube, W/mK.
            epst (float): Effective roughness height of tube, m.
            lt (float): Length of tube section, m.
            do (float): Tube outer diameter, m.
            di (float): Tube inner diameter, m.
            Nt (int): Number of tubes.
            Y (float): Tube pitch, m.
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

    def hydraulicAnalysisTube(self, mdot, verbose=False):
        """Perform pressure drop analysis on tube flow path for given mdot.

        Args:
            mdot (float): Tube mass flow rate, kg/s. 
            verbose (bool): Whether or not to print intermediate values (velocities, areas, etc).       
        """

        # Tube friction
        self.Attot = self.Nt * self.di ** 2 * np.pi / 4  # Tube total flowpath area, m^2
        Vt = mdot / (self.hotStream["rho"] * self.Attot)  # Bulk tube velocity, m/s
        Ret = self.hotStream["rho"] * Vt * self.di / self.hotStream["mu"]  # Tube Reynolds

        # Haaland approximation of Colebrook-White for Darcy friction factor
        self.fTube = (-1.8 * np.log10((self.epst / (self.di * 3.7)) ** 1.11 + (6.9 / Ret))) ** (-2)
        dpFric = self.fTube * (self.lt / self.di) * 0.5 * self.hotStream["rho"] * (Vt**2)

        # End loss
        self.Apipe = self.ds ** 2 * np.pi / 4
        self.sigma = self.Attot * self.Np / self.Apipe # Note scaling with number of passes Np

        # Look up total loss factor, Kc+Ke, and calculate friction
        self.Ktot = K(self.sigma, Ret)
        dpEnds = self.Ktot * 0.5 * self.hotStream["rho"] * Vt ** 2

        # Nozzle loss
        self.An = self.dn ** 2 * np.pi / 4
        Vn = mdot / (self.hotStream["rho"] * self.An)
        dpNozzles = 2 * 0.5 * self.hotStream["rho"] * (Vn**2)

        tubeTotaldP = dpFric + dpEnds + dpNozzles

        if verbose:
            print(f"\nTUBE HYDRAULIC ANALYSIS SUMMARY FOR mdot = {mdot:.2f} kg/s\n")
            print(f"Tube flow area: {self.Attot:.6f} m^2")
            print(f"Tube bulk velocity: {Vt:.2f} m/s")
            print(f"Tube Reynolds number: {Ret:.0f}")
            print(f"Tube friction factor: {self.fTube:.4f}")
            print(f"Total inlet/exit loss factor: {self.Ktot:.3f}")
            print(f"Total pressure drop: {tubeTotaldP:.0f} Pa (friction {dpFric:.0f},"\
                  f" ends {dpEnds:.0f}, nozzles {dpNozzles:.0f})\n")
            
        return tubeTotaldP 

    def hydraulicAnalysisShell(self, mdot, verbose=False):
        """Perform pressure drop analysis on shell flow path for given mdot.

        Args:
            mdot (float): Shell mass flow rate, kg/s. 
            verbose (bool): Whether or not to print intermediate values (velocities, areas, etc).       
        """
        
        # Shell loss
        self.As = self.ds*(self.Y - self.do)*self.B/self.Y # Approximation of shell flow area
        Vs = mdot/(self.coldStream["rho"]*self.As)

        self.dseff = self.ds*(self.As/self.Apipe)
        Res = self.coldStream["rho"]*Vs*self.dseff/self.coldStream["mu"]

        if self.isSquare:
            a = 0.34
        else:
            a = 0.2
        
        shelldp1 = 4*a*(Res**(-0.15))*self.Nt*self.coldStream["rho"]*(Vs**2)

        # Nozzle loss
        self.An = self.dn ** 2 * np.pi / 4
        Vn = mdot / (self.hotStream["rho"] * self.An)
        dpNozzles = 2 * 0.5 * self.coldStream["rho"] * (Vn**2)

        shellTotaldp = shelldp1 + dpNozzles

        if verbose:
            print(f"\nSHELL HYDRAULIC ANALYSIS SUMMARY FOR mdot = {mdot:.2f} kg/s\n")
            print(f"Shell flow area: {self.As:.6f} m^2")
            print(f"Shell characteristic velocity: {Vs:.2f} m/s")
            print(f"Shell Reynolds number: {Res:.0f}")
            print(f"Total pressure drop: {shellTotaldp:.0f} Pa (shell {shelldp1:.0f},"\
                  f"  nozzles {dpNozzles:.0f})\n")

        return shellTotaldp

    def plotHXChics(self, mdotMin=0, mdotMax=1, n=100):
        """Plot HX flow characteristics.

        Args:
            mdotMin (float): Minimum mass flow to evaluate pressure drop at, kg/s.
            mdotMax (float): Maximum mass flow to evaluate pressure drop at, kg/s.
            n (int): Number of mass flow rates to evaluate between mdotMin and mdotMax.
        """

        mdots = np.linspace(mdotMin, mdotMax, n+1)
        dpsTube, dpsShell = np.zeros_like(mdots), np.zeros_like(mdots)

        for i in range(n+1):
            dpsTube[i] = self.hydraulicAnalysisTube(mdot=mdots[i])
            dpsShell[i] = self.hydraulicAnalysisShell(mdot=mdots[i])
        
        plt.plot(mdots, dpsTube/1000, label="Tube flow path", color="red")
        plt.plot(mdots, dpsShell/1000, label="Shell flow path", color="blue")
        plt.xlim(mdotMin, mdotMax)
        plt.ylim(0)
        plt.xlabel("Mass flow rate, $kg/s$")
        plt.ylabel("Pressure drop, $kPa$")
        plt.title("Heat exchanger characteristics")
        plt.legend()
        plt.grid()
        plt.show()

    def thermalAnalysis(self):
        pass


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

        if (flowrate > self.flowMax) or (flowrate < self.flowMin):
            raise ValueError(f"Flowrate {flowrate:.5f} m^3/s lies outside of "\
                             f"pump curve domain ({self.flowMin:.5f} to {self.flowMax:.5f})")

        return np.clip(self.poly(flowrate), 0, None)
