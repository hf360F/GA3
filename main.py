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

import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import warnings

# Inlet properties
coldStream = {"Tci": 20,  # C
              "cp": 4179,  # J/kgK
              "rho": 990.1,  # kg/m^3
              "k": 0.632,  # W/mK
              "mu": 6.51E-4  # kg/ms
              }

hotStream = {"Tci": 60,  # C
             "cp": 4179,  # J/kgK
             "rho": 990.1,  # kg/m^3
             "k": 0.632,  # W/mK
             "mu": 6.51E-4  # kg/ms
             }

def K(sigma, Ret, verbose = False):
    """Entrance and exit loss coefficient for a tube.  Assumes turbulent flow in tube.
       Will warn for out of range Reynolds number, but match to closest curve.

    Args:
        sigma (float): Ratio of total tube area to shell area.
        Ret (float): Tube Reynolds number.

    Returns:
        float: Sum of Kc and Ke loss coefficients.
    """

    polyOrder = 2

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


class HX:
    def __init__(self, coldStream, hotStream, kt, epst, lt, do, di, Nt, Y, isSquare, Np, Nb, B, G, ds, dn):
        """Heat exchanger design class.

        Args:
            coldStream (dict): Dictionary of cold inlet temperature 'Tci' C, specific
                               heat capacity 'cp' J/kgK, fluid density 'rho' kg/m^3,
                               thermal conductivity 'k' W/mK, viscocity 'mu' kg/ms.
            hotStream (dict): Dictionary of hot stream properties, as for cold stream.

            kt (float): Thermal conducitivty of tube, W/mK.
            epst (float): Effective roughness height of tube, m.
            lt (float): Length of tube section, m.
            do (float): Tube outer diameter, m
            di (float): Tube inner diameter, m
            Nt (int): Number of tubes.
            Y (float): Tube pitch, m.
            isSquare (bool): Are tubes arranged in square pattern, if False tubes are
                             asssumed to be in triangular pattern.
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

    def hydraulicAnalysisTube(self, mdt):
        """Perform pressure drop analysis on tube flow path for given mdot.

        Args:
            mdt (float): Tube mass flow rate, kg/s.        
        """

        # Tube analysis
        self.Attot = self.Nt * self.di ** 2 * np.pi / 4  # Tube total flowpath area, m^2
        self.Vt = self.mdt / (self.coldStream["rho"] * self.Attot)  # Bulk tube velocity, m/s
        self.Ret = self.coldStream["rho"] * self.Vt * self.di / self.coldStream["mu"]  # Tube Reynolds

        # Haaland approximation of Colebrook-White for Darcy friction factor
        self.fTube = (-1.8 * np.log10((self.epst / (self.di * 3.7)) ** 1.11 + (6.9 / self.Ret))) ** (-2)
        dpTube = self.fTube * (self.lt / self.di) * 0.5 * self.coldStream["rho"] * self.Vt ** 2

        # End analysis
        self.As = self.ds ** 2 * np.pi / 4
        self.sigma = self.Attot * self.Np / self.As

        # Implement kc, ke lookup based on Ret, sigma

        dpEnds = K(self.sigma, self.Ret) * 0.5 * self.coldStream["rho"] * self.Vt ** 2

        # Nozzle analysis
        self.An = self.dn ** 2 * np.pi / 4
        self.Vn = self.mdt / (self.coldStream["rho"] * self.An)

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
        self.dp_data = data[:, 1]
        self.poly = np.poly1d(np.polyfit(self.flowrate_data, self.dp_data, 2))

    def dp(self, flowrate: np.ndarray | float) -> np.ndarray | float:
        """
        Return the pressure rise over the pump as a function of the required flowrate.

        :param flowrate: The required flowrate, m^3/s
        :return: The total pressure rise, Pa.
        """

        return np.clip(self.poly(flowrate), 0, None)
