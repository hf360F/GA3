"""
EXTERNAL CONSTRAINTS:
Total mass < 1 kg
Total length < 0.35 m
Total copper tube length < 3.5 m

Compressor characteristics fixed (see csv)
Cold stream inlet fixed as 20 C
Hot stream inlet fixed as 60 C
"""

import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Inlet properties
coldStream = {"Tci": 20, # C
              "cp": 4179, # J/kgK
              "rho": 990.1, # kg/m^3
              "k": 0.632, # W/mK
              "mu": 6.51E-4 # kg/ms
              }

hotStream =  {"Tci": 60, # C
              "cp": 4179, # J/kgK
              "rho": 990.1, # kg/m^3
              "k": 0.632, # W/mK
              "mu": 6.51E-4 # kg/ms
              }

class HX:
    def__init__(self, coldStream, hotStream, kt, lt, do, di, Nt, Y, isSquare, Np, Nb, B, G):
    """Heat exchanger design class.
    
    Args:
        coldStream (dict): Dictionary of cold inlet temperature 'Tci' C, specific
                           heat capacity 'cp' J/kgK, fluid density 'rho' kg/m^3,
                           thermal conductivity 'k' W/mK, viscocity 'mu' kg/ms.
        hotStream (dict): Dictionary of hot stream properties, as for cold stream.
        
        kt (float): Thermal conducitivty of tube, W/mK.
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

        dn (float): Nozzle diameter for both streams, m.

    """

    self.coldStream = coldStream
    self.hotStream = hotStream
    self.kt = kt
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

    def hydraulicAnalysisTube(self, mdt):
        """Perform pressure drop analysis on tube flow path for given mdot.

        Args:
            mdt (float): Tube mass flow rate, kg/s.        
        """

        # Tube analysis
        self.Attot = self.Nt * self.di**2 * np.pi/4 # Tube total flowpath area, m^2
        self.Vt = self.mdt / (self.coldStream["rho"] * self.Attot) # Bulk tube velocity, m/s
        self.Ret = self.coldStream["rho"]*self.Vt*self.di/self.coldStream["mu"] # Tube Reynolds

        # Nozzle analysis
        self.Vn


    def thermalAnalysis(self):