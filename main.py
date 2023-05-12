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

epst = 0.0000015 # Pipe effective roughness height, mm

class HX:
    def__init__(self, coldStream, hotStream, kt, epst, lt, do, di, Nt, Y, isSquare, Np, Nb, B, G):
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
        self.Attot = self.Nt * self.di**2 * np.pi/4 # Tube total flowpath area, m^2
        self.Vt = self.mdt / (self.coldStream["rho"] * self.Attot) # Bulk tube velocity, m/s
        self.Ret = self.coldStream["rho"]*self.Vt*self.di/self.coldStream["mu"] # Tube Reynolds

        # Haaland approximation of Colebrook-White for Darcy friction factor 
        self.fTube = (-1.8*np.log10((self.epst/(self.di*3.7))**1.11 + (6.9/self.Ret)))**(-2)
        dpTube = self.fTube*(self.lt/self.di)*0.5*self.coldStream["rho"]*self.Vt**2

        # End analysis
        kc = 0.5
        ke = 0.8
        self.As = self.ds**2 * np.pi/4
        self.sigma = self.Attot/self.As

        # Implement kc, ke lookup based on Ret, sigma

        dpEnds = (kc + ke)*0.5*self.coldStream["rho"]*self.Vt**2

        # Nozzle analysis
        self.An = self.dn**2 * np.pi/4
        self.Vn = self.mdt / (self.coldStream["rho"] * self.An)

    def thermalAnalysis(self):
        pass