import main as ga3
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

# INLET FLUID PROPERTIES
coldStream = {"Tci": 20,  # C
              "cp": 4179,  # J/kgK
              "rho": 990.1,  # kg/m^3
              "k": 0.632,  # W/mK
              "mu": 6.51E-4}  # kg/ms

hotStream = {"Tci": 60,  # C
             "cp": 4179,  # J/kgK
             "rho": 990.1,  # kg/m^3
             "k": 0.632,  # W/mK
             "mu": 6.51E-4}  # kg/ms

# TUBE PROPERTIES
epst = 0.0000015  # Effective roughness height for new drawn copper, mm
kt = 386 # Thermal conductivity of copper, W/mK
do, di = 0.008, 0.006 # Tube outer and inner diameter, m

# SHELL PROPERTIES
ds = 0.064 # Shell inner diameter, m

# NOZZLE PROPERTIES
dn = 0.02 # Nozzle diameter, m

# FREE VARIABLES
lt = 0.35 # Tube length, m
Nt = 13 # Tube count
Y = 0.014 # Tube pitch, m
isSquare = True # Square or triangular tube pattern

Np = 1 # Number of passes
Nb = 9 # Number of shell baffle
B = lt/(Nb + 1) # Baffle pitch

G = 0.2*ds # Baffle cut, m NOT YET USED

# Test heat exchanger using worked example values
HXtest = ga3.HX(coldStream, hotStream, kt, epst, lt, do, di, Nt, Y, isSquare, Np, Nb, B, G, ds, dn)
HXtest.hydraulicAnalysisTube(mdot=0.45, verbose=True) # Worked example tube pressure drop 4360 Pa
HXtest.hydraulicAnalysisShell(mdot=0.50, verbose=True) # Worked example shell pressure drop 3700 Pa

#HXtest.plotHXChics()
#print(ga3.chicSolver(HXtest, ga3.Pump(ga3.Pump.HOT)))
#print(ga3.chicSolver(HXtest, ga3.Pump(ga3.Pump.COLD)))