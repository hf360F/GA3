import main as ga3
from GA3_CONSTS import *

# FREE VARIABLES
lt = 0.35  # Tube length, m
Nt = 13  # Tube count
Y = 0.014  # Tube pitch, m
isSquare = True  # Square or triangular tube pattern

Np = 1  # Number of passes
Nb = 9  # Number of shell baffle
B = lt / (Nb + 1)  # Baffle pitch

G = 0.2 * DS  # Baffle cut, m NOT YET USED

# Test heat exchanger using worked example values
HXtest = ga3.HX(COLDSTREAM, HOTSTREAM, KT, EPST, lt, DO, DI, Nt, Y, isSquare, Np, Nb, B, G, DS, DN)
HXtest.hydraulicAnalysisTube(mdot=0.45, verbose=True)  # Worked example tube pressure drop 4360 Pa
HXtest.hydraulicAnalysisShell(mdot=0.50, verbose=True)  # Worked example shell pressure drop 3700 Pa
HXtest.thermalAnalysis(mdot_t=0.45, mdot_s=0.50)

# HXtest.plotHXChics()
# print(ga3.chicSolver(HXtest, ga3.Pump(ga3.Pump.HOT)))
# print(ga3.chicSolver(HXtest, ga3.Pump(ga3.Pump.COLD)))
