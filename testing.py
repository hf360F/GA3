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

mdot_h, dp_h = ga3.chicSolver(HXtest, ga3.Pump(ga3.Pump.HOT))
mdot_c, dp_c = ga3.chicSolver(HXtest, ga3.Pump(ga3.Pump.COLD))
Q = HXtest.thermalAnalysis(mdot_h, mdot_c)

print(f'Q: {Q}\n'
      f'mdot_h: {mdot_h}\n'
      f'mdot_c: {mdot_c}\n'
      f'dp_h: {dp_h}\n'
      f'dp_c: {dp_c}\n')

