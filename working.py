import numpy as np
import pandas as pd

import main as ga3
from GA3_CONSTS import *

coldPump = ga3.Pump(ga3.Pump.COLD)
hotPump = ga3.Pump(ga3.Pump.HOT)

Nt = 13 # number of tubes
lt = 3.5 / Nt  # tube lengths
Y = 0.014 # tube pitch

Nb = 9  # Number of shell baffle
B = lt / (Nb + 1)  # Baffle pitch

isSquare = False
Nps = 1
Npt = 1

G = 0.2 * DS

hx = ga3.HX(COLDSTREAM, HOTSTREAM, KT, EPST, lt, DO, DI, Nt, Y, isSquare, Nps, Npt, Nb, G, DS, DN)

mdotc = ga3.chicSolver(hx, coldPump)[0] # Tube
mdoth = ga3.chicSolver(hx, hotPump)[0] # Shell

Q = hx.thermalAnalysis_LMTD(mdot_t = mdotc, mdot_s = mdoth, verbose = True)

#hx.plotHXChics()
#ga3.chicSolver()