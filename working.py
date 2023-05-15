import numpy as np
import pandas as pd

import main as ga3
from GA3_CONSTS import *

Nt = 13 # number of tubes
lt = 3.5 / Nt  # tube lengths
Y = 0.014 # tube pitch

Nb = 9  # Number of shell baffle
B = lt / (Nb + 1)  # Baffle pitch

isSquare = False
Nps = 1
Npt = 1

G = 0.2 * DS

hx = ga3.HX(COLDSTREAM, HOTSTREAM, KT, EPST, lt, DO, DI, Nt, Y, isSquare, Nps, Npt, Nb, B, G, DS, DN)

hx.hydraulicAnalysisTube(mdot=0.5, verbose=True)