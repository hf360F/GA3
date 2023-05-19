import numpy as np
from scipy.interpolate import interp1d
import pandas as pd

import main as ga3
from GA3_CONSTS import *

Npop = 1000

# Free variables
# lt, Nt, Y, isSquare, Nps, Npt, Nb

isSquare = True

NtValues = np.arange(1, 48, 1)
NpsValues = np.arange(1, 11, 1)
NptValues = np.arange(1, 11, 1)
Nb = np.arange(1, 31, 1)


tube_arr = pd.read_csv("data/tube-arrangements.csv")
tube_interp = interp1d(tube_arr["Nt"].to_numpy(), tube_arr["Y"].to_numpy())

def randomHX():
    Nt = np.random.choice(NtValues)
    Nps = np.random.choice(NpsValues)
    Npt = np.random.choice(NptValues)

    Ymax = (DS / 0.064) * tube_interp(Nt) # Maximum Y for given Nt
    Y = np.random.choice(np.linspace(0, Ymax, 100))

    return ga3.HX(COLDSTREAM, HOTSTREAM, KT, EPST, lt, DO, DI, Nt, Y, isSquare, Nps, Npt, Nb, G, DS, DN)