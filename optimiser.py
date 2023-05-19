import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

import main as ga3
from GA3_CONSTS import *


def num_passes():
    Npt = 1
    Nps = 1
    Nt = 13
    lt = np.min((LT_TOTAL / (Nt * Npt), LT_MAX_1P))  # tube lengths

    tube_arr = pd.read_csv("data/tube-arrangements.csv")
    tube_interp = interp1d(tube_arr["Nt"].to_numpy(), tube_arr["Y"].to_numpy())
    Y = (DS / 0.064) * tube_interp(Nt * Npt) / 1000  # tube pitch (reference diameter for data weas 0.064)

    Nb = 9

    isSquare = False
    G = 0.2 * DS

    Q = []

    hx = ga3.HX(COLDSTREAM, HOTSTREAM, KT, EPST, lt, DO, DI, Nt, Y, isSquare, Nps, Npt, Nb, G, DS, DN)
    hpump = ga3.Pump(ga3.Pump.HOT)
    cpump = ga3.Pump(ga3.Pump.COLD)

    mdot_h, dp_h = hx.chicSolver(hpump)
    mdot_c, dp_c = hx.chicSolver(cpump)

    hx.hydraulicAnalysisTube(mdot_h, True)
    hx.hydraulicAnalysisShell(mdot_c, True)

    print(mdot_h, mdot_c)
    print(hx.thermalAnalysis(mdot_h, mdot_c, True))


if __name__ == "__main__":
    num_passes()
