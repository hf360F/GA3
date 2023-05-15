import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

import main as ga3
from GA3_CONSTS import *


def num_passes():
    Nt = 13
    lt = np.clip(LT_TOTAL / Nt, None, LT_MAX_1P)  # tube lengths

    tube_arr = pd.read_csv("data/tube-arrangements.csv")
    tube_interp = interp1d(tube_arr["Nt"].to_numpy(), tube_arr["Y"].to_numpy())
    Y = (DS / 0.064) * tube_interp(Nt)  # tube pitch (reference diameter for data weas 0.064)

    Nb = np.arange(0, 10)

    isSquare = False
    Np = 1
    G = 0.2 * DS
    Q = []
    hx = ga3.HX(COLDSTREAM, HOTSTREAM, KT, EPST, lt, DO, DI, Nt, Y, isSquare, Np, Nb[0], G, DS, DN)
    hpump = ga3.Pump(ga3.Pump.HOT)
    cpump = ga3.Pump(ga3.Pump.COLD)

    for nb in Nb:
        hx.Nb = nb
        mdot_h, dp_h = hx.chicSolver(hpump)
        mdot_c, dp_c = hx.chicSolver(cpump)
        Q.append(hx.thermalAnalysis(mdot_h, mdot_c))

    Q = np.array(Q)

    plt.plot(Nb, Q)
    plt.ylim(0, 11000)

    plt.show()


if __name__ == "__main__":
    num_passes()
