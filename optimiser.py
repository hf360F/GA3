import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import main as ga3
from GA3_CONSTS import *


def num_tubes():
    # define solution space
    Nt = np.arange(13, 48)  # number of tubes
    lt = 0.27 / Nt  # tube lengths
    tube_arr = pd.read_csv("data/tube-arrangements.csv")
    Y = (DS / 0.064) * np.interp(Nt, tube_arr["Nt"].to_numpy(), tube_arr["Y"].to_numpy())  # tube pitch

    Nb = 9  # Number of shell baffle
    B = lt / (Nb + 1)  # Baffle pitch

    isSquare = False
    Np = 1

    G = 0.2 * DS

    hx = ga3.HX(COLDSTREAM, HOTSTREAM, KT, EPST, lt, DO, DI, Nt, Y, isSquare, Np, Nb, B, G, DS, DN)
    mdot_h, dp_h = ga3.chicSolver(hx, ga3.Pump(ga3.Pump.HOT))
    mdot_c, dp_c = ga3.chicSolver(hx, ga3.Pump(ga3.Pump.COLD))
    Q = hx.thermalAnalysis(mdot_h, mdot_c)
    plt.plot(Nt, Q)
    plt.show()


if __name__ == "__main__":
    num_tubes()
