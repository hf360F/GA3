import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

import main as ga3
from GA3_CONSTS import *


def num_passes():
    Npt = 2
    Nps = 1
    Nt = 12

    # lt = np.min((LT_TOTAL / Nt, LT_MAX))-2*END_WASTAGE
    lts = np.linspace(0.95 * 0.251, 1.05 * 0.251, 10)

    tube_arr = pd.read_csv("data/tube-arrangements.csv")
    tube_interp = interp1d(tube_arr["Nt"].to_numpy(), tube_arr["Y"].to_numpy())
    # Y = lambda Nt: (DS / 0.064) * tube_interp(Nt) / 1000  # tube pitch (reference diameter for data weas 0.064)
    Y = 0.014

    Nb = 12
    isSquare = False
    G = 0.2 * DS

    hx = ga3.HX(COLDSTREAM, HOTSTREAM, KT, EPST, lts[0], DO, DI, Nt, Y, isSquare, Nps, Npt, Nb, 0.0167, G, DS, DN)

    hpump_23 = ga3.Pump(ga3.Pump.HOT, 2023)
    cpump_23 = ga3.Pump(ga3.Pump.COLD, 2023)

    Q_NTU = np.zeros_like(lts)
    Q_LMTD = np.zeros_like(lts)

    for i, Nb in enumerate(lts):
        hx.lt
        hx.B = hx.lt - BAFFLE_END_SPACE / (Nb - 1)

        mdot_h, dp_h = hx.chicSolver(hpump_23)
        mdot_c, dp_c = hx.chicSolver(cpump_23)

        Q_NTU[i] = hx.thermalAnalysis_NTU(mdot_h, mdot_c)
        Q_LMTD[i] = hx.thermalAnalysis_LMTD(mdot_h, mdot_c)

    plt.plot(Nbs, Q_NTU, 'gx-', label='$\\epsilon$-NTU')
    plt.plot(Nbs, Q_LMTD, 'bx-', label='F-LMTD')
    plt.vlines(13, 0, 16.000, 'r')
    plt.annotate('$N_t=13$', (13.5, 1.000), color='red')
    plt.legend()
    plt.ylabel('$\dot Q$')
    plt.xlabel('$N_t$')
    plt.show()


if __name__ == "__main__":
    num_passes()
