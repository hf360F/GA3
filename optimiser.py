import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

import main as ga3
from GA3_CONSTS import *

def num_tubes():
    # define solution space
    Nt = np.arange(7, 20)  # number of tubes
    lt = np.clip(LT_TOTAL / Nt, None, LT_MAX_1P) - 2 * END_WASTAGE  # tube lengths
    tube_arr = pd.read_csv("data/tube-arrangements.csv")
    tube_interp = interp1d(tube_arr["Nt"].to_numpy(), tube_arr["Y"].to_numpy())
    Y = (DS / 0.064) * tube_interp(Nt)  # tube pitch (reference diameter for data weas 0.064)

    Nb = 9  # Number of shell baffle

    isSquare = False
    Nps = 1
    Npt = 1

    G = 0.2 * DS

<<<<<<< HEAD
    hx = ga3.HX(COLDSTREAM, HOTSTREAM, KT, EPST, lt, DO, DI, Nt, Y, isSquare, Nps, Npt, Nb, B, G, DS, DN)
=======
    hx = ga3.HX(COLDSTREAM, HOTSTREAM, KT, EPST, lt, DO, DI, Nt, Y, isSquare, Np, Nb, G, DS, DN)
>>>>>>> dbd3188fe97cbf808e6c9e137e26b1a7a7b5d1d5
    mdot_h, dp_h = ga3.chicSolver(hx, ga3.Pump(ga3.Pump.HOT))
    mdot_c, dp_c = ga3.chicSolver(hx, ga3.Pump(ga3.Pump.COLD))
    Q = hx.thermalAnalysis(mdot_h, mdot_c)
    print(np.max(Q))
    plt.plot(Nt, Q)
    plt.show()


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
    MC, MH = [], []
    PC, PH = [], []
    hx = ga3.HX(COLDSTREAM, HOTSTREAM, KT, EPST, lt, DO, DI, Nt, Y, isSquare, Np, Nb[0], G, DS, DN)

    for nb in Nb:
        hx.Nb = nb
        mdot_h, dp_h = ga3.chicSolver(hx, ga3.Pump(ga3.Pump.HOT))
        mdot_c, dp_c = ga3.chicSolver(hx, ga3.Pump(ga3.Pump.COLD))
        print(mdot_c)
        print(mdot_c / (hx.coldStream["rho"] * hx.As))
        MC.append(mdot_c)
        MH.append(mdot_h)
        PC.append(dp_c)
        PH.append(dp_h)
        Q.append(hx.thermalAnalysis(mdot_h, mdot_c))

    Q = np.array(Q)
    MC = np.array(MC)
    MH = np.array(MH)
    PC = np.array(PC)
    PH = np.array(PH)

    plt.plot(Nb, PC)
    plt.ylim(0, None)

    plt.show()


if __name__ == "__main__":
<<<<<<< HEAD
    num_tubes()
    #print(hx.hydraulicAnalysisShe;l(mdot = 0.5, verbose=True))
=======
    num_passes()
>>>>>>> dbd3188fe97cbf808e6c9e137e26b1a7a7b5d1d5
