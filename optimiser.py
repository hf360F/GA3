import pandas as pd
from scipy.interpolate import interp1d

import main as ga3
from GA3_CONSTS import *


def num_passes():
    Npt = 2
    Nps = 1
    Nt = 12
    lt = 0.25 - 2 * END_WASTAGE  # tube lengths

    tube_arr = pd.read_csv("data/tube-arrangements.csv")
    tube_interp = interp1d(tube_arr["Nt"].to_numpy(), tube_arr["Y"].to_numpy())
    Y = (DS / 0.064) * tube_interp(Nt) / 1000  # tube pitch (reference diameter for data weas 0.064)

    Nb = 12
    B = lt / (Nb + 1)

    isSquare = False
    G = 0.2 * DS

    hx = ga3.HX(COLDSTREAM, HOTSTREAM, KT, EPST, lt, DO, DI, Nt, Y, isSquare, Nps, Npt, Nb, B, G, DS, DN)

    hpump_22 = ga3.Pump(ga3.Pump.HOT, 2022)
    cpump_22 = ga3.Pump(ga3.Pump.COLD, 2022)

    mdot_h, dp_h = hx.chicSolver(hpump_22)
    mdot_c, dp_c = hx.chicSolver(cpump_22)

    hx.thermalAnalysis_NTU(mdot_h, mdot_c, verbose=True)
    hx.thermalAnalysis_LMTD(mdot_h, mdot_c, verbose=True)


if __name__ == "__main__":
    num_passes()
