import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import main as ga3
from GA3_CONSTS import *


def main():
    historical_data = pd.read_csv("data/tests.csv")
    Qs = np.zeros(len(historical_data.index))
    Q_preds = np.zeros_like(Qs)
    var = np.zeros_like(Qs)
    G = 0.2 * DS

    hx = ga3.HX(COLDSTREAM, HOTSTREAM, KT, EPST, 0.25, DO, DI, 12, 0.01, False, 1, 1, 12, 0.015, G, DS, DN)

    for i, row in historical_data.iterrows():
        group, year, N, Y, Nb, B, L, Npt, Nps, flow_h, flow_c, Q, eff = row
        hx.Nt = N
        hx.Y = Y
        hx.Nb = Nb
        hx.B = B
        hx.lt = L
        hx.Npt = Npt
        hx.Nps = Nps

        hpump = ga3.Pump(ga3.Pump.HOT, year)
        cpump = ga3.Pump(ga3.Pump.COLD, year)

        mdot_h, dp_h = hx.chicSolver(hpump)
        mdot_c, dp_c = hx.chicSolver(cpump)

        Q_pred = hx.thermalAnalysis(mdot_h, mdot_c, False)

        Qs[i] = Q
        var[i] = Nps
        Q_preds[i] = Q_pred

    dQs = (Q_preds - Qs) / Qs

    plt.scatter(var, dQs)

    for i in historical_data.index:
        plt.annotate(i,(var[i],dQs[i]))

    plt.show()


if __name__ == "__main__":
    main()
