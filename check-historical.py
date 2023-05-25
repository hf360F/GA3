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
    year_last = 0

    for i, row in historical_data.iterrows():
        hx.Nt = row.N
        hx.Y = row.Y
        hx.Nb = row.Nb
        hx.B = row.B
        hx.lt = row.L
        hx.Npt = row.Npt
        hx.Nps = row.Nps

        hpump = ga3.Pump(ga3.Pump.HOT, row.year)
        cpump = ga3.Pump(ga3.Pump.COLD, row.year)

        mdot_h, dp_h = hx.chicSolver(hpump)
        mdot_c, dp_c = hx.chicSolver(cpump)

        Q_pred = hx.thermalAnalysis_LMTD(mdot_h, mdot_c, False)
        Qs[i] = row.Q
        var[i] = i
        Q_preds[i] = Q_pred
        if row.year != year_last:
            plt.vlines(float(i) - 0.5, 0.500, 16.000, 'r')
            plt.annotate(row.year, (i, 1.000))
            year_last = row.year
        plt.plot([i, i], [row.Q/1000, Q_pred/1000], 'gray', zorder=-1)

    dQs = (Q_preds - Qs) / Qs

    plt.scatter(var, Qs/1000, label='data')
    plt.scatter(var, Q_preds/1000, label='predicted')
    plt.legend()
    plt.ylim(0,None)
    # for i in historical_data.index:
    # plt.annotate(i,(var[i],dQs[i]))
    plt.title('Comparison Against Historical Data')
    plt.xlabel('test id')
    plt.ylabel('$Q$ (kW)')
    plt.show()


if __name__ == "__main__":
    main()
