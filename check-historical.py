import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import main as ga3
from GA3_CONSTS import *


def main():
    historical_data = pd.read_csv("data/tests.csv")
    out = np.zeros(len(historical_data.index))
    out_preds = np.zeros_like(out)
    var = np.zeros_like(out)
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

        Thi = 60
        Tci = 20

        Tho = Thi - row.Q / mdot_h / HOTSTREAM["cp"]
        Tco = Tci + row.Q / mdot_c / COLDSTREAM["cp"]

        T1 = Thi - Tco
        T2 = Tho - Tci
        LMTD = T1 if np.isclose(T1, T2) else (T1 - T2) / np.log(T1 / T2)

        Q = hx.thermalAnalysis_NTU(mdot_h, mdot_c)

        var[i] = i
        out[i] = row.Q if row.year != 2023 else Q
        out_preds[i] = Q
        # plt.annotate(f'{row.Nps}{row.Npt}', (i, max(out[i], out_preds[i]) * 1.05))

        if row.year != year_last:
            plt.vlines(float(i) - 0.5, 0, 0.8, 'r')
            plt.annotate(row.year, (i, 0))
            year_last = row.year

        plt.plot([i, i], [out[i], out_preds[i]], 'gray', zorder=-1)

    print(np.sum((out - out_preds) ** 2))

    plt.scatter(var, out, label='data')
    plt.scatter(var, out_preds, label='predicted')
    plt.legend()
    # plt.ylim(-0.5, 0.5)
    # for i in historical_data.index:
    # plt.annotate(i,(var[i],dQs[i]))
    plt.title('Historical Data')
    plt.xlabel('id')
    plt.ylabel('$Q$')
    plt.show()


if __name__ == "__main__":
    main()
