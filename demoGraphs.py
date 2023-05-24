import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import numpy as np
import pandas as pd

import main as ga3
from GA3_CONSTS import *

coldPump = ga3.Pump(ga3.Pump.COLD)
hotPump = ga3.Pump(ga3.Pump.HOT)

Nt = 13 # number of tubes
lt = 3.5  # tube lengths
Y = 0.014 # tube pitch

Nb = 9  # Number of shell baffle
B = lt / (Nb + 1)  # Baffle pitch

isSquare = False
Nps = 1
Npt = 1

G = 0.2*DS

tube_arr = pd.read_csv("data/tube-arrangements.csv")
tube_interp = interp1d(tube_arr["Nt"].to_numpy(), tube_arr["Y"].to_numpy())

# Vary tube count
NtRange = np.arange(6, 20, 1)
Nhx = len(NtRange)

YRange = (DS / 0.064) * tube_interp(NtRange * Npt) / 1000 # Corresponding maximum pitch
ltRange = np.clip(3.5 / (NtRange * Npt), 0, LT_MAX_1P) # Limited by total copper length and HX overall length

hxList = []
QList = np.zeros(Nhx)
mdothList = np.zeros_like(QList)
mdotcList = np.zeros_like(QList)
#dPList = np.zeros_like(hxList)

for i in range(Nhx):
    # Initialise HX object
    hxList.append(ga3.HX(COLDSTREAM, HOTSTREAM, KT, EPST, ltRange[i], DO, DI, NtRange[i], YRange[i], isSquare, Nps, Npt, Nb, G, DS, DN))

    # Solve for mass flows with given pumps
    mdotcList[i] = ga3.chicSolver(hxList[i], coldPump)[0] # Tube
    mdothList[i] = ga3.chicSolver(hxList[i], hotPump)[0] # Shell

    QList[i]= hxList[i].thermalAnalysis_LMTD(mdot_t = mdotcList[i], mdot_s = mdothList[i], verbose = False)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.plot(NtRange, QList/1000, color="orange", label = "$Q$")
ax2.plot(NtRange, mdotcList, color="blue", label="$\dot{m}_{cold}$")
ax2.plot(NtRange, mdothList, color="red", label="$\dot{m}_{hot}$")

ax1.set_xlabel("Tube count $Nt$")
ax1.set_ylabel("Heat transfer, kW")
ax2.set_ylabel("Mass flow rate, kg/s")
ax1.set_title("Reference design modified to fit maximum copper length constraint, with varying tube count.\n Tube pitch maximised in all cases.")

ax1.set_xlim(0)
ax1.set_xticks(np.arange(0, np.max(NtRange)+1, 1))
ax1.set_ylim(0)
ax2.set_ylim(0)
ax1.grid()
ax1.legend(loc=0)
ax2.legend(loc=2)
plt.show()