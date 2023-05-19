import matplotlib.pyplot as plt
import numpy as np

from main import Pump


def plot_pump_data(year):
    cold = Pump(Pump.COLD,year)
    hot = Pump(Pump.HOT,year)

    n = 100
    rhoHot = rhoCold = 940

    coldmdots = np.linspace(cold.flowMin * rhoCold, cold.flowMax * rhoCold, n)
    hotmdots = np.linspace(hot.flowMin * rhoHot, hot.flowMax * rhoHot, n)

    colddps = np.zeros_like(coldmdots)
    hotdps = np.zeros_like(hotmdots)

    for i in range(n):
        colddps[i] = cold.dp(coldmdots[i] / rhoCold)
        hotdps[i] = hot.dp(hotmdots[i] / rhoHot)

    plt.plot(coldmdots, colddps / 1000, label="Cold fit", color="blue")
    plt.scatter(cold.flowrate_data * rhoCold, cold.dp_data / 1000, label="Cold data", marker="x", color="blue")
    plt.plot(hotmdots, hotdps / 1000, label='Hot fit', color="red")
    plt.scatter(hot.flowrate_data * rhoHot, hot.dp_data / 1000, label="Hot data", marker="x", color="red")
    plt.xlabel("Pump mass flow rate, $kg/s$")
    plt.ylabel("Pump pressure rise, $kPa$")
    plt.legend()
    plt.grid()
    plt.xlim([np.amin((coldmdots, hotmdots)), np.amax((coldmdots, hotmdots))])
    plt.title("Pump characteristics")


if __name__ == "__main__":
    [plot_pump_data(y) for y in (2017,2018,2020,2022)]
    plt.show()
