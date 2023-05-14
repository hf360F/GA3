import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


def main():
    # define solution space
    Nt = np.arange(13, 48)
    lt = 0.27 / Nt
    tube_arr = pd.read_csv("data/tube-arrangements.csv")
    Y = np.interp(Nt, tube_arr["Nt"].to_numpy(), tube_arr["Y"].to_numpy())




if __name__ == "__main__":
    main()
