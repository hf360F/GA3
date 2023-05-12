from main import Pump
import matplotlib.pyplot as plt
import numpy as np

cold = Pump(Pump.COLD)
hot = Pump(Pump.HOT)

Qs = np.linspace(0, 0.0009, 1000)
plt.plot(Qs, cold.dp(Qs), label='cold')
plt.plot(Qs, hot.dp(Qs), label='hot')
plt.legend()
plt.xlim([0, None])
plt.ylim([0, None])
plt.show()
