# INLET FLUID PROPERTIES
COLDSTREAM = {"Ti": 20,  # C
              "cp": 4179,  # J/kgK
              "rho": 990.1,  # kg/m^3
              "k": 0.632,  # W/mK
              "mu": 6.51E-4}  # kg/ms

HOTSTREAM = {"Ti": 60,  # C
             "cp": 4179,  # J/kgK
             "rho": 990.1,  # kg/m^3
             "k": 0.632,  # W/mK
             "mu": 6.51E-4}  # kg/ms

# TUBE PROPERTIES
EPST = 300e-6  # Effective roughness height for new drawn copper, mm
KT = 386  # Thermal conductivity of copper, W/mK
DO, DI = 0.008, 0.006  # Tube outer and inner diameter, m
LT_MAX_1P = 0.27  # max length of individual cut tubes, m
LT_MAX_2P = 0.34  # max length of individual cut tubes, m
LT_TOTAL = 3.5  # max aggregate length of all copper tubing, m
END_WASTAGE = 0.005  # wasted tube sitting in seals on endcap, 5mm per endcap

# SHELL PROPERTIES
DS = 0.064  # Shell inner diameter, m

# NOZZLE PROPERTIES
DN = 0.02  # Nozzle diameter, m

M_TUBE = 0.2  # kg/m
M_SHELL = 0.65  # kg/m
M_NOZZLE = 0.025  # kg
M_BAFFLE = 2.39  # kg/m2
M_ENDS = 1150  # kg/m3
M_ORING = 0.0053  # kg
