import main as ga3

# INLET FLUID PROPERTIES
coldStream = {"Tci": 20,  # C
              "cp": 4179,  # J/kgK
              "rho": 990.1,  # kg/m^3
              "k": 0.632,  # W/mK
              "mu": 6.51E-4}  # kg/ms

hotStream = {"Tci": 60,  # C
             "cp": 4179,  # J/kgK
             "rho": 990.1,  # kg/m^3
             "k": 0.632,  # W/mK
             "mu": 6.51E-4}  # kg/ms

# TUBE PROPERTIES
epst = 0.0000015  # Effective roughness height for new drawn copper, mm
kt = 386 # Thermal conductivity of copper, W/mK
do, di = 0.008, 0.006 # Tube outer and inner diameter, m

# SHELL PROPERTIES
ds = 0.064 # Shell inner diameter, m

# NOZZLE PROPERTIES
dn = 0.02 # Nozzle diameter, m

# FREE VARIABLES
lt = 0.35 # Tube length, m
Nt = 13 # Tube count
Y = 0.014 # Tube pitch, m NOT YET USED
isSquare = True # Square or triangular tube pattern NOT YET USED

Np = 1 # Number of passes
Nb = 9 # Number of shell baffles NOT YET USED
B = lt/(Nb + 1) # Baffle pitch NOT YET USED

G = 0.2*ds # Baffle cut, m NOT YET USED

# Test heat exchanger using worked example values
HXtest = ga3.HX(coldStream, hotStream, kt, epst, lt, do, di, Nt, Y, isSquare, Np, Nb, B, G, ds, dn)
HXtest.hydraulicAnalysisTube(mdot=0.45, verbose=True) # Worked example tube pressure drop 4360 Pa