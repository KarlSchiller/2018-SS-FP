import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.optimize import curve_fit


# Umrechnung R/Ohm zu T/K
def temperature(R):
    return 0.00134*R**2 + 2.296*R - 243.02 + 273.15


# T_zylinder1, T_zylinder2, T_puppe1, T_puppe2 = np.genfromtxt('rohdaten/data3.txt', unpack=True)

# ascii.write(
#             [T_zylinder1, noms(T_z1), stds(T_z1), T_zylinder2, noms(T_z2), stds(T_z2)],
#             'build/table_zylinder.tex', format='latex')
