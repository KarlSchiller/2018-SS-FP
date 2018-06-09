import numpy as np
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.constants import codata


# Volumen eines Zylinders
# INPUT Durchmesser d in m
#       Länge l in m
# OUTPUT Volumen in m^3
def VolZyl(d, l):
    return np.pi*(d/2)**2*l


# Umrechnung Liter in Kubikmeter
def liter2qmeter(V_liter):
    return V_liter*1e-3


# selbst vermessen:
print('Selbst vermessene Bauteile')
# Schlauch S2
d_S2 = ufloat(24.75, 0.1)*1e-3  # m
l_S2 = ufloat(1210, 20)*1e-3  # m
V_S2 = VolZyl(d_S2, l_S2)  # m^3
print('V_S2     {:} m^3'.format(V_S2))
# Verbindungsstück zwischen S2 und B5
d_E1 = ufloat(11.95, 0.1)*1e-3  # m
l_E1 = ufloat(59.90, 0.1)*1e-3  # m
V_E1 = VolZyl(d_E1, l_E1)  # m^3
print('V_E1     {:} m^3'.format(V_E1))
# Kreuzstück zwischen B2 und B3
d_B6 = ufloat(40.5, 0.1)*1e-3  # m
l1_B6 = ufloat(130.0, 0.1)*1e-3  # m
l2_B6 = ufloat(45.0, 5.0)*1e-3  # m
V_B6 = VolZyl(d_B6, l1_B6)+VolZyl(d_B6, l2_B6)  # m^3
print('V_B6     {:} m^3'.format(V_B6))
# Schlauch S1
d_S1 = ufloat(16, 1)*1e-3  # m
l_S1 = ufloat(430, 10)*1e-3  # m
V_S1 = VolZyl(d_S1, l_S1)  # m^3
print('V_S1     {:} m^3'.format(V_S1))

# Aus Anleitung:
V_B1 = liter2qmeter(ufloat(9.5, 0.8))  # m^3
V_B3 = liter2qmeter(ufloat(0.25, 0.01))  # m^3
V_B5 = liter2qmeter(ufloat(0.016, 0.002))  # m^3
V_B2 = liter2qmeter(ufloat(0.177, 0.09))  # m^3
V_V4 = liter2qmeter(ufloat(0.025, 0.005))  # m^3
V_V1auf = liter2qmeter(ufloat(0.044, 0.004))  # m^3
V_V1zu = liter2qmeter(ufloat(0.022, 0.002))  # m^3
V_B4 = liter2qmeter(ufloat(0.067, 0.004))  # m^3
V_V5zu = liter2qmeter(ufloat(0.005, 0.001))  # m^3
V_V5auf = liter2qmeter(ufloat(0.015, 0.002))  # m^3 nicht verwendet
V_V2auf = V_V5auf
V_V2zu = V_V5zu
V_V3auf = V_V5auf
V_V3zu = V_V5zu

# Berechnung der Volumina
V_D_leck = V_B1 + V_S2 + V_S1 + V_B3 + V_B5 + V_B2 + V_V5zu + V_V4 + \
            V_V1zu + V_E1 + V_V2auf + V_B6 + V_V3auf
V_D_evak = V_B1 + V_S2 + V_S1 + V_B3 + V_B5 + V_B2 + V_V5zu + V_V4 + \
            V_V1zu + V_E1 + V_V2auf + V_B6 + V_V3zu
V_T_leck = V_B1 + V_B3 + V_B2 + V_V1auf + V_B4 + V_V2zu + V_B6 + V_V2auf
V_T_evak = V_B1 + V_B3 + V_B2 + V_V1auf + V_B4 + V_V2zu + V_B6 + V_V2zu
print('Berechnete Volumina')
print('Vol D Leck     {:} m^3'.format(V_D_leck))
print('Vol D Evak     {:} m^3'.format(V_D_evak))
print('Vol T Leck     {:} m^3'.format(V_T_leck))
print('Vol T Evak     {:} m^3'.format(V_T_evak))

# Speichere für weitere Auswertung
with open("build/volumina.txt", 'w+') as file:
    file.write("{:} {:} {:} {:}".format(noms(V_D_leck),
                                        stds(V_D_leck),
                                        noms(V_T_leck),
                                        stds(V_T_leck)))
