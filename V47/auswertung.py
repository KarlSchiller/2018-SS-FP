import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.optimize import curve_fit
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import uncertainties.unumpy as unp
from uncertainties import ufloat
import pandas as pd
from scipy.constants import codata

# Umrechnung R/Ohm zu T/grd
def temperature(R):
    return 0.00134*R**2 + 2.296*R - 243.02


# Fitfunktion für den linearen Ausdehnungskoeffizienten
def func_alpha(T, a, b, c, d, e):
    return a + b*T + c*T**2 + d*T**3 + e*T**4


# Masse der Probe
m = 0.342                      # in kg
# Kompressionsmodul Kupfer
kappa = 139e9                   # in N/m^2
# Molvolumen Kupfer
V0 = 7.11e-6                    # in m^3/mol
# Molare Masse Kupfer
M = 63.55*1e-3                 # in kg/mol
# Stoffmenge der Probe
n = m / M                    # in mol
# Loschmidtsche Zahl von CODATA
Nl = ufloat(2.6516467, 0.0000015)*1e25
# longitudinale Phasengeschwindigkeit in Kupfer
vlong = 4.7*1e3                 # in m/s
# transversale Phasengeschwindigkeit in Kupfer
vtrans = 2.26*1e3               # in m/s


# Einlesen der gemessenen Werte
# Notation: r steht für read, also eingelesen
# Endung p steht für probe, z für Zylinder
rdt, rU, rI, riRp, riRz, rfRp, rfRz = np.genfromtxt('rohdaten/messung.txt', unpack=True)
dt = unp.uarray(rdt, 3)         # in s
U = unp.uarray(rU, 0.01)      # in V
I = unp.uarray(rI, 0.3)*1e-3  # in A
iRp = unp.uarray(riRp, 0.1)     # in ohm
iRz = unp.uarray(riRz, 0.1)     # in ohm
fRp = unp.uarray(rfRp, 0.1)     # in ohm
fRz = unp.uarray(rfRz, 0.1)     # in ohm

# Berechne die Temperaturen aus den Widerständen
iTp = temperature(iRp)          # in grd
iTz = temperature(iRz)          # in grd
fTp = temperature(fRp)          # in grd
fTp = temperature(fRz)          # in grd


###
# 5.a) Berechne Cp(T)
###


# Berechne elektrische Arbeit pro Zeitintervall
Wel = U * I * dt                # in J

# Berechne Cp aus den Messwerten
Cp = Wel / (np.absolute(fTp - iTp) * n)    # in J/(grd*mol)


###
# 5.b) Berechne Cv(T)
###


# Einlesen der Werte von alpha, aus der Tabelle der Anleitung
Talpha, alpha = np.genfromtxt('rohdaten/alpha.txt', unpack=True)
Talpha -= 273.15                # Umrechnen K in grd
alpha *= 1e-6                    # in 1/grd

# Bestimmung einer allgemeinen Funktion von alpha
print('Regression für alpha')
params, covariance = curve_fit(func_alpha, Talpha, alpha)
# covariance is the covaniance matrix
errors = np.sqrt(np.diag(covariance))
# print('a = ', params[0], ' +/- ', errors[0], ' 1/grd')
# print('b = ', params[1], ' +/- ', errors[1])
print('a = ', params[0], ' +/- ', errors[0], ' 1/grd')
print('b = ', params[1], ' +/- ', errors[1], ' 1/grd^2')
print('c = ', params[2], ' +/- ', errors[2], ' 1/grd^3')
print('d = ', params[3], ' +/- ', errors[3], ' 1/grd^4')
print('e = ', params[4], ' +/- ', errors[4], ' 1/grd^5')

# Plotten von alpha
plt.figure()
plt.plot(Talpha, alpha*1e6, 'bx', label='Stützstellen')
Tplot = np.linspace(Talpha[0], Talpha[-1], 500)
plt.plot(Tplot, func_alpha(Tplot, *params)*1e6, 'k-', label='Regression')
plt.xlabel(r'$T\;\mathrm{in}\;\mathrm{K}$')
plt.ylabel(r'$\alpha\;\mathrm{in}\;10^{-6}\mathrm{grd}^{-1}$')
# plt.xlim(60, 310)
# plt.ylim(6.70, 16.95)
plt.grid()
plt.tight_layout()
plt.savefig('build/alpha.pdf')
plt.clf()

# Berechne Cv mittels Korrekturformel
Tmittel = (iTp + fTp)/2             # in grd
Cv = Cp - 9 * func_alpha(Tmittel, *params)**2 * kappa * V0 * Tmittel

# Plotten von Cv
Tmax = 170 - 273.15                 # in grd
plt.figure()
plt.errorbar(x=noms(Tmittel), xerr=stds(Tmittel), y=noms(Cv), yerr=stds(Cv), color='b', fmt='x', label='Stützstellen')
plt.axvline(x=Tmax, color='k')
plt.xlabel(r'$T\;\mathrm{in}\;\mathrm{K}$')
plt.ylabel(r'$C_{\mathrm{V}}\;\mathrm{in}\;\mathrm{J/(mol\;grd)}$')
# plt.xlim(60, 310)
# plt.ylim(6.70, 16.95)
plt.grid()
plt.tight_layout()
plt.savefig('build/cv.pdf')
plt.clf()


###
# 5.c) Experimentelle Bestimmung von Cv
###


# Print (Cv, T)
table = pd.DataFrame({'Cv': noms(Cv),
                      'T': noms(Tmittel)})
print(table)

# TODO


###
# Theoretische Bestimmung von omega-debye und theta-debye
###
print('Theoretische Berechnung')

# Volumen der Probe
Vp = V0 * m / M          # in m^3

omega_theo = 18*np.pi**2*Nl / (Vp * (1/vlong**3 + 2/vtrans**3))
omega_theo = omega_theo**(1/3)
print('Omega-Debye  ', omega_theo)
theta_theo = omega_theo * codata.value('Planck constant over 2 pi') / codata.value('Boltzmann constant')
print('Theta-Debye  ', theta_theo)

# ascii.write(
#             [T_zylinder1, noms(T_z1), stds(T_z1), T_zylinder2, noms(T_z2), stds(T_z2)],
#             'build/table_zylinder.tex', format='latex')
