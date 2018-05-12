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
print('Stoffmenge Probe     ', n, ' mol')
# Loschmidtsche Zahl von CODATA
Nl = ufloat(2.6516467, 0.0000015)*1e25  # in 1/m^3
# longitudinale Phasengeschwindigkeit in Kupfer
vlong = 4.7*1e3                 # in m/s
# transversale Phasengeschwindigkeit in Kupfer
vtrans = 2.26*1e3               # in m/s
# Volumen der Probe
Vp = V0 * n               # in m^3
print('Volumen Probe        ', Vp, 'm^3')
# Avogadro-Konstante
Na = codata.value('Avogadro constant')
print('Avogadro-Konstante   ', Na, ' 1/mol')

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
print('a = ', params[0], ' +/- ', errors[0], ' 1/grd')
print('b = ', params[1], ' +/- ', errors[1], ' 1/grd^2')
print('c = ', params[2], ' +/- ', errors[2], ' 1/grd^3')
print('d = ', params[3], ' +/- ', errors[3], ' 1/grd^4')
print('e = ', params[4], ' +/- ', errors[4], ' 1/grd^5')
a = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])
c = ufloat(params[2], errors[2])
d = ufloat(params[3], errors[3])
e = ufloat(params[4], errors[4])
# par = np.array([a, b, c, d, e])

# Plotten von alpha
plt.figure()
plt.plot(Talpha, alpha*1e6, 'bx', label='Stützstellen')
Tplot = np.linspace(Talpha[0], Talpha[-1], 500)
plt.plot(Tplot, func_alpha(Tplot, *params)*1e6, 'k-', label='Regression')
plt.xlabel(r'$T\;\mathrm{in}\;\mathrm{°C}$')
plt.ylabel(r'$\alpha\;\mathrm{in}\;10^{-6}\mathrm{grd}^{-1}$')
# plt.xlim(60, 310)
# plt.ylim(6.70, 16.95)
plt.grid()
plt.tight_layout()
plt.savefig('build/alpha.pdf')
plt.clf()

# Berechne Cv mittels Korrekturformel
Tmittel = (iTp + fTp)/2             # in grd
Cv = Cp - 9 * func_alpha(Tmittel, a, b, c, d, e)**2 * kappa * V0 * Tmittel

# Plotten von Cv
plt.figure()
plt.errorbar(x=noms(Tmittel), xerr=stds(Tmittel), y=noms(Cv), yerr=stds(Cv), color='b', fmt='x', label='Stützstellen')
# Tmax = 170 - 273.15                 # in grd
# plt.axvline(x=Tmax, color='k')
plt.xlabel(r'$T\;\mathrm{in}\;\mathrm{°C}$')
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
print("Tabelle zum Bestimmen der Debye-Temperatur aus Tabelle in Anleitung")
table = pd.DataFrame({'Cv': noms(Cv),
                      'T': noms(Tmittel)})
print(table)

# theta-debye/T für die Messwerte bis 170 K
abgelesen = np.array([3.2, 0.8, 5.5, 8.3, 9.2, 8.8, 7.9, 6.8])
theta_exp = abgelesen * (Tmittel[0:len(abgelesen)]+273.15) - 273.15
print('Experimentelle Bestimmung')
print('    Theta-Debye  ', np.mean(theta_exp), '°C')


###
# Theoretische Bestimmung von omega-debye und theta-debye
###
print('Theoretische Berechnung')

omega_theo = 18*(np.pi**2)*Nl / (Vp * (1/(vlong**3) + 2/(vtrans**3)))
omega_theo = omega_theo**(1/3)
print('    Omega-Debye  ', omega_theo, ' 1/s')
theta_theo = omega_theo * codata.value('Planck constant over 2 pi') / codata.value('Boltzmann constant') # in K
theta_theo -= 273.15        # in grd
print('    Theta-Debye  ', theta_theo, ' °C')

print('Verwende Na*n anstelle von Nl')
omega_verzweiflung = 18*(np.pi**2)*Na*n / (Vp * (1/(vlong**3) + 2/(vtrans**3)))
omega_verzweiflung = omega_verzweiflung**(1/3)
theta_verzweiflung = omega_verzweiflung * codata.value('Planck constant over 2 pi') / codata.value('Boltzmann constant') # in K
theta_verzweiflung -= 273.15    # in °C
print('    Omega-Debye  ', omega_verzweiflung, ' 1/s')
print('    Theta-Debye  ', theta_verzweiflung, ' °C')


# Erstelle Tabellen
ascii.write(
            [np.round(noms(dt), 0),
             noms(U),
             np.round(noms(I)*1e3, 1), # in mA
             noms(iRp),
             noms(fRp),
             np.round(noms(iTp), 1),
             np.round(stds(iTp), 1),
             np.round(noms(fTp), 1),
             np.round(stds(fTp), 1)],
            'build/table_messwerte.tex',
            format='latex',
            overwrite=True)

ascii.write([
                np.round(noms(Cp), 2),
                np.round(stds(Cp), 2),
                np.round(noms(Tmittel), 1),
                np.round(stds(Tmittel), 1),
                np.round(noms(func_alpha(Tmittel, a, b, c, d, e)*1e6), 2),
                np.round(stds(func_alpha(Tmittel, a, b, c, d, e)*1e6), 2),
                np.round(noms(Cv), 2),
                np.round(stds(Cv), 2),
            ],
            'build/table_cv.tex',
            format='latex',
            overwrite=True)

ascii.write([
                np.round(noms(Tmittel[0:len(abgelesen)]), 1),
                np.round(stds(Tmittel[0:len(abgelesen)]), 1),
                np.round(noms(Cv[0:len(abgelesen)]), 2),
                np.round(stds(Cv[0:len(abgelesen)]), 2),
                abgelesen,
                np.round(noms(theta_exp), 1),
                np.round(stds(theta_exp), 1)
            ],
            'build/table_theta.tex',
            format='latex',
            overwrite=True)

# Nochmalige Berechnung von Theta-Debye durch Anpassung der C_V Werte.
# Für die Diskussion
shift = 16
Cv_neu = Cv[2:8] + shift
print('Wirds bei einem shift von Cv um ', shift, 'besser?')
print(Cv_neu)
abgelesen_neu = np.array([0.9, 2.4, 2.5, 2.5, 2.3, 1.9])
# shift 13 abgelesen_neu = np.array([1.9, 3.1, 3.3, 3.2, 3.0, 3.6])
theta_neu = abgelesen_neu * (Tmittel[2:8]+273.15)
theta_neu -= 273.15 # Umrechnen kelin in °C
print(theta_neu)
print('Es ergibt sich für Theta-Debye ', np.mean(theta_neu), ' °C')

ascii.write([
                np.round(noms(Tmittel[2:8]), 1),
                np.round(stds(Tmittel[2:8]), 1),
                np.round(noms(Cv_neu), 2),
                np.round(stds(Cv_neu), 2),
                abgelesen_neu,
                np.round(noms(theta_neu), 1),
                np.round(stds(theta_neu), 1)
            ],
            'build/table_neu.tex',
            format='latex',
            overwrite=True)
