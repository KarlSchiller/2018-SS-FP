import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.optimize import curve_fit, minimize_scalar
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.constants import codata


def feldfit(x, a, b, c, d, e):
    return a*x**4 + b*x**3 + c*x**2 + d*x + e


def toRad(phi):
    return 2*np.pi*phi/(360)


def linear(x, a, b):
    return a*x + b


# INPUT     N Elektronendichte in (cm)^(-3)
def theoriemasse(N):
    return 0.0635 + 2.06e-22*N + 1.16e-40*N**2


def feldmessung():
    print('Auswertung Feldmessung')
    # Volumina der Rezipienten
    B, z = np.genfromtxt("rohdaten/feldmessung2.txt", unpack=True)
    z0 = 100  # Überlapp in mm
    z -= z0  # Nullposition nun bei z=0

    print(' quadratische Regression an Messwerte')
    params, covariance = curve_fit(f=feldfit, xdata=z, ydata=B)
    # covariance is the covaniance matrix
    errors = np.sqrt(np.diag(covariance))
    print(' a = ', params[0], ' +/- ', errors[0], ' mT/(mm)^4')
    print(' b = ', params[1], ' +/- ', errors[1], ' mT/(mm)^3')
    print(' c = ', params[2], ' +/- ', errors[2], ' mT/(mm)^2')
    print(' d = ', params[3], ' +/- ', errors[3], ' mT/mm')
    print(' e = ', params[4], ' +/- ', errors[4], ' mT')

    # Finde Maximum des Feldes
    fm = lambda x: -feldfit(x, *params)
    r = minimize_scalar(fm, bounds=(-10, 10))
    print(" Erfolg der numerischen Maximumbestimmung:", r["success"])
    print(" Maximum bei", r["x"], "mm und", -r["fun"], "mT")

    # Plotte Messwerte und Regression
    zplot = np.linspace(-10, 10, 1000)
    plt.plot(zplot, feldfit(zplot, *params), 'b-', label='Lineare Regression')
    plt.plot(z, B, 'kx', label='Messpunkte')
    plt.xlabel(r'$z\;/\;\mathrm{mm}$')
    plt.ylabel(r'$B\;/\;\mathrm{mT}$')
    # plt.xlim(60, 310)
    # plt.ylim(6.70, 16.95)
    plt.grid()
    plt.tight_layout()
    # No whitespace around plots
    plt.savefig('build/feldmessung2.pdf', bbox_inches='tight')
    plt.clf()

    # Erstelle Tabelle
    ascii.write(
                [z,
                B],
                output='build/feldmessung2.tex',
                format='latex',
                overwrite=True)

    # Gebe maximale Flussdichte zurück
    return -r["fun"]


# INPUT:    B anliegende magnetische Flussdichte in mT
def effMasse(B):
    print("Bestimmung der effektiven Massen")
    # Probe 1
    L1 = 1.296e-3  # Dicke in m
    N1 = 2.8e18  # in (cm)^(-3)
    wavel, t11, t11m, t12, t12m = np.genfromtxt("rohdaten/probe1.txt", unpack=True)
    wavel *= 1e-3  # Umrechnung in µm
    theta1 = toRad((t11 + t11m / 60) - (t12 + t12m / 60))  # Winkel in rad
    theta1 = theta1 / L1  # auf Probendicke normierte Winkel

    # Probe 2
    L2 = 1.36e-3  # Dicke in m
    N2 = 1.2e18  # in (cm)^(-3)
    wavel, t21, t21m, t22, t22m = np.genfromtxt("rohdaten/probe1.txt", unpack=True)
    wavel *= 1e-3  # Umrechnung in µm
    theta2 = toRad((t21 + t21m / 60) - (t22 + t22m / 60))  # Winkel in rad
    theta2 = theta2 / L2  # auf Probendicke normierte Winkel

    # Probe 3
    L3 = 5.11e-3  # Dicke in m
    wavel, t31, t31m, t32, t32m = np.genfromtxt("rohdaten/probe1.txt", unpack=True)
    wavel *= 1e-3  # Umrechnung in µm
    theta3 = toRad((t31 + t31m / 60) - (t32 + t32m / 60))  # Winkel in rad
    theta3 = theta3 / L3  # auf Probendicke normierte Winkel

    # Alles gegen das Quadrat der Wellenlänge auftragen
    wavesquare = wavel**2  # in µm²
    # Differenzen der gewichteten Winkel
    diff1 = theta1 - theta3  # Für Probe 1
    diff2 = theta2 - theta3  # Für Probe 2

    print(' Regression Probe1')
    params1, covariance1 = curve_fit(f=linear, xdata=wavesquare[:-2], ydata=diff1[:-2])
    errors1 = np.sqrt(np.diag(covariance1))
    print(' m = ', params1[0], ' +/- ', errors1[0], ' *1e12 °/m^3')
    print(' b = ', params1[1], ' +/- ', errors1[1], ' *1e12 °/m')
    a1 = ufloat(params1[0], errors1[0])*1e12 # 1e12, da wavesquare in µm²
    m1_eff = unp.sqrt(e**3*N1*B/(8*np.pi**2*epsilon0*c**3*a1*m0**2))
    m1_theo = theoriemasse(N1)
    print(' m1_eff/m0 = ', m1_eff)
    print(' Theoriewert = ', m1_theo)
    print(' Abweichung ', (m1_eff-m1_theo)/m1_eff)

    print(' Regression Probe2')
    params2, covariance2 = curve_fit(f=linear, xdata=wavesquare[:-2], ydata=diff2[:-2])
    errors2 = np.sqrt(np.diag(covariance2))
    print(' m = ', params2[0], ' +/- ', errors2[0], ' ')
    print(' b = ', params2[1], ' +/- ', errors2[1], ' ')
    a2 = ufloat(params2[0], errors2[0])*1e12  # 1e12, da wavesquare in µm²
    m2_eff = unp.sqrt(e**3*N2*B/(8*np.pi**2*epsilon0*c**3*a2*m0**2))
    m2_theo = theoriemasse(N2)
    print(' m2_eff/m0 = ', m2_eff)
    print(' Theoriewert = ', m2_theo)
    print(' Abweichung ', (m2_eff-m2_theo)/m2_eff)

    # Plotten der Faraday-Winkel einzelnd
    plt.plot(wavesquare, theta1, 'rx', label='Probe 1')
    plt.plot(wavesquare, theta2, 'bx', label='Probe 2')
    plt.plot(wavesquare, theta3, 'kx', label='Probe 3')
    plt.xlabel(r'$\lambda^2\;/\;(\mathrm{µm})^2$')
    plt.ylabel(r'$\frac{\theta}{L}\;/\;\mathrm{rad\:m}^{-1}$')
    # plt.xlim(60, 310)
    # plt.ylim(6.70, 16.95)
    plt.grid()
    plt.tight_layout()
    plt.legend(loc='best')
    # No whitespace around plots
    plt.savefig('build/winkel.pdf', bbox_inches='tight')
    plt.clf()

    # Plotte Messwerte und Regression
    waveplot = np.linspace(wavesquare[0], wavesquare[-2], 10)
    plt.plot(waveplot, linear(waveplot, *params1), 'r-', label='Regression 1')
    plt.plot(waveplot, linear(waveplot, *params2), 'b-', label='Regression 2')
    plt.plot(wavesquare, diff1, 'rx', label='Probe 1')
    plt.plot(wavesquare, diff2, 'bx', label='Probe 2')
    plt.xlabel(r'$\lambda^2\;/\;(\mathrm{µm})^2$')
    plt.ylabel(r'$\frac{\theta}{L}\;/\;\mathrm{rad\:m}^{-1}$')
    # plt.xlim(60, 310)
    # plt.ylim(6.70, 16.95)
    plt.grid()
    plt.tight_layout()
    plt.legend(loc='best')
    # No whitespace around plots
    plt.savefig('build/differenzen.pdf', bbox_inches='tight')
    plt.clf()

    # Speichere Messwerte
    ascii.write(
                [np.round(wavel, 2),
                np.round(wavesquare, 2),
                np.round((theta1*L1)*360/np.pi, 2),
                np.round(theta1, 2),
                np.round(diff1, 2)],
                output='build/probe1.tex',
                format='latex',
                overwrite=True)
    ascii.write(
                [np.round(wavel, 2),
                np.round(wavesquare, 2),
                np.round((theta2*L2)*360/np.pi, 2),
                np.round(theta2, 2),
                np.round(diff2, 2)],
                output='build/probe2.tex',
                format='latex',
                overwrite=True)
    ascii.write(
                [np.round(wavel, 2),
                np.round(wavesquare, 2),
                np.round((theta3*L3)*360/np.pi, 2),
                np.round(theta3, 2)],
                output='build/probe3.tex',
                format='latex',
                overwrite=True)


if __name__ == '__main__':
    # Physikalischer Konstante
    e = codata.value('elementary charge')
    print('e', e)
    epsilon0 = codata.value('electric constant')
    print('eps0', epsilon0)
    c = codata.value('speed of light in vacuum')
    print('c', c)
    m0 = codata.value('electron mass')
    print('m0', m0)
    # Brechungsindex GaAs
    n = 3.397

    if not os.path.isdir('build'):
        os.mkdir('build')

    # Anliegende magnetische Flussdichte
    B = feldmessung()  # in mT
    # B = 376.17136964908855

    # Bestimmung der effektiven Masse
    effMasse(B)
