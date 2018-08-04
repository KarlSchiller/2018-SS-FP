import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.optimize import curve_fit
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.constants import codata


def g(x, m, b):
    return m*x + b

def g_1g_2(L,r_1,r_2):
    return (1-L/r_1)*(1-L/r_2)

def quadrat(x, a, b, c):
    return a*x**2+b*x+c

r_e = 1e15
r_k1 = 100
r_k2 = 140

#-------------------konkav-konkav
def evak_Tkk():
    L, I = np.genfromtxt("rohdaten/spiegel_rund.txt", unpack=True)

    c = g_1g_2(47.3, r_k2, r_k2)
    test = (I*c/max(I))

    params_kk, covariance_kk = curve_fit(quadrat, L, test)
    errors_kk = np.sqrt(np.diag(covariance_kk))

    print('Skalierungsfaktor: ', c)
    print('Quadrat: ', ufloat(params_kk[0], errors_kk[0]))
    print('Linear: ', ufloat(params_kk[1], errors_kk[1]))
    print('Knostante: ', ufloat(params_kk[2], errors_kk[2]))

    # Plot
    x = np.linspace(L[0]-1, L[-1]+1)
    plt.plot(L, test, 'kx', label='Messwerte')
    plt.plot(x, quadrat(x, *params_kk), label=r'Fit')
    plt.plot(x, g_1g_2(x, r_k2, r_k2), label=r'$\mathrm{Theorie:} r_{k2}, r_{k2}$')
    plt.xlabel(r'$d /\ mm$')
    plt.ylabel(r'$g_1g_2$')
    # plt.xlim(60, 310)
    # plt.ylim(6.70, 16.95)
    plt.grid()
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('build/stab-rund.pdf')
    plt.clf()

#-------------------konkav-flach
def evak_Tkf():
    L2, I2 = np.genfromtxt("rohdaten/spiegel_flach.txt", unpack=True)

    c_kf = g_1g_2(92.3, r_e, r_k2)
    test = (I2*c_kf)/max(I2)

    params_kf, covariance_kf = curve_fit(g, L2, test)
    errors_kf = np.sqrt(np.diag(covariance_kf))

    print('Skalierungsfaktor:', c_kf)
    print('linear: ', ufloat(params_kf[0], errors_kf[0]))
    print('konstante: ', ufloat(params_kf[1], errors_kf[1]))

    # Plot
    x = np.linspace(L2[0]-1, L2[-1]+1)
    plt.plot(L2, test, 'kx', label='Messwerte')
    plt.plot(x, g(x, *params_kf), '.', label=r'Fit')
    plt.plot(x, g_1g_2(x, r_e, r_k2), label=r'$\mathrm{Theorie:} \, flach, r_{k1}$')
    plt.xlabel(r'$d /\ mm$')
    plt.ylabel(r'$g_1g_2$')
    # plt.xlim(60, 310)
    # plt.ylim(6.70, 16.95)
    plt.grid()
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('build/stab-flach.pdf')
    plt.clf()

if __name__ == '__main__':

    evak_Tkk()
    evak_Tkf()
