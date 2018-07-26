import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.optimize import curve_fit
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.constants import codata
import os

if plt.rcParams["text.usetex"] is False:
    plt.rcParams["text.usetex"] = True

if plt.rcParams["text.latex.unicode"] is False:
    plt.rcParams["text.latex.unicode"] = True


# Optischer Resonator
# Ein Spiegel ist eben, also r1 -> inf
def eben(L, r2):
    return 1-(L/r2)


# Optischer Resonator
# Beide Spiegel gekrümmt, Radien r1 und r2
def rund(L, r1, r2):
    return 1 - L/r1 - L/r2 + L**2/(r1*r2)


def stabilitaet():

    # Out Coupling ist fest
    r2 = 1.4  # in m

    Lplot = np.linspace(0, 3, 1000)
    # Planarer Spiegel
    plt.plot(Lplot, eben(Lplot, r2), 'k-', label=r'Flacher Spiegel $r_1 = \infty$')
    # Gekrümmter Spiegel
    r1 = 1.4  # in m
    plt.plot(Lplot, rund(Lplot, r1, r2), 'b-', label=r'$r_1 = 1,4\:\mathrm{m}$')
    # Gekrümmter Spiegel
    r1 = 1  # in m
    plt.plot(Lplot, rund(Lplot, r1, r2), 'r-', label=r'$r_1 = 1\:\mathrm{m}$')
    # Stabiler Bereich
    plt.barh(0.5 , 3, 1, 0, color="green", alpha=0.2, edgecolor="gray", label="Stabiler Bereich")
    plt.xlabel(r'$L\;/\;\mathrm{m}$')
    plt.ylabel(r'$g_1 \cdot g_2$')
    plt.xlim(0, 3)
    # plt.ylim(6.70, 16.95)
    plt.grid()
    plt.tight_layout()
    plt.legend(loc='best')
    # No whitespace around plots
    plt.savefig('build/stabilitaet.pdf', bbox_inches='tight')
    plt.clf()


if __name__ == '__main__':

    if not os.path.isdir('build'):
        os.mkdir('build')

    stabilitaet()
