import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.optimize import curve_fit
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.constants import codata


def linear(x, m, b):
    return m*x + b


def evak_T():
    phi, I = np.genfromtxt("rohdaten/polarisation.txt", unpack=True)
    # Plot
    plt.plot(phi, I, 'kx', label='Messwerte')
    # plt.xlabel(r'$\overline{t_\mathrm{1..6}}\;/\;\mathrm{s}$')
    # plt.ylabel(r'$\ln(\frac{p(t)-p_\mathrm{e}}{p_\mathrm{0}-p_\mathrm{e}})$')
    # plt.xlim(60, 310)
    # plt.ylim(6.70, 16.95)
    plt.grid()
    plt.tight_layout()
    plt.savefig('build/polarisation.pdf')
    plt.clf()




if __name__ == '__main__':

    evak_T()
