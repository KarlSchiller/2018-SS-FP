import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.optimize import curve_fit
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.constants import codata


def t00(d, I0, d0, w):
    return I0*np.exp( -2*( (d-d0)/w )**2 )


def evak_T():
    d, I = np.genfromtxt("rohdaten/T00.txt", unpack=True)

    x = np.linspace(-30, 30)
    params, covariance = curve_fit(t00, d, I)
    errors = np.sqrt(np.diag(covariance))
    print("d0 : ", params[0])
    print("I0 : ", params[1])
    print("w : ", params[2])

    # Plot
    plt.plot(d, I, 'kx', label='Messwerte')
    plt.plot(x, t00(x, *params), label='Fit')
    plt.xlabel(r'$d/mm$')
    plt.ylabel(r'$I/nA$')
    # plt.xlim(60, 310)
    # plt.ylim(6.70, 16.95)
    plt.legend(loc='best')
    plt.grid()
    plt.tight_layout()
    plt.savefig('build/T00.pdf')
    plt.clf()




if __name__ == '__main__':

    evak_T()
