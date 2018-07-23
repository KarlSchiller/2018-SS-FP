import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.optimize import curve_fit
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.constants import codata


def polar(x, I0, phi0):
    return I0*np.sin(x+phi0)**2


def evak_T():
    phi, I = np.genfromtxt("rohdaten/polarisation.txt", unpack=True)
    phi = phi/360*2*np.pi
    x = np.linspace(0, 2*np.pi)

    params, covariance = curve_fit(f=polar, xdata=phi, ydata=I)
    errors = np.sqrt(np.diag(covariance))
    print("I0 : ",params[0])
    print("phi0 : ", params[1])
    # Plot
    plt.plot(phi, I, 'kx', label='Messwerte')
    plt.plot(x, polar(x, *params), label='Fit')
    plt.xlabel(r'$\phi/rad$')
    plt.ylabel(r'$I/nA$')
    plt.xlim(0, 2*np.pi)
    # plt.ylim(6.70, 16.95)
    plt.legend(loc='best')
    plt.grid()
    plt.tight_layout()
    plt.savefig('build/polarisation.pdf')
    plt.clf()



if __name__ == '__main__':
    evak_T()
