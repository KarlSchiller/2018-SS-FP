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
    phi_deg = phi
    phi_rad = np.deg2rad(phi)
    x = np.linspace(0, 2*np.pi)

    params, covariance = curve_fit(f=polar, xdata=phi, ydata=I)
    errors = np.sqrt(np.diag(covariance))

    I0 = ufloat(params[0], errors[0])
    phi0 = ufloat(params[1], errors[1])

    print("I0 : ",I0)
    print("phi0 : ", phi0)
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

    ascii.write([phi_deg[1:18], np.round(phi_rad[1:18], 2), I[1:18], phi_deg[19:-1], np.round(phi_rad[1:18], 2), I[1:18]],
            'table/polar.tex',
            format='latex',
            overwrite=True
            )

if __name__ == '__main__':
    evak_T()
