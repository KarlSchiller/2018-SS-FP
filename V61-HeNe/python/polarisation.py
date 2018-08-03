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
    x = np.linspace(phi_rad[0], phi_rad[-1])

    params, covariance = curve_fit(f=polar, xdata=phi_rad, ydata=I)
    errors = np.sqrt(np.diag(covariance))

    I0 = ufloat(params[0], errors[0])
    phi0 = ufloat(params[1], errors[1])
    phi0_2 = ufloat(np.rad2deg(params[1]), errors[1])

    print("I0 : ",I0)
    print("phi0_rad : ", phi0)
    print("phi0_deg: ", phi0_2/360)

    # Plot
    plt.plot(phi_rad, I, 'kx', label='Messwerte')
    plt.plot(x, polar(x, *params), label='Fit')
    plt.xlabel(r'$\phi/rad$')
    plt.ylabel(r'$I/nA$')
    # plt.yscale(0, 2*np.pi)
    # plt.ylim(6.70, 16.95)
    plt.legend(loc='best')
    plt.xticks([0,0.25*np.pi,0.5*np.pi,0.75*np.pi,np.pi,1.25*np.pi,1.5*np.pi,1.75*np.pi, 2*np.pi],['0','$\\frac{1}{4}\,\\pi$', '$\\frac{1}{2}\,\\pi$','$\\frac{3}{4}\,\\pi$' ,'$\\pi$','$\\frac{5}{4}\,\\pi$','$\\frac{3}{2}\, \\pi$','$\\frac{7}{4}\, \\pi$', '$2\, \\pi$'])
    plt.grid()
    plt.tight_layout()
    plt.savefig('build/polarisation.pdf')
    plt.clf()

    #ascii.write([phi_deg[1:18], np.round(phi_rad[1:18], 2), I[1:18], phi_deg[19:-1], np.round(phi_rad[1:18], 2), I[1:18]],
    #        'table/polar.tex',
    #        format='latex',
    #        overwrite=True
    #        )

if __name__ == '__main__':
    evak_T()
