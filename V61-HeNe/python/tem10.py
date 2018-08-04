import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.optimize import curve_fit
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.constants import codata


def t10(d, I1, d1, w1, I2, d2, w2):
    return I1 * np.exp( -2 * ( (d-d1) / w1 )**2 )+ I2 * np.exp( -2 * ( ( d-d2) / w2 )**2 )


def evak_T():
    d, I = np.genfromtxt("rohdaten/T10.txt", unpack=True)
    x = np.linspace(-30, 30)

    Imaxleft=max(I[1:15])
    Imaxright=max(I[17:-1])
    rmaxleft=d[np.where(I==Imaxleft)[0]]
    rmaxright=d[np.where(I==Imaxright)[0]]

    print('I_1 ',Imaxleft)
    print('d_1 ', rmaxleft)
    print('I_2 ',Imaxright)
    print('d_2 ', rmaxright)

    params, covariance = curve_fit(t10, d, I, p0=[Imaxleft, -10, 1, Imaxright, 12, 1])
    errors = np.sqrt(np.diag(covariance))

    I1 = ufloat(params[0], errors[0])
    d1 = ufloat(params[1], errors[1])
    w1 = ufloat(params[2], errors[2])
    I2 = ufloat(params[3], errors[3])
    d2 = ufloat(params[4], errors[4])
    w2 = ufloat(params[5], errors[5])

    print("I1 : ", I1)
    print("d1 : ", d1)
    print("w1 : ", w1)
    print("I2 : ", I2)
    print("d2 : ", d2)
    print("w2 : ", w2)

    # Plot
    plt.plot(d, I, 'kx', label='Messwerte')
    plt.plot(x, t10(x, *params), label='Fit')
    plt.xlabel(r'$d/\ mm$')
    plt.ylabel(r'$I/\ nA$')
    # plt.xlim(60, 310)
    # plt.ylim(6.70, 16.95)
    plt.legend(loc='best')
    plt.grid()
    plt.tight_layout()
    plt.savefig('build/T10.pdf')
    plt.clf()

    #ascii.write([d[1:15], I[1:15], d[16:-1], I[16:-1]],
    #            'table/T_10.tex',
    #            format='latex',
    #            overwrite=True)

if __name__ == '__main__':

    evak_T()
