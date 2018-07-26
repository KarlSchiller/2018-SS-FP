import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.optimize import curve_fit
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.constants import codata

def ber_lambda(d, g, L, n):
	return g*np.sin(np.arctan(d/L)/n)


o, d = np.genfromtxt("rohdaten/wellenlaenge.txt", unpack=True)
d /= 100
lamb = []
for n in range(len(d)):
	if(n>4):
		m = (n-4)*(-1)
	else:
		m=(n+1)
	print(n, m,  d[n])
	lamb.append(ber_lambda(d[n], 1e-3/1000, 0.113, m))

for i in range(len(lamb)):
	lamb[i] = lamb[i]*10**9

print(lamb)
print(np.mean(lamb[1:5]))
print(np.mean(lamb[6:-1]))

ascii.write([o, d, lamb],
            'table/wellenlaenge.tex',
            format='latex',
            overwrite=True)

lambda_theo = 632.8e-9
delta_lambda = (np.mean(lamb)/lambda-1)*100
print('Abweichung: ', delta_lambda)
