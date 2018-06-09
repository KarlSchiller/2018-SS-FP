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


# Volumina der Rezipienten
volumes = np.genfromtxt("build/volumina.txt", unpack=True)

# Evakuierungskurve der Drehschieberpumpe
print("Evakuierungskurve Drehschieberpumpe\n")

# Einlesen und bearbeiten der Messwerte
V_D = ufloat(volumes[0], volumes[1])
print('Vol D     {:} m^3'.format(V_D))
p_D, t1_D, t2_D, t3_D, t4_D, t5_D, t6_D = np.genfromtxt("rohdaten/druck_dreh.txt", unpack=True)
druckmessfehler_D = 0.2  # Fehler des Druckmessers
pend_D = ufloat(3e-2, druckmessfehler_D * 3e-2)  # Enddruck in mbar
pstart_D = ufloat(1000, druckmessfehler_D * 1000)  # Startdruck in mbar
print("Startdruck   {:} mbar".format(pstart_D))
print("Enddruck     {:} mbar\n".format(pend_D))
p_D = unp.uarray(p_D, druckmessfehler_D * p_D)
lnp_D = unp.log((p_D - pend_D)/(pstart_D - pend_D))
times_D = np.array([t1_D, t2_D, t3_D, t4_D, t5_D, t6_D])
times_D = np.cumsum(times_D, axis=1)  # wirkliche Zeiten sind die Kumulante
tmean_D = unp.uarray(np.mean(times_D, axis=0), np.std(times_D, axis=0, ddof=1))

# Fitten an einzelne Druckbereiche
# Bereich 1: Messwerte 2-4
print('Regression Bereich 1')
params1_D, covariance1_D = curve_fit(f=linear, xdata=noms(tmean_D[1:4]), ydata=noms(lnp_D[1:4]))
# covariance is the covaniance matrix
errors1_D = np.sqrt(np.diag(covariance1_D))
print('m = ', params1_D[0], ' +/- ', errors1_D[0], ' 1/s')
print('b = ', params1_D[1], ' +/- ', errors1_D[1])
m1_D = ufloat(params1_D[0], errors1_D[0])
S1_D = -m1_D * V_D
print('Saugvermögen 1       {:} m^3/s\n'.format(S1_D))

# Bereich 2: Messwerte 5-12
print('Regression Bereich 2')
params2_D, covariance2_D = curve_fit(f=linear, xdata=noms(tmean_D[4:12]), ydata=noms(lnp_D[4:12]))
# covariance is the covaniance matrix
errors2_D = np.sqrt(np.diag(covariance2_D))
print('m = ', params2_D[0], ' +/- ', errors2_D[0], ' 1/s')
print('b = ', params2_D[1], ' +/- ', errors2_D[1])
m2_D = ufloat(params2_D[0], errors2_D[0])
S2_D = -m2_D * V_D
print('Saugvermögen 2       {:} m^3/s\n'.format(S2_D))

# Bereich 3: Messwerte 13-16
print('Regression Bereich 3')
params3_D, covariance3_D = curve_fit(f=linear, xdata=noms(tmean_D[12:17]), ydata=noms(lnp_D[12:17]))
# covariance is the covaniance matrix
errors3_D = np.sqrt(np.diag(covariance3_D))
print('m = ', params3_D[0], ' +/- ', errors3_D[0], ' 1/s')
print('b = ', params3_D[1], ' +/- ', errors3_D[1])
m3_D = ufloat(params3_D[0], errors3_D[0])
S3_D = -m3_D * V_D
print('Saugvermögen 3       {:} m^3/s\n'.format(S3_D))

# Plot
tplot1_D = np.linspace(noms(tmean_D[1]), noms(tmean_D[3]), 10)
plt.plot(tplot1_D, linear(tplot1_D, *params1_D), 'b-', label='Regression 1')
tplot2_D = np.linspace(noms(tmean_D[4]), noms(tmean_D[11]), 10)
plt.plot(tplot2_D, linear(tplot2_D, *params2_D), 'g-', label='Regression 2')
tplot3_D = np.linspace(noms(tmean_D[12]), noms(tmean_D[16]), 10)
plt.plot(tplot3_D, linear(tplot3_D, *params3_D), 'r-', label='Regression 3')
plt.errorbar(x=noms(tmean_D), xerr=stds(tmean_D),
                y=noms(lnp_D), yerr=stds(lnp_D),
                color='k', fmt='x', label='Messpunkte')
plt.xlabel(r'$t_\mathrm{m}\;\mathrm{in}\;\mathrm{s}$')
plt.ylabel(r'$\ln(\frac{p(t)-p_\mathrm{e}}{p_\mathrm{0}-p_\mathrm{e}})$')
# plt.xlim(60, 310)
# plt.ylim(6.70, 16.95)
# plt.grid()
plt.tight_layout()
plt.savefig('build/evak_D.pdf')
plt.clf()


# Evakuierungskurve der Turbomolekularpumpe
print("Evakuierungskurve Turbomolekularpumpe\n")

V_T = ufloat(volumes[2], volumes[3])
print('Vol T     {:} m^3'.format(V_T))
