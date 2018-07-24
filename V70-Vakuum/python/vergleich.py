import numpy as np
import matplotlib.pyplot as plt
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import uncertainties.unumpy as unp

print("Drehschieberpumpe")

druckmessfehler = 0.2
# Aus Leckratenmessung
p_leck = np.array([0.1, 0.4, 0.8, 1.0])
p_leck = unp.uarray(p_leck, p_leck * druckmessfehler)
S_leck = unp.uarray([0.47, 0.80, 1.26, 1.31],
                    [0.1, 0.17, 0.27, 0.28])
for i in range(len(noms(p_leck))):
    plt.errorbar(x=noms(p_leck[i]),
                 xerr=stds(p_leck[i]),
                 y=noms(S_leck[i]),
                 yerr=stds(S_leck[i]),
                 label=r'Leck {}'.format(i+1))
# Aus Evakuierungskurve
# Bereich 1 zwischen 1000 und 40 Milli Bar
plt.errorbar(x=520,
             xerr=480,
             y=1.06,
             yerr=0.23,
             label=r'Evak 1')
# Bereich 2 zwischen 20 und 0,8 Milli Bar
plt.errorbar(x=10.4,
             xerr=9.6,
             y=1.08,
             yerr=0.08,
             label=r'Evak 2')
# Bereich 3 zwischen 0,6 und 0,06 Milli Bar
plt.errorbar(x=0.33,
             xerr=0.27,
             y=0.52,
             yerr=0.05,
             label=r'Evak 3')
# Aus Anleitung
plt.axhline(y=1.1, label='Hersteller')
plt.xlabel(r'$p\;/\;\mathrm{mbar}$')
plt.ylabel(r'$S\;/\;\mathrm{L\;s}^{-1}$')
# plt.xlim(60, 310)
# plt.ylim(6.70, 16.95)
plt.xscale('log')
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
# No whitespace around plots
plt.savefig('build/vergleich_D.pdf', bbox_inches='tight')
plt.clf()


print("Turbomolekularpumpe")

druckmessfehler = 0.1
# Aus Leckratenmessung
p_leck = np.array([5, 10, 15, 20])*1e-5
p_leck = unp.uarray(p_leck, p_leck * druckmessfehler)
S_leck = unp.uarray([14.6, 23.6, 23.7, 25.0],
                    [1.8, 3.0, 3.0, 3.2])
for i in range(len(noms(p_leck))):
    plt.errorbar(x=noms(p_leck[i]),
                 xerr=stds(p_leck[i]),
                 y=noms(S_leck[i]),
                 yerr=stds(S_leck[i]),
                 label=r'Leck {}'.format(i+1))
# Aus Evakuierungskurve
# Bereich 1 zwischen 5e-3 und 2e-4 Milli Bar
plt.errorbar(x=2.6e-3,
             xerr=2.4e-3,
             y=9.2,
             yerr=0.9,
             label=r'Evak 1')
# Bereich 2 zwischen 3e-5 und 8e-5 Milli Bar
plt.errorbar(x=5.5e-5,
             xerr=2.5e-5,
             y=3.3,
             yerr=0.5,
             label=r'Evak 2')
# Aus Anleitung
plt.axhline(y=77, label='Hersteller')
plt.xlabel(r'$p\;/\;\mathrm{mbar}$')
plt.ylabel(r'$S\;/\;\mathrm{L\;s}^{-1}$')
# plt.xlim(60, 310)
# plt.ylim(6.70, 16.95)
plt.xscale('log')
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
# No whitespace around plots
plt.savefig('build/vergleich_T.pdf', bbox_inches='tight')
plt.clf()
