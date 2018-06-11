import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.optimize import curve_fit
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import uncertainties.unumpy as unp
from uncertainties import ufloat


def linear(x, m, b):
    return m*x + b


# Auswertung einer Leckratenmessung
# Erstellt an eine lineare Regression an die gemittelten Zeiten und
# berechnet aus der Steigung das Saugvermögen. Speichert die Ergebnisse
# als Plot und Tabelle
# INPUT:    Volumen V des Rezipienten
#           Unsicherheit p_err des Druckmessers
#           String datafile mit dem Pfad zu der auszuwertenden Datei
#           Gleichgewichtsdruck pg in mbar
#           String plotfile mit dem Pfad des zu speichernden Plots
#           String tablefile mit dem Pfad der zu speichernden Tabelle
def leckrate(datafile, V=10, p_err=0.1, pg=1, plotfile='plot.pdf', tablefile='table.tex'):
    p, t1, t2, t3 = np.genfromtxt(fname=datafile, unpack=True)
    druckmessfehler = p_err  # Fehler des Druckmessers
    p = unp.uarray(p, druckmessfehler * p)
    pg = ufloat(pg, druckmessfehler * pg)
    times = np.array([t1, t2, t3])
    times = np.cumsum(times, axis=1)  # wirkliche Zeiten sind die Kumulante
    tmean = unp.uarray(np.mean(times, axis=0), np.std(times, axis=0, ddof=1))

    print('Regression Gleichgewichtsdruck {:} mbar'.format(pg))
    params, covariance = curve_fit(f=linear, xdata=noms(tmean), ydata=noms(p))
    # covariance is the covaniance matrix
    errors = np.sqrt(np.diag(covariance))
    print('m = ', params[0], ' +/- ', errors[0], ' mbar/s')
    print('b = ', params[1], ' +/- ', errors[1], ' mbar')
    m = ufloat(params[0], errors[0])

    # Berechne Saugvermögen
    S = m * V_D / pg
    print('Saugvermögen       {:} m^3/s\n'.format(S))

    # Plotte Messwerte und Regression
    tplot = np.linspace(noms(tmean[0]), noms(tmean[-1]), 10)
    plt.plot(tplot, linear(tplot, *params), 'b-', label='Lineare Regression')
    plt.errorbar(x=noms(tmean), xerr=stds(tmean),
                    y=noms(p), yerr=stds(p),
                    color='k', fmt='x', label='Messpunkte')
    plt.xlabel(r'$\overline{t_\mathrm{1..3}}\;/\;\mathrm{s}$')
    plt.ylabel(r'$p\;/\;\mathrm{mbar}$')
    # plt.xlim(60, 310)
    # plt.ylim(6.70, 16.95)
    plt.grid()
    plt.tight_layout()
    # No whitespace around plots
    plt.savefig(plotfile, bbox_inches='tight')
    plt.clf()

    # Erstelle Tabelle
    ascii.write(
                [noms(p),
                stds(p),
                np.round(times[0], 2),
                np.round(times[1], 2),
                np.round(times[2], 2),
                np.round(noms(tmean), 2),
                np.round(stds(tmean), 2)],
                output=tablefile,
                format='latex',
                overwrite=True)


# Auswertung der Leckratenmessung der Drehschieberpumpe
def Drehschieberpumpe():
    print("Leckratenmessung Drehschieberpumpe\n")
    druckmessfehler_D = 0.2  # Pirani Messgerät

    leckrate(datafile="rohdaten/leck_dreh_0,1mbar.txt",
             V=V_D,
             p_err=druckmessfehler_D,
             pg=0.1,
             plotfile='build/leck/D_0,1.pdf',
             tablefile='build/tab/D_0,1.tex')

    leckrate(datafile="rohdaten/leck_dreh_0,4mbar.txt",
             V=V_D,
             p_err=druckmessfehler_D,
             pg=0.4,
             plotfile='build/leck/D_0,4.pdf',
             tablefile='build/tab/D_0,4.tex')

    leckrate(datafile="rohdaten/leck_dreh_0,8mbar.txt",
             V=V_D,
             p_err=druckmessfehler_D,
             pg=0.8,
             plotfile='build/leck/D_0,8.pdf',
             tablefile='build/tab/D_0,8.tex')

    leckrate(datafile="rohdaten/leck_dreh_1,0mbar.txt",
             V=V_D,
             p_err=druckmessfehler_D,
             pg=1.0,
             plotfile='build/leck/D_1,0.pdf',
             tablefile='build/tab/D_1,0.tex')


# Auswertung der Leckratenmessung der Turbomolekularpumpe
def Turbomolekularpumpe():
    print("Leckratenmessung Turbomolekularpumpe\n")
    druckmessfehler_T = 0.1  # Glühkathoden Messgerät

    leckrate(datafile="rohdaten/leck_turbo_5e-5mbar.txt",
             V=V_T,
             p_err=druckmessfehler_T,
             pg=5e-5,
             plotfile='build/leck/T_5e-5.pdf',
             tablefile='build/tab/T_5e-5.tex')

    leckrate(datafile="rohdaten/leck_turbo_1e-4mbar.txt",
             V=V_T,
             p_err=druckmessfehler_T,
             pg=1e-4,
             plotfile='build/leck/T_1e-4.pdf',
             tablefile='build/tab/T_1e-4.tex')

    leckrate(datafile="rohdaten/leck_turbo_1,5e-4mbar.txt",
             V=V_T,
             p_err=druckmessfehler_T,
             pg=1.5e-4,
             plotfile='build/leck/T_1,5e-4.pdf',
             tablefile='build/tab/T_1,5e-4.tex')

    leckrate(datafile="rohdaten/leck_turbo_2e-4mbar.txt",
             V=V_T,
             p_err=druckmessfehler_T,
             pg=2e-4,
             plotfile='build/leck/T_2e-4.pdf',
             tablefile='build/tab/T_2e-4.tex')


if __name__ == '__main__':

    # Volumina der Rezipienten
    volumes = np.genfromtxt("build/volumina.txt", unpack=True)
    V_D = ufloat(volumes[0], volumes[1])
    print('Vol D     {:} m^3'.format(V_D))
    V_T = ufloat(volumes[2], volumes[3])
    print('Vol D     {:} m^3\n'.format(V_T))

    if not os.path.isdir('build'):
        os.mkdir('build')
    if not os.path.isdir('build/leck'):
        os.mkdir('build/leck')
    if not os.path.isdir('build/tab'):
        os.mkdir('build/tab')

    # Leckratenmessung der Drehschieberpumpe
    Drehschieberpumpe()
    # Leckratenmessung der Turbomolekularpumpe
    Turbomolekularpumpe()
