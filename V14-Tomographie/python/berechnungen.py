import numpy as np
import numpy.linalg as lina
import pandas as pd
import matplotlib.pyplot as plt
import os

plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 14
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['xtick.labelsize'] = 13
plt.rcParams['ytick.labelsize'] = 13
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['legend.fontsize'] = 13
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

# make build directory if not existant
if not os.path.isdir('build'):
    os.mkdir('build')
#---------------------------
I_0_1 = 9520
t_I_1 = 60

I_0_2 = 6196
t_I_2 = 60

#   Wegstrecke
d = 0.01 #m

# --------Daten der Würfel(gemessen)
t_1, c_1, f_1 = np.genfromtxt('rohdaten/Würfel_1.txt', unpack=True)
t_2, c_2, f_2 = np.genfromtxt('rohdaten/Würfel_2.txt', unpack=True)
t_3, c_3, f_3 = np.genfromtxt('rohdaten/Würfel_3.txt', unpack=True)
t_4, c_4, f_4 = np.genfromtxt('rohdaten/Würfel_4.txt', unpack=True)


#--------- Plot zur Eingangsintensität
I0 = np.genfromtxt('rohdaten/kein_wuerfel.txt', unpack=True)
I0_error = np.sqrt(I0)

hist = np.linspace(0, len(I0), 512)
hist=hist+448
#print(len(I0), len(I0_error), len(hist))
plt.bar(hist, I0)
plt.xlim(33+448, 250+448)
plt.xlabel('Energie [keV]')
plt.ylabel('Ereignisse')
plt.savefig('build/Nullmessung1.pdf')
plt.clf()

#--------- Plot zur Eingangsintensität
I0 = np.genfromtxt('rohdaten/W4_I0.txt', unpack=True)
I0_error = np.sqrt(I0)

hist = np.linspace(0, len(I0), 512)
hist+=448
#print(len(I0), len(I0_error), len(hist))
plt.bar(hist, I0)
plt.xlim(33+448, 250+448)
plt.xlabel('Energie [keV]')
plt.ylabel('Ereignisse')
plt.savefig('build/Nullmessung2.pdf')
plt.clf()

#--------- Plot zum ersten Würfel W1_pos1
W1_pos1 = np.genfromtxt('rohdaten/W1_pos1.txt', unpack=True)
W1_pos1_error = np.sqrt(W1_pos1)

hist = np.linspace(0, len(W1_pos1), 512)
hist += 448
plt.bar(hist, W1_pos1)
plt.xlim(33+448, 250+448)
plt.xlabel('Energie [keV]')
plt.ylabel('Ereignisse')
plt.savefig('build/W1_pos1.pdf')
plt.clf()

##--------- Plot zum zweiften Würfel W2_pos1
W2_pos1 = np.genfromtxt('rohdaten/W2_pos1.txt', unpack=True)
W2_pos1_error = np.sqrt(W2_pos1)

hist = np.linspace(0, len(W2_pos1), 512)

plt.bar(hist, W2_pos1)
plt.xlim(33, 442)
plt.title('Messung des zweiten Würfels')
plt.xlabel('Kanal')
plt.ylabel('Ereignisse')
plt.savefig('build/W2_pos1.pdf')
plt.clf()

##--------- Plot zum zweiften Würfel W3_pos1
W3_pos1 = np.genfromtxt('rohdaten/W3_pos1.txt', unpack=True)
W3_pos1_error = np.sqrt(W3_pos1)
hist = np.linspace(0, len(W3_pos1), 512)

for i in range(1,len(hist)):
    hist[i]=hist[i]/(i*662)

plt.bar(hist, W3_pos1)
plt.xlim(33, 700)
plt.title('Messung des dritten Würfels')
plt.xlabel('Kanal')
plt.ylabel('Ereignisse')
plt.savefig('build/W3_pos1.pdf')
plt.clf()

# Initialisiere Würfelmatrix
i = np.sqrt(2)
A = np.array([[1, 0, 0, 1, 0, 0, 1, 0, 0],  # I1
              [0, 1, 0, 0, 1, 0, 0, 1, 0],
              [0, 0, 1, 0, 0, 1, 0, 0, 1],
              [1, 1, 1, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 1, 1, 1, 0, 0, 0],  # I5
              [0, 0, 0, 0, 0, 0, 1, 1, 1],
              [0, i, 0, i, 0, 0, 0, 0, 0],
              [0, 0, i, 0, i, 0, i, 0, 0],
              [0, 0, 0, 0, 0, i, 0, i, 0],
              [0, 0, 0, i, 0, 0, 0, i, 0],  # I10
              [i, 0, 0, 0, i, 0, 0, 0, i],
              [0, i, 0, 0, 0, i, 0, 0, 0]])

# Geometriematrix für Würfel 2 und 3 (da gleich durchstrahlt)
A_2 = np.matrix([[3],
                 [3],
                 [3*i],
                 [3*i]])

def kleinsteQuadrate(y, W, A):
     temp = np.dot(np.linalg.inv(np.dot(A.T, np.dot(W, A))), A.T)
     a = np.dot(temp, np.dot(W, y))
     a_err = np.linalg.inv(np.dot(A.T, np.dot(W, A)))
     return a, np.sqrt(np.diag(a_err))


#----------- Raten (counts/time)
# Der Nullmessungen
rate_null =  I_0_1/t_I_1
err_rate_null = np.sqrt(I_0_1)/t_I_1

rate_null_2 = I_0_2/t_I_2
err_rate_null_2 = np.sqrt(I_0_2)/t_I_2

# Des leeren Würfels
rate_leer = c_1/t_1
err_rate_leer = np.sqrt(c_1)/t_1

# Der drei Würfel
rate_2 = c_2/t_2
err_rate_2 = np.sqrt(c_2)/t_2

rate_3 = c_3/t_3
err_rate_3 = np.sqrt(c_3)/t_3

rate_4 = c_4/t_4
err_rate_4 = np.sqrt(c_4)/t_4


#------------Zuordnen der Projekstionen von den ersten Messungen zu der letzten
I_leer = np.array(np.zeros(len(rate_4)))
err_I_leer = np.array(np.zeros(len(rate_4)))

for i in range(len(rate_4)):
    if(i==1) or (i==4):
        I_leer[i] = rate_leer[0]
        err_I_leer[i] = err_rate_leer[0]
    if(i==0) or (i==2) or (i==3) or(i==5):
        I_leer[i] = rate_leer[1]
        err_I_leer[i] = err_rate_leer[1]
    if(i==6) or (i==8) or (i==9) or(i==11):
        I_leer[i] = err_I_leer[3]
        err_I_leer[i] = err_rate_leer[3]
    if(i==7) or (i==10):
        I_leer[i] = rate_leer[2]
        err_I_leer[i] = err_rate_leer[2]


# ln(I_0/N_j): Die j-ten Projektionen bei Würfel 4

I_2 = np.log(rate_leer/rate_2)
I_3 = np.log(rate_leer/rate_3)
I_4 = np.log(I_leer/rate_4)

err_I_2 = np.sqrt((np.sqrt(rate_2)/ rate_2)**2 + (err_rate_2/rate_2)**2)
err_I_3 = np.sqrt((np.sqrt(rate_3)/ rate_3)**2 + (err_rate_3/rate_3)**2)
err_I_4 = np.sqrt((np.sqrt(rate_4)/ rate_4)**2 + (err_I_leer/I_leer)**2)

print('''
    ###############################################################
    ~~~ Raten der verschiedenen Würfel ~~~
    -----------------------------------------------
    Werte = {}
    Fehler = {}
    Alu_leer mit Mittelung: (alle 12 Projektionen)
    -----------------------------------------------
    Werte = {}
    Fehler = {}
    Würfel 2: (Projektionen 2, 9, 10, 11)
    -----------------------------------------------
    Werte = {}
    Fehler = {}
    Würfel 3: (Projektionen 2, 10, 11)
    -----------------------------------------------
    Werte = {}
    Fehler = {}
    Würfel 4: (alle {} Projektionen)
    -----------------------------------------------
    Werte = {}
    Fehler = {}
    ###############################################################
    '''
      .format(rate_leer, err_rate_leer, I_leer, err_I_leer, rate_2, err_rate_2,
              rate_3, err_rate_3, len(rate_4), rate_4, err_rate_4))

print('''
~~~ Logarithmen der Raten der verschiedenen Würfel ~~~
Würfel 2:
-----------------------------------------------
Werte = {}
Fehler = {}
Würfel 3:
-----------------------------------------------
Werte = {}
Fehler = {}
Würfel 4:
-----------------------------------------------
Werte = {}
Fehler = {}
###############################################################
'''.format(I_2, err_I_2, I_3, err_I_3, I_4, err_I_4)
)

# Gewichtungsmatrizen
W_2 = np.diag(1/err_I_2**2)
W_3 = np.diag(1/err_I_3**2)
W_4 = np.diag(1/err_I_4**2)

# Über kleinsteQuadrate mu berechnen
mu_2, err_mu_2 = kleinsteQuadrate(I_2, W=W_2, A=A_2)
mu_3, err_mu_3 = kleinsteQuadrate(I_3, W=W_3, A=A_2)
mu_4, err_mu_4 = kleinsteQuadrate(I_4, W=W_4, A=A)

print(
'''
~~~ Absorptionskoeffizienten der verschiedenen Würfel ~~~
Würfel 2:
-----------------------------------------------
Werte = {}
Fehler = {}
Würfel 3:
-----------------------------------------------
Werte = {}
Fehler = {}
Würfel 4:
-----------------------------------------------
Werte = {}
Fehler = {}
'''.format(mu_2, err_mu_2, mu_3, err_mu_3, mu_4, err_mu_4)
)
