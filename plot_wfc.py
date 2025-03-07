import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

colors=sns.color_palette('tab20b')

energies = []
wavefunctions = []
coordinates = []
hs = []

zeta = 1.6875  # STO parameter for He 1s

gto3_exp = [
    0.6362421394e+01,
    0.1158922999e+01,
    0.3136497915e+00]

gto3_coef = [
    0.1543289673e+00,
    0.5353281423e+00,
    0.4446345422e+00]

gto4_exp = [
    0.1489982967E+02,
    0.2726485258E+01,
    0.7574474599E+00,
    0.2513900027E+00]

gto4_coef = [
    0.5675242080E-01,
    0.2601413550E+00,
    0.5328461143E+00,
    0.2916254405E+00]

def normalize(arr, h):
    nor = np.sum(arr**2) * h
    return 1/(np.sqrt(nor*4*np.pi)) * arr

def sstate(r, coefs, exps):
    result = [] 
    for i in range(len(coefs)):
        coef = coefs[i]
        ex = exps[i]
        result.append(coef * (2*ex/(np.pi)) ** (3/4) * np.exp(-ex*(r.copy()**2)))
    result = np.array(result)
    return np.sum(result, axis=0)

def snorm(state, coor, h):
    summe = 0
    for i in range(len(state)-1):
        summe += (((state[i+1] + state[i]) * 0.5 ) ** 2 * h) * 4 * np.pi * ((coor[i+1]+coor[i])/2)**2
    state = (1/np.sqrt(summe)) * state
    return state

with open('helium.txt', 'r') as f:
    for line in f:
        h,t,e,wave,coords = line.split('\t')
        energies.append(e)
        hs.append(h)
        wavefunctions.append(list(map(float, wave.split())))
        coordinates.append(list(map(float, coords.split())))

testfunc = gto3_coef[2] * (2*gto3_exp[2]/(np.pi)) ** (3/4) * np.exp(-gto3_exp[2]*(np.array(coordinates[-1])**2))
testfunc1 = gto3_coef[1] * (2*gto3_exp[1]/(np.pi)) ** (3/4) * np.exp(-gto3_exp[1]*(np.array(coordinates[-1])**2))
testfunc0 = gto3_coef[0] * (2*gto3_exp[0]/(np.pi)) ** (3/4) * np.exp(-gto3_exp[0]*(np.array(coordinates[-1])**2))

sto3g = sstate(np.array(coordinates[-1]), gto3_coef, gto3_exp)
sto4g = sstate(np.array(coordinates[-1]), gto4_coef, gto4_exp)
sto4gu = sstate(np.array(coordinates[-1]), gto4_coef, gto4_exp)
basisbased = sstate(np.array(coordinates[-1]), gto3_coef, gto3_exp)
basisbased = snorm(basisbased, np.array(coordinates[-1]), float(hs[-1]))
sto4gu = snorm(sto4gu, np.array(coordinates[-1]), float(hs[-1]))
basisbased = basisbased * (np.array(coordinates[-1]))
sto4gu = sto4gu * (np.array(coordinates[-1]))

sto = (zeta**3 / np.pi)**0.5 * np.exp(-zeta * np.array(coordinates[-1]))
wavefunction = (zeta**3 / np.pi)**0.5 * np.exp(-zeta * np.array(coordinates[-1]))
wavefunction = snorm(wavefunction, np.array(coordinates[-1]), float(hs[-1]))
wavefunction = wavefunction * np.array(coordinates[-1])

fullphi = np.array(wavefunctions[-1][100:])/np.array(coordinates[-1][100:])


fig, axs = plt.subplots(1,2, figsize=(6,3))
for i in range(1):
    wavefunctions[-1]=np.array(wavefunctions[-1])
    axs[1].plot(coordinates[-1], (4 * np.pi) * wavefunctions[-1]**2, label=f"$\Delta_r =0.0001$", color=colors[1])
    axs[1].set_xlim(-0.1,5)
axs[1].plot(coordinates[-1], np.pi * 4 * basisbased ** 2, label='STO-3G', color=colors[5])
#axs[1].plot(coordinates[-1], np.pi * 4 * wavefunction ** 2, label='STO', color=colors[9])
axs[1].set_xlabel('$r$ (a.u.)')
axs[1].set_ylabel('$4\pi r^2 |\psi(r)|^2$')
axs[1].legend(fontsize='small')
axs[0].plot(coordinates[-1][100:], fullphi, label=f"$\Delta_r =0.0001$", color=colors[1])
axs[0].plot(coordinates[-1], sto3g, label='STO-3G', color=colors[5])
#axs[0].plot(coordinates[-1], sto, label='STO', color=colors[9])
axs[0].set_xlabel('$r$ (a.u.)')
axs[0].set_ylabel('$\psi(r)$')
axs[0].legend(fontsize='small')
axs[0].set_xlim(-0.1,5)
plt.subplots_adjust(wspace=0.35)
plt.subplots_adjust(right=0.99)
plt.subplots_adjust(left=0.1)
plt.subplots_adjust(bottom=0.15)
plt.subplots_adjust(top = 0.99)
plt.savefig('helium.pdf')

