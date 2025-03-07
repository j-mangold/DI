import numpy as np
import matplotlib.pyplot as plt
import time

def hamiltonian(V, xs, dx):
    diag = V + dx ** (-2)
    off = -np.ones(len(xs)-1) * 0.5 * dx ** (-2)
    return diag, off

def get_bound_states(V, xs, dx, n=1):
    main, off = hamiltonian(V, xs, dx)
    E, psi = eigh_tridiagonal(main, off, select='i', select_range=(0, n-1))
    if n==1:
        psi = psi.ravel()
    return E, psi.T

def downscale(psi, r, h):
    start = r[0]
    end = r[-1]
    hfine = r[1]-r[0]
    step = int(h/hfine)
    indices = np.arange(0, len(psi), step).astype(int)
    rnew = np.arange(start, end, h)
    downscaled = psi[indices]
    return downscaled, rnew


def corepot(r):
    return - 2/r

def boundarymaker(fun, coordinates):
    end = fun[-1]
    endcoordinate = coordinates[-1]
    return fun - (end-1)/endcoordinate * coordinates

def normalize(arr, h):
    nor = np.sum(arr**2) * h
    return 1/(np.sqrt(nor*4*np.pi)) * arr

def poisson(u, epsilon):
    inpot = np.zeros(len(coord))
    u = normalize(u, h)
    a = np.zeros(len(coord))
    a[1:] = -4 * np.pi  * u[1:] ** 2 / (coord[1:])
    a[0] = 0

    q = np.zeros(len(coord))

    q[0] = 0
    q[1] = h

    for j in ind[2:]:
        q[j] = 2*q[j-1] - q[j-2] + a[j] * h**2 

    q = boundarymaker(q, coord)
    
    inpot[1:] = q[1:] / coord[1:]
    inpot[0] = inpot[1]

    Eee = 4*np.pi*np.sum((u**2) * inpot)*h

    print(f"Eee: {Eee}, iteration: {iteration}, epsilon: {epsilon}")

    return inpot, Eee

hs = np.logspace(-2,-6,5)
times = []
energies = []
wavefunctions = []
coordinates = []

for h in hs:

    start = time.time()
    Z = 2
    up = 40 
    low = 0
    deltae = 100
    cutoff = 1e-5
    
    eps = - Z ** 2 * 0.5
    epsilons = []
    
    coord = np.arange(low, up, h)
    ind = np.arange(0, len(coord), 1).astype(int)
    iteration = 0
    
    
    pot = np.zeros(len(coord))
    pot[1:] = corepot(coord[1:])
    pot[0] = pot[1]
    E, psi = get_bound_states(pot, coord, h)
    
    epsilons.append(E[0])
    
    pot2, Eee = poisson(psi, E)
    
    while deltae > cutoff and iteration < 100:
    
        iteration += 1
    
        E, psi = get_bound_states(pot + pot2, coord, h)
        
        epsilons.append(E[0])
    
        deltae = abs(epsilons[-1] - epsilons[-2])
        print(f"deltaE: {deltae}")
        print(np.sum(psi **2))
    
        pot2, Eee = poisson(psi, E)
    
    print(f"final result: {2 * E - Eee}")

    psi = normalize(psi, h) 
    psi, coord = downscale(psi, coord, hs[0])
    energies.append(2 * E - Eee)
    wavefunctions.append(psi)
    coordinates.append(coord)
    end = time.time()
    times.append(end-start)

with open('helium.txt', 'w') as f:
    for i in range(len(hs)):
        f.write(str(hs[i]))
        f.write('\t')
        f.write(str(times[i]))
        f.write('\t')
        f.write(str(energies[i][0]))
        f.write('\t')
        f.write(" ".join(map(str, wavefunctions[i])))
        f.write('\t')
        f.write(" ".join(map(str, coordinates[i])))
        f.write('\n')
