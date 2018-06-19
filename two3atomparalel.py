import pylab
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import scipy
from qutip import *
from scipy import linalg
plt.style.use('ggplot')
from two3levelHamiltonian import *
from nndist import *

#Code about two three-level atoms

def ind(i,j, dim):
    #map of i and j into a single index, runs through columns first, and rows second
    return i*dim + j

t = np.linspace(0., 400., 2000.)  #vector for time evolution

Delta_1 = 0.0001
#Delta_2 = 0
Omeg_1 = 0.2
#coupling from g to e
Omeg_2 = 1.    #coupling from e to r
gamma_eg = 0.4  #decay from e to g
gamma_re = 0.01  #decay from r to e
gamma_rg = 0.0  #decay r to g
gamma_egcol = 0.0  #collective decay, e to g
V = 1000000. #constant characterizing the dipole interaction strength

rho = 1/1000. #density of the medium

dim = 3 #three-level system
idx = range(dim)  #vector, goes from to 0 to dim-1

kets_1atom = qutrit_basis(); #basis for a three-level system
bras_1atom = []
for i in idx:
    bras_1atom.append((kets_1atom[i]).dag()) #dual basis for one atom
     
#Sigma matrices for one atom only.
sig_1atom = []
for i in idx:
    for j in idx:
        sig_1atom.append(kets_1atom[i]*bras_1atom[j])

#Sigma matrices for the expanded space of 2 atoms
#The indices go from 1 to dim*dim. sig1[0] = sig1_gg; sig1[1] = sig1_ge, etc.
sig1 = []
sig2 = []
for sig_aux in sig_1atom:
    sig1.append(tensor(sig_aux, identity(dim)))
    sig2.append(tensor(identity(dim), sig_aux))


trace = []

pop_r = []
pop_rr = []

coherence = []

#operators whose expectation values are of interest
popinr = sig1[ind(2,2,dim)] + sig2[ind(2,2,dim)]
popine = sig1[ind(1,1,dim)] + sig2[ind(1,1,dim)]
poping = sig1[ind(0,0,dim)] + sig2[ind(0,0,dim)]
    
obs = sig1[ind(1,0,dim)] + sig2[ind(1,0,dim)]
    
popinrr = sig1[ind(2,2,dim)]*sig2[ind(2,2,dim)]

num_step = 100
step_x = 1000
Delta2vec = np.linspace(-10.,10., num_step)
xvec = np.linspace(0.01, 25, step_x)
dx = (25. - 0.01)/(1.*step_x)



paveraged_e = []
paveraged_r = []  
for Delta_2 in Delta2vec:
    pop_e = []
    pop_r = []
    for x in xvec:
        Vdip = V/(x**6.)
        H = Hamiltonian(Delta_1, Delta_2, Omeg_1, Omeg_2, Vdip)
        #Heff = Hamiltonianeff(Delta_1, Delta_2, Omeg_1, Omeg_2, Vdip)
        decay_ops = decay(gamma_eg, gamma_re, gamma_rg, gamma_egcol)

        psi0 = tensor(kets_1atom[0], kets_1atom[0]) #initial state

        rho_ss = steadystate(H, decay_ops, method = 'svd'); #steady state, given H and decay_ops

        pop_e.append((rho_ss*popine).tr())  #population in e, as a function of x, for a given Delta_2
        pop_r.append((rho_ss*popinr).tr())

    #Now calculate integral of pop_e (x)  *nndist (x) delta x
    
    paveraged_e.append(dx*np.inner(np.asarray(pop_e), nndist(xvec, rho)))
    paveraged_r.append(dx*np.inner(np.asarray(pop_r), nndist(xvec, rho)))

vector_ones = np.ones(len(nndist(xvec, rho)))    
        

    
#result = mesolve(H, psi0, t, decay_ops, [exp1, exp2,exp3,exp4])

#plt.figure()
#plt.plot(Delta2vec, trace, '-', color = 'blue')

plt.figure()
plt.plot(Delta2vec, paveraged_r, '-')

plt.figure()
dist = np.linspace(0.0,50,1000)

y = nndistV(1000000, dist, rho)
plt.plot(dist, y)

#plt.figure()
#plt.plot(xvec, Vdip)


plt.show()
