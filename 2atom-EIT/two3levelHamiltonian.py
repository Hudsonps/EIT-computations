import pylab
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import scipy
from qutip import *
from scipy import linalg
plt.style.use('ggplot')

def ind(i,j, dim):
    #map of i and j into a single index
    #this index runs through columns first, and rows second
    return i*dim + j

dim = 3 #each atom is a three-level system
idx = range(dim)

kets_1atom = qutrit_basis(); #basis for a three-level system
bras_1atom = []
for i in idx:
    bras_1atom.append((kets_1atom[i]).dag()) #generates the dual basis for one atom
     

#Sigma matrices for one atom only.
sig_1atom = []
for i in idx:
    for j in idx:
        sig_1atom.append(kets_1atom[i]*bras_1atom[j])

#Sigma matrices operating on atom 1, while atom 2 does nothing, and vice-versa.
sig1 = []
sig2 = []
for sig in sig_1atom:
    sig1.append(tensor(sig, identity(dim)))
    sig2.append(tensor(identity(dim), sig))
#The indices go from 1 to dim*dim. sig1[0] = sig1_gg; sig1[1] = sig1_ge, etc.

def Hamiltonian(Delta_1, Delta_2, Omeg_1, Omeg_2, Vdip):
    #atomic Hamiltonian
    H_1 = Delta_1*sig1[ind(1,1,dim)] + (Delta_2 + Delta_1)*sig1[ind(2,2,dim)] #atom 1  (* identity in atom 2)
    H_2 = Delta_1*sig2[ind(1,1,dim)] + (Delta_2 + Delta_1)*sig2[ind(2,2,dim)] #atom 2  (* identity in atom 1)
    #coupling with light
    V_1 =  Omeg_1 * (sig1[ind(1,0,dim)] + sig1[ind(0,1,dim)]) + Omeg_2 * (sig1[ind(2,1,dim)] + sig1[ind(1,2,dim)])
    V_2 = Omeg_1 * (sig2[ind(1,0,dim)] + sig2[ind(0,1,dim)]) + Omeg_2 * (sig2[ind(2,1,dim)] + sig2[ind(1,2,dim)])
    #interaction between the Rydberg levels
    U = Vdip*sig1[ind(2,2,dim)]*sig2[ind(2,2,dim)]
    #Full Hamiltonian
    H = H_1 + H_2 + V_1 + V_2 + U
    return H

def Hamiltonianeff(Delta_1, Delta_2, Omeg_1, Omeg_2, Vdip):
    dim_eff = 2
    #Effective 4x4 Hamiltonian
    Omeg_eff = -Omeg_1*Omeg_2/Delta_1
    Heff_1 =  -((Omeg_1**2 - Omeg_2**2)/(2.*Delta_1))*tensor(sigmaz(), identity(dim_eff)) + Omeg_eff*tensor(sigmax(),identity(2))
    Heff_1 = Heff_1 - (Omeg_1**2 + Omeg_2**2)/(2.*Delta_1)
    Heff_2 = -((Omeg_1**2 - Omeg_2**2)/(2.*Delta_1))*tensor(identity(dim_eff), sigmaz()) + Omeg_eff*tensor(identity(2), sigmax())
    Heff_2 = Heff_2 - (Omeg_1**2 + Omeg_2**2)/(2.*Delta_1)
    Heff = Heff_1 + Heff_2
    return Heff

def decay(gamma_eg, gamma_re, gamma_rg, gamma_egcol):
    decay_eg1 = np.sqrt(gamma_eg)*sig1[ind(0,1,dim)]
    decay_eg2 = np.sqrt(gamma_eg)*sig2[ind(0,1,dim)]
    decay_egcol = np.sqrt(gamma_egcol)*(sig1[ind(0,1,dim)]+sig2[ind(0,1,dim)])
    decay_re1 = np.sqrt(gamma_re)*sig1[ind(1,2,dim)]
    decay_re2 = np.sqrt(gamma_re)*sig2[ind(1,2,dim)]
    return [decay_eg1, decay_eg2, decay_egcol, decay_re1, decay_re2]

