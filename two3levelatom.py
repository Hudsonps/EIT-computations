import pylab
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import scipy
from qutip import *
#from scipy import linalg
#from two3levelHamiltonian import *
#from nndist import *
plt.style.use('bmh')
from itertools import cycle







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

















#plt.rcParams.update({'font.size': 14})

#Code about two three-level atoms

def ind(i,j, dim):
    #map of i and j into a single index, runs through columns first, and rows second
    return i*dim + j

t = np.linspace(0., 40., 2000.)  #vector for time evolution

#Delta_1 = 1.
Delta_2 = 0.
Omeg_1 = 1.
#coupling from g to e
Omeg_2 = 2.    #coupling from e to r
gamma_eg = 1.0  #decay from e to g
gamma_re = 0.0000000  #decay from r to e
gamma_rg = 0.0  #decay r to g
gamma_egcol = 0.0  #collective decay, e to g
#Vdip = 0.   #dipole interaction betweel levels r

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


#operators whose expectation values are of interest
popinr = sig1[ind(2,2,dim)] + sig2[ind(2,2,dim)]
popine = sig1[ind(1,1,dim)] + sig2[ind(1,1,dim)]
poping = sig1[ind(0,0,dim)] + sig2[ind(0,0,dim)]

coe_eg = sig1[ind(1,0,dim)] + sig2[ind(1,0,dim)]
coe_er = sig1[ind(2,1,dim)] + sig2[ind(2,1,dim)]
    
#obs = sig1[ind(1,0,dim)] + sig2[ind(1,0,dim)]
    
popinrr = sig1[ind(2,2,dim)]*sig2[ind(2,2,dim)]

num_step = 200
Delta1vec = np.linspace(-8,8., num_step)

Vdip = np.arange(0, 8, 1.)


plt.figure()
fig = plt.subplot()
fig.set_facecolor('white')

plt.xlabel(r'$\Delta_1$')
plt.ylabel("Population in "+ r'$|e \rangle$')
#plt.ylabel(r"$\langle \sigma_{rr}^{(1)} \otimes \sigma_{rr}^{(2)} \rangle$") 
#plt.title("r'$\langle \sigma_{eg}^{(1)} +  \sigma_{eg}^{(2) \rangle}$)


lines = ["-", "--", "-.", ":","-"]
linecycler = cycle(lines)

for V in Vdip:
    trace = []
    pop_g = []
    pop_e = []
    pop_r = []
    pop_rr = []
    coherence = []
    coherence_er = []
    for Delta_1 in Delta1vec:
        H = Hamiltonian(Delta_1, Delta_2, Omeg_1, Omeg_2, 1.*V)
        #Heff = Hamiltonianeff(Delta_1, Delta_2, Omeg_1, Omeg_2, Vdip)
        decay_ops = decay(gamma_eg, gamma_re, gamma_rg, gamma_egcol)
        
        psi0 = tensor(kets_1atom[0], kets_1atom[0]) #initial state

        rho_ss = steadystate(H, decay_ops, method="eigen") #steady state, given H and decay_ops
        trace.append((rho_ss).tr())

        pop_rr.append((rho_ss*popinrr).tr())
        pop_r.append((rho_ss*popinr).tr())
        pop_g.append((rho_ss*poping).tr())
        pop_e.append((rho_ss*popine).tr())

        coherence.append(-1j*(rho_ss*coe_eg).tr())
        coherence_er.append(-1j*(rho_ss*coe_er).tr())

    plt.plot(Delta1vec, pop_e, next(linecycler), label= '$V_0=$'+str(V))

    #plt.plot(Delta1vec, pop_rr, next(linecycler), label= '$V_0=$'+str(V))
    
    #plt.plot(Delta1vec, coherence, '-.')
    #plt.plot(Delta1vec, coherence_er, '--')
    #plt.plot(Delta1vec, pop_r, label= 'r$V=$'+str(V))
    #plt.plot(Delta1vec, pop_g, label= 'r$V=$'+str(V))

#plt.text(x=4.1, y=0.08, s = '$\Omega_2 = 2$', fontsize = 14)
#plt.text(x=4.1, y=0.06, s = '$\Omega_1 = \gamma = 1$', fontsize = 14)
#plt.text(x=4.1, y=0.04, s = '$\Delta_2 = 0$', fontsize=14)
plt.text(x=-8, y = 0.48, s= '$\Omega_2 = 2$', fontsize=14)    
plt.text(x=-8, y = 0.44, s= '$\Omega_1 = \gamma = 1$', fontsize=14)
plt.text(x=-8, y = 0.40, s= '$\Delta_2 = 0$', fontsize=14)
plt.title("Fluorescence profile")

#plt.title("Double occupation of the Rydberg state")

plt.legend()

plt.savefig("doubleoccupation.pdf")


    
exp1 = popine
H = Hamiltonian(0., 0., Omeg_1, Omeg_2, 5.)
decay_ops = decay(gamma_eg, gamma_re, gamma_rg, gamma_egcol)
result = mesolve(H, psi0, t, decay_ops, [exp1])

plt.figure()
plt.plot(t, result.expect[0], color = 'orange')


############################



#plt.figure()
#plt.plot(Delta1vec, trace, '-', color = 'blue')

#plt.figure()
#plt.plot(Delta1vec, pop_r, '-', label = r'$|r \rangle \langle r|_1 + | r \rangle \langle r |_2$')
#plt.plot(Delta1vec, pop_rr, '--', label = r'$|r \rangle \langle r|_1 \otimes | r \rangle \langle r |_2$')
#plt.plot(Delta2vec, pop_g)
#plt.plot(Delta1vec, pop_e, label= r'$|e \rangle \langle e|_1 + | e \rangle \langle e |_2$')


#gamma_eg = 2.
#rhoeg_teor = 2j*(gamma_re + 2j*(Delta1vec + Delta_2))*Omeg_1
#rhoeg_teor = rhoeg_teor/((2j*Delta1vec + gamma_eg)*(2j*(Delta_2 + Delta1vec) + gamma_re) + 4.*Omeg_2*Omeg_2 )

#plt.plot(Delta1vec, 2.*rhoeg_teor.imag, '-.', label = r'weakly-probing limit, $V=0$', color='black' )
plt.legend()

#plt.figure()
#dist = np.linspace(0,20,1000)
#C = 1000000
#n = 1/1000.
#y = nndistV(C, dist, 1/1000.)
#plt.plot(dist, y)


plt.show()
