import pylab
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import scipy
from qutip import *
from scipy import linalg
from scipy import linalg
plt.style.use('bmh')

plt.rcParams.update({'font.size': 14})

Delta_1 = 1.
Delta_2 = 0.
Omeg_1 = 0.1
Omeg_2 = 5.
gamma_eg = 1.
gamma_re = 2.
gamma_rg = 0.

param_step = 300

Delta2vec = np.linspace(-40.,40., param_step)
Delta1vec = np.linspace(-10.,10., param_step)
t = np.linspace(0., 40., 2000.)

ket_g, ket_e, ket_r = qutrit_basis(); #basis for a three-level system
bra_g = ket_g.dag()
bra_e = ket_e.dag()
bra_r = ket_r.dag()

psi0 = ket_g

#All the relevant sigma matrices can be written as products of kets with bras.
#for example, |e > < g | = ket_e * bra_g.
sig_gg = ket_g * bra_g
sig_ge = ket_g * bra_e
sig_gr = ket_g * bra_r
sig_eg = ket_e * bra_g
sig_ee = ket_e * bra_e
sig_er = ket_e * bra_r
sig_rg = ket_r * bra_g
sig_re = ket_r * bra_e
sig_rr = ket_r * bra_r


rho_ss_list = []

rho_ee_steady = []

rho_eg_steady = []
rho_ge_steady = []

rho_eg_steady_aux = []
rho_ge_steady_aux = []

rho_gg_steady = []

rho_ee_steady_aux = []

for Delta_1 in Delta1vec:
    #Diagonals: Delta_1 |e > < e | + Delta_2 |r > < r|
    H_0 = Delta_1*sig_ee + (Delta_2 + Delta_1)*sig_rr
    #Couplings: Omeg_1 |e > < g| + Omeg_2 | r > < e| + H.c.
    V = Omeg_1 * (sig_eg + sig_ge) + Omeg_2 * (sig_re + sig_er)
    H = H_0 + V

    #two-level Hamiltonian for comparison
    Haux = (Delta_1/2.)*sig_ee + (Omeg_1/2.) * (sig_eg + sig_ge) 
    #Decays
    decay_eg = np.sqrt(gamma_eg)*sig_ge  #e decays into g
    decay_re = np.sqrt(gamma_re)*sig_er
    decay_rg = np.sqrt(gamma_rg)*sig_rg #r decays into e or g
    #The decay into g is included for completeness, but it is often set to 0.

    #result = mesolve(H, psi0, t, [decay_eg, decay_re, decay_rg], [sig_ee])

    rho_ss = steadystate(H, [decay_eg, decay_re, decay_rg], method = 'svd');
    rho_ss_aux = steadystate(Haux, [decay_eg], method = 'svd')
    print((rho_ss).tr())

    #rho_ss_list.append(rho_ss)
    rho_ee_steady.append(rho_ss[1,1])

    rho_eg_steady.append(rho_ss[1,0])
    rho_ge_steady.append(rho_ss[0,1])

    rho_gg_steady.append(rho_ss[0,0])

    rho_eg_steady_aux.append(rho_ss_aux[1,0])
    rho_ge_steady_aux.append(rho_ss_aux[0,1])
    rho_ee_steady_aux.append(rho_ss_aux[1,1])
    print(rho_ss_aux[1,1])
    

rho_eg_steady = np.asarray(rho_eg_steady)
rho_ge_steady = np.asarray(rho_ge_steady)
rho_eg_steady_aux = np.asarray(rho_eg_steady_aux)
rho_ge_steady_aux = np.asarray(rho_ge_steady_aux)
rho_ee_steady_aux = np.asarray(rho_ee_steady_aux)  

print("Here are the matrices")
#print rho_ss_list[0][0,1]

#rhoeg_teor = 2j*Delta1vec*Omeg_1/(gamma_eg*Delta1vec + 2j*(Delta1vec*Delta1vec - Omeg_2*Omeg_2))

rhoeg_teor = 2j*(gamma_re + 2j*(Delta1vec + Delta_2))*Omeg_1
rhoeg_teor = rhoeg_teor/((2j*Delta1vec + gamma_eg)*(2j*(Delta_2 + Delta1vec) + gamma_re) + 4.*Omeg_2*Omeg_2 )

#rhoeg_teorim = 2.*(Delta1vec*Delta1vec*Omeg_1/gamma_eg)/(Delta1vec*Delta1vec + 4*((Delta1vec*Delta1vec - Omeg_2*Omeg_2)/gamma_eg)**2.)

#plt.figure()
#plt.plot(Delta2vec, rho_gg_steady)

plt.figure()
#plt.plot(Delta1vec, rho_ee_steady)
fig = plt.subplot()
fig.set_facecolor('white')
plt.plot(Delta1vec, rho_ge_steady.imag, label = 'exact')
#plt.plot(Delta1vec, rho_ee_steady, label = 'EIT')

#plt.plot(Delta1vec, rho_ge_steady_aux.imag, '--', label = r'$\Im \{ \rho_{ge} \}$')

#plt.plot(Delta1vec, rho_ee_steady_aux, '-.', label = r'2-level atom')
plt.plot(Delta1vec, rhoeg_teor.imag, '-.', label = 'weakly-probing limit' )
#plt.plot(Delta1vec, rhoeg_teorim, '-.')
#plt.xticks([0], [r'$\Delta_1 = 0$'])
#plt.yticks([],[])
plt.legend()
plt.xlabel(r'$\Delta_1$', fontsize = '20')
plt.ylabel(r'$\Im \{ \rho_{ge} \}$', fontsize = '20')
#plt.ylabel(r'$\rho_{ee}$', fontsize = '20')
#plt.title('Fluorescence')
fig.patch.set_facecolor('white')
plt.tight_layout()

#plt.plot(result.times, result.expect[0])
plt.savefig('absorbmoregeneral.pdf')
plt.show()

