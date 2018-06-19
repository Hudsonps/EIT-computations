import pylab
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import scipy
from qutip import *
from scipy import linalg
plt.style.use('ggplot')

t = np.linspace(0., 4., 2000.)

psi0 = basis(2,0) #that means ground state. (2,0) means excited.

#Hamiltonian in the rotating wave frame
# -Delta/2 sigma_z  + Omega/2 sigma_x

Delta = 0.
Omega = 1.

H = -(Delta/2.)*sigmaz() + (Omega/2.)*sigmax()

#Dissipation
gamma = 2.
diss = np.sqrt(gamma)*sigmam()

result = mesolve(H, psi0, t, [diss], [sigmaz()])

env = np.exp(-gamma*t)

plt.plot(result.times, result.expect[0])
plt.plot(t, 2.*env - 1, '--')
plt.show()
