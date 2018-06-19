import pylab
import random
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import linalg
plt.style.use('bmh')

omeg = np.arange(-0.6,0.6, 0.01)



C2 = 1.
c1 = 1.
D = 0.5

k = omeg*(1+ (c1*c1/(C2*C2))*(1/(1 + ((D*omeg-omeg*omeg)/(C2*C2)))))


lambda0 = []
lambda1 = []
lambda2 = []
for w in np.arange(-10,10,0.01):
    
    H = np.matrix([[w, c1, 0],[c1,D,C2],[0,C2,0]])
    
    eigenval, eigenveg = scipy.linalg.eig(H)
    eigenval = np.sort(eigenval)
    
    lambda0.append(eigenval[0])  
    lambda1.append(eigenval[1])  
    lambda2.append(eigenval[2])
    

w = np.arange(-10,10,0.01)

plt.xticks([-1, -0.5, 0, 0.5, 1], [-1, -0.5, 0, 0.5, 1])
#plt.yticks([],[])

plt.plot(w, lambda2, linewidth = 2.0, ls = '--', color = 'blue', label = r'$\epsilon_+ (k)$')
plt.plot(w, lambda1, linewidth = 2.0, ls = '-', color = 'magenta',  label = r'$\epsilon_D (k)$' )
plt.plot(w, lambda0, linewidth = 2.0, ls = '-.', color = 'red', label = r'$\epsilon_- (k)$' )
plt.legend(fontsize = '12')
plt.xlabel(r'$k$', fontsize = '20')
plt.ylabel(r'$\epsilon$', fontsize = '20')

plt.plot(k, omeg, '--', color = 'orange')
#plt.savefig('darkpoldisp.pdf')
plt.show()











