import pylab
import random
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import linalg


C2 = 1.
c1 = 1.
D = 1

lambda0 = []
lambda1 = []
lambda2 = []
lambda0eff = []
lambda1eff = []
for w in np.arange(-1,1,0.01):
    
    H = np.matrix([[w, c1, 0],[c1,D,C2],[0,C2,0]])
    
    Heff = np.matrix([[w-c1*c1/D, -c1*C2/D],[-c1*C2/D, -C2*C2/D]])
    
    eigenval, eigenveg = scipy.linalg.eig(H)
    eigenvaleff, eigenvegeff = scipy.linalg.eig(Heff)
    
    eigenval = np.sort(eigenval)
    eigenvaleff = np.sort(eigenvaleff)
    
    lambda0.append(eigenval[0])  
    lambda1.append(eigenval[1])  
    lambda2.append(eigenval[2])
    
    lambda0eff.append(eigenvaleff[0])
    lambda1eff.append(eigenvaleff[1])
    

w = np.arange(-1,1,0.01)  

plt.plot(w, lambda0, linewidth = 2.0, color = 'red', )
plt.plot(w, lambda1, linewidth = 2.0, color = 'magenta')
plt.plot(w, lambda0eff, '--', color = 'black')
plt.plot(w, lambda1eff, '--', color = 'black')
plt.plot(w,lambda2)
plt.show()










