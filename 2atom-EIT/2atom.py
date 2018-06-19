import pylab
import random
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import linalg
plt.style.use('ggplot')



from matplotlib import rc
rc('text', usetex=True)

vec_idx = 2

mpl.use('pdf')

c2 = 1.
c1 = 0.1
Delta1 = 5.
V = 5.
ceff = -c1*c2/Delta1

dim = 6

lambda0 = []
lambda1 = []
lambda2 = []
lambda3 = []
lambda4 = []
lambda5 = []

lambda0eff = []
lambda1eff = []
lambda2eff = []

gg = []
ge = []
gr = []
ee = []
er = []
rr = []

ggeff = []
greff = []
rreff = []

#Delta2min = c1*c1/Delta1
Delta2particular = -V/2 + (c2*c2-c1*c1)/Delta1
#Delta2particular = 0.
Delta2min = Delta2particular - 1.
Delta2max = Delta2particular + 1.
Delta2step = 0.01

for Delta2 in np.arange(Delta2min, Delta2max, Delta2step):
    H = np.zeros((dim,dim))
    H[0,1] = np.sqrt(2)*c1
    H[1,1] = Delta1
    H[1,2] = c2
    H[1,3] = np.sqrt(2)*c1
    H[2,2] = Delta2
    H[2,4] = c1
    H[3,3] = 2.*Delta1
    H[3,4] = np.sqrt(2)*c2
    H[4,4] = Delta1 + Delta2
    H[4,5] = np.sqrt(2)*c2
    H[5,5] = 2.*Delta2 + V

    for i in range(0, dim):
        for j in range(0, i):
            H[i,j] = H[j,i]

    
    Heff = np.zeros((3,3))
    Heff[0,0] = -2.*c1*c1/Delta1
    Heff[0,1] = np.sqrt(2.)*ceff
    Heff[1,0] = Heff[0,1]
    Heff[1,1] = Delta2 - c2*c2/Delta1 - c1*c1/Delta1
    Heff[1,2] = np.sqrt(2.)*ceff
    Heff[2,1] = Heff[1,2]
    Heff[2,2] = 2.*Delta2 - 2.*c2*c2/Delta1 + V

    
    
    eigenval, eigenvec = scipy.linalg.eig(H)
    eigenvaleff, eigenveceff = scipy.linalg.eig(Heff)

    sort_perm = eigenval.argsort()
    eigenval.sort()     # <-- This sorts the list in place.
    eigenvec = eigenvec[:, sort_perm]


    sort_perm = eigenvaleff.argsort()
    eigenvaleff.sort()
    eigenveceff = eigenveceff[:, sort_perm]
    

    lambda0.append(eigenval[0])
    lambda3.append(eigenval[3])
    lambda4.append(eigenval[4])
    lambda5.append(eigenval[5])  
    lambda1.append(eigenval[1])  
    lambda2.append(eigenval[2])

    gg.append(eigenvec[0,vec_idx])
    ge.append(eigenvec[1,vec_idx])
    gr.append(eigenvec[2,vec_idx])
    ee.append(eigenvec[3,vec_idx])
    er.append(eigenvec[4,vec_idx])
    rr.append(eigenvec[5,vec_idx])
    

    lambda0eff.append(eigenvaleff[0])
    lambda1eff.append(eigenvaleff[1])
    lambda2eff.append(eigenvaleff[2])

    ggeff.append(eigenveceff[0,vec_idx])
    greff.append(eigenveceff[1,vec_idx])
    rreff.append(eigenveceff[2,vec_idx])

    
Delta2 = np.arange(Delta2min, Delta2max, Delta2step)
ee = np.asarray(ee)
er = np.asarray(er)
ge = np.asarray(ge)
gg = np.asarray(gg)
gr = np.asarray(gr)
rr = np.asarray(rr)

ggeff = np.asarray(ggeff)
greff = np.asarray(greff)
rreff = np.asarray(rreff)


plt.plot(Delta2, lambda0, linewidth = 2.0)
plt.plot(Delta2, lambda1, linewidth = 2.0)
plt.plot(Delta2, lambda2, linewidth = 2.0)
plt.xlabel(r'$\Delta_2$')
plt.ylabel(r'$\epsilon$')
#plt.plot(Delta2, lambda3, linewidth = 2.0)
#plt.plot(Delta2, lambda4, linewidth = 2.0)
#plt.plot(Delta2, lambda5, linewidth = 2.0)

#plt.plot(Delta2, lambda0eff, '--', color = 'black')
#plt.plot(Delta2, lambda1eff, '--', color = 'black')
#plt.plot(Delta2, lambda2eff, '--', color = 'black')

plt.figure()
plt.plot(Delta2, ee*ee, '-', linewidth = 2.0, label = r'$| 2_e \rangle$')
plt.plot(Delta2, er*er, '--', linewidth = 2.0, label = r'$|1_e, 1_r \rangle$')
plt.plot(Delta2, ge*ge, '-.', linewidth = 2.0, label = r'$|1_g, 1_e \rangle$')
plt.xlabel(r'$\Delta_2$')
plt.ylabel('Population')
plt.legend()

plt.figure()
plt.plot(Delta2, gg*gg, '-', linewidth = 2.0, label = r'$| 2_g \rangle$')
plt.plot(Delta2, gr*gr, '--', linewidth = 2.0, label = r'$|1_g, 1_r \rangle$')
plt.plot(Delta2, rr*rr, '-.', linewidth = 2.0, label = r'$|2_r \rangle$')
plt.xlabel(r'$\Delta_2$')
plt.ylabel('Population')
#plt.plot(Delta2, ggeff*ggeff, '--', linewidth = 1.0, color = 'black')
#plt.plot(Delta2, greff*greff, '--', linewidth = 1.0, color = 'black')
#plt.plot(Delta2, rreff*rreff, '--', linewidth = 1.0, color = 'black')
plt.legend()

print eigenvec[:,1]
print eigenvec[:,2]

#print eigenval
print eigenvaleff


plt.show()










