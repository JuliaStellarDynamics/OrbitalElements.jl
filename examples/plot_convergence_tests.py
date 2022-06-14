"""
plot_convergence_tests.py
    visualisation of some of the convergence tests for orbital elements

"""

import numpy as np

# enable hdf5 reading
import h5py

# plotting basics
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# override some parameters for a 'house style'
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.width'] = 0.75
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.width'] = 0.75
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True

A = np.genfromtxt("ELconvergence.txt",delimiter=',')

plt.figure(figsize=(4.5,3))

# make predictions in grey, real values in black
#plt.plot(np.log10(A[1:,0]),np.log10(1/(A[1:,0]**2))-0.7,color='grey')
plt.plot(np.log10(A[1:,0]),np.log10(np.abs(A[1:,1])),color='black')
plt.plot(np.log10(A[1:,0]),np.log10(np.abs(A[1:,2])),color='black')

#plt.plot((A[:,0]),np.log10(np.abs(A[:,1])),color='black')
#plt.plot((A[:,0]),np.log10(np.abs(A[:,2])),color='black')

plt.xlabel("log $e$")
plt.ylabel("log $\epsilon_{E}$ (solid)")
#plt.title("Isochrone mode test")
plt.tight_layout()
plt.savefig('figures/el_convergence.png')



A = np.genfromtxt("NINTconvergence.txt",delimiter=',')

plt.figure(figsize=(4.5,3))

# make predictions in grey, real values in black
plt.text(0.,-2.7,'1/$N_\Omega^2$',color='grey',size=10)
plt.plot(np.log10(A[1:,0]),np.log10(1/(A[1:,0]**2))-0.7,color='grey')
plt.plot(np.log10(A[1:,0]),np.log10(np.abs(A[1:,1]-A[0,1])),color='black')

plt.text(0.,-11.0,'1/$N_\Omega^0$',color='grey',size=10)
plt.plot(np.log10(A[1:,0]),np.log10(1/(A[1:,0]**0))-10.,color='grey',linestyle='dashed')
plt.plot(np.log10(A[1:,0]),np.log10(np.abs(A[1:,2]-A[0,2])),color='black',linestyle='dashed')

plt.text(0.,-5.9,'1/$N_\Omega^4$',color='grey',size=10)
plt.plot(np.log10(A[1:,0]),np.log10(1/(A[1:,0]**4))-5.9,color='grey',linestyle='dotted')
plt.plot(np.log10(A[1:,0]),np.log10(np.abs(A[1:,3]-A[0,3])),color='black',linestyle='dotted')

plt.xlabel("log N$_\Omega$")
plt.ylabel("log $\epsilon_{\Omega_1}$ (solid)\nlog $\epsilon_{\Omega_2}$ (dashed)\nlog $\epsilon_{J_r}$ (dotted)")
#plt.title("Isochrone mode test")
plt.tight_layout()
plt.savefig('figures/nint_convergence.png')


A = np.genfromtxt("NINTarray.txt",delimiter=',')

na,ne = 250,100
avals = np.log10(A[:,0]).reshape(na,ne)
evals = A[:,1].reshape(na,ne)
nintv = A[:,2].reshape(na,ne)

plt.figure(figsize=(5,4))
plt.pcolormesh(avals,evals,np.log10(nintv))
plt.xlabel('log a')
plt.ylabel('e')
plt.colorbar(label='log N')
plt.title('iterations to 1e-5 accuracy\n(max 256)')
plt.tight_layout()
plt.savefig('figures/nint_rate.png')
