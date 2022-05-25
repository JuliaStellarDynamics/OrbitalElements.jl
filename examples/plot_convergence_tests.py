"""
plasma_plot.py
    visualise some of the plasma calculations

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


A = np.genfromtxt("NINTconvergence.txt",delimiter=',')


plt.plot(A[1:,0],np.log10(np.abs(A[1:,1]-A[0,1])),color='grey')
plt.plot(A[1:,0],np.log10(np.abs(A[1:,2]-A[0,2])),color='grey',linestyle='dashed')

plt.xlabel("NINT")
plt.ylabel("Omega1")
#plt.title("Isochrone mode test")
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
