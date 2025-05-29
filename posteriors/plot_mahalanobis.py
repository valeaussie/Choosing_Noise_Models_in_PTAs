import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

array_of_additions = np.arange(3, 31, 3)

data = np.load('/fred/oz002/vdimarco/P3_spectral_index/posteriors/varying_chains_mahalanobis.npy', allow_pickle=True)


print('plotting')
plt.figure()
plt.xlabel('Number of pulsars with chromatic noise', fontsize=13)
plt.ylabel('Distance', fontsize=13)
plt.plot(array_of_additions, data[0], color="blue")
for j in range(14):
    #if j!=13 and j!=8:
        plt.plot(array_of_additions, data[j], color="grey",
                linewidth=0.5, label='Shuffle ' + str(j))
#plt.legend(fontsize=12)
plt.savefig('/fred/oz002/vdimarco/P3_spectral_index/posteriors/varying_chains_mahalanobis.pdf')