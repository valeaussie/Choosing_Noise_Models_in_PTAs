import numpy as np
import chainconsumer
from chainconsumer import ChainConsumer
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

outdir_1 = '/fred/oz005/users/vdimarco/Choosing_Noise_Models_in_PTAs/sims_sanity_check/sim_3/chains/commonNoise_pl_nocorr_freegam_DE421/1'
outdir_2 = '/fred/oz005/users/vdimarco/Choosing_Noise_Models_in_PTAs/sims_sanity_check/sim_3/chains/commonNoise_pl_nocorr_freegam_DE421/2'
all_data_1 = np.load(f'/fred/oz005/users/vdimarco/Choosing_Noise_Models_in_PTAs/sims_sanity_check/sim_3/chains/commonNoise_pl_nocorr_freegam_DE421/1/chain.npy')
all_data_2 = np.load(f'/fred/oz005/users/vdimarco/Choosing_Noise_Models_in_PTAs/sims_sanity_check/sim_3/chains/commonNoise_pl_nocorr_freegam_DE421/2/chain.npy')


pars_1 = np.loadtxt(outdir_1 + '/pars.txt', dtype=np.unicode_)
pars_2 = np.loadtxt(outdir_2 + '/pars.txt', dtype=np.unicode_)

cpars = ['gw']

all_indices_1 = []
for cp in cpars:
    indices = [cp in p for p in pars_1]
    all_indices_1.append(indices)
indices_1 = np.any(all_indices_1, axis=0)

all_indices_2 = []
for cp in cpars:
    indices = [cp in p for p in pars_2]
    all_indices_2.append(indices)
indices_2 = np.any(all_indices_2, axis=0)

corner_pars = pars_1[indices_1]

data_1 = all_data_1[:, :-4]
data_2 = all_data_2[:, :-4]

data_1 = data_1[:, indices_1]
data_2 = data_2[:, indices_2]


cc = ChainConsumer()
parameters = [r"$\gamma_{\mathrm{CRN}}$", r"$\log_{10}A_{\mathrm{CRN}}$"]
#parameters = list(corner_pars.squeeze())

linestyles = [
    "-",
    "--",
]
linewidths = [
    1.0,
    1.5,
]
shade = [
    True,
    True,
]
shade_gradient = [
    2.0,
    0.5,
]
shade_alpha = [
    0.2,
    0.2,
]
colors = [
    mcolors.CSS4_COLORS['darkorange'],
    mcolors.CSS4_COLORS['mediumblue'],
]
sigmas = [
    1,
    2
]

tfs = 18 #tick
lfs = 27 #label

cc.add_chain(data_1, parameters=parameters, name='Misspecified model').add_chain(data_2, name='Correctly specified model')

cc.configure(summary=False, linestyles=linestyles, linewidths=linewidths, sigmas=sigmas,
            shade=shade, bar_shade=shade, shade_alpha=shade_alpha, shade_gradient=shade_gradient, 
            legend_color_text=True, legend_artists=True, colors=colors, smooth=False,
            max_ticks=7, tick_font_size=tfs, label_font_size=lfs, spacing=1.0, legend_kwargs={'loc': 'lower right'})#, plot_hists=False)

#cfig = cc.plotter.plot(extents=[(3.9, 7), (-10.00, -4.00)], truth=[13./3.,-14.7], figsize=(10, 10))
cfig = cc.plotter.plot(truth=[13./3.,-14.7], figsize=(10, 10))
plt.savefig('sanity_check_sim_3.pdf', bbox_inches='tight')
