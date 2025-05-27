import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import os
import sys
import json
#from scipy.spatial import distance
from scipy.spatial.distance import mahalanobis

########################
# Defintion of functions
########################

# function to find data for specific parameters
def find_data(pars, cpars, data):
    all_indices = []
    for cp in cpars:
        indices = [cp in p for p in pars]
        all_indices.append(indices)
    indices = np.any(all_indices, axis=0)

    data = data[:, :-4]
    data = data[:, indices]

    return data

def find_mahalanobis(data1, data2):
    # find means
    mean1_gamma = np.mean(data1[:,0], axis=0)
    mean1_logA = np.mean(data1[:,1], axis=0)
    mean2_gamma = np.mean(data2[:,0], axis=0)
    mean2_logA = np.mean(data2[:,1], axis=0)
    mean1 = [mean1_gamma, mean1_logA]
    mean2 = [mean2_gamma, mean2_logA]
    #find covariance matrices
    covar1 = np.cov(data1.T)
    covar2 = np.cov(data2.T)
    covar = (np.linalg.inv(covar1 + covar2))/2
    return mahalanobis(mean1, mean2, covar)

# parameters
cpars = ['gw']
array_of_additions = np.arange(3, 31, 3)

mean_distances_1 = {}
mean_distances_2 = {}
for i in range(0, 100):
    rootdir = f"/fred/oz005/users/vdimarco/Choosing_Noise_Models_in_PTAs/sims_sanity_check/sim_{i}/chains/commonNoise_pl_nocorr_freegam_DE421"

    #print('loading chains')
    ### import data for simple model with no CH
    dir_data_no_CH = rootdir + '/3'
    dir_data_with_CH_1 = rootdir + '/4'
    dir_data_with_CH_2 = rootdir + '/5'
    outdir = "/fred/oz005/users/vdimarco/Choosing_Noise_Models_in_PTAs/posteriors/mahalanobis_distances"
    os.makedirs(outdir, exist_ok=True)

    # print(dir_data_no_CH)
    # print(dir_data_with_CH)
    # print(outdir)

    #dir_data_no_CH = '/fred/oz002/vdimarco/P3_spectral_index/sims/sims_for_paper/WN_TN_DM/chains/commonNoise_pl_nocorr_freegam_DE421/1'
    all_data_no_CH = np.load(dir_data_no_CH + '/chain.npy')
    all_data_with_CH_1 = np.load(dir_data_with_CH_1 + '/chain.npy')
    all_data_with_CH_2 = np.load(dir_data_with_CH_2 + '/chain.npy')
    pars_no_CH = np.loadtxt(dir_data_no_CH + '/pars.txt', dtype=np.unicode_)
    pars_with_CH_1 = np.loadtxt(dir_data_with_CH_1 + '/pars.txt', dtype=np.unicode_)
    pars_with_CH_2 = np.loadtxt(dir_data_with_CH_1 + '/pars.txt', dtype=np.unicode_)
    data_no_CH = find_data(pars_no_CH, cpars, all_data_no_CH)
    data_with_CH_1 = find_data(pars_with_CH_1, cpars, all_data_with_CH_1)
    data_with_CH_2 = find_data(pars_with_CH_2, cpars, all_data_with_CH_2)


    distance_1 = find_mahalanobis(data_no_CH, data_with_CH_1)
    distance_2 = find_mahalanobis(data_with_CH_1, data_with_CH_2)
    mean_distances_1[i] = distance_1
    mean_distances_2[i] = distance_2
    print(f"mean distance {i}: ", distance_1, distance_2)


json_outfile_1 = f"{outdir}/mean_mahalanobis_distances_sanity_check_1.json"
json_outfile_2 = f"{outdir}/mean_mahalanobis_distances_sanity_check_2.json"
with open(json_outfile_1, 'w') as f:
    json.dump(mean_distances_1, f, indent=4)
with open(json_outfile_2, 'w') as f:
    json.dump(mean_distances_2, f, indent=4)


