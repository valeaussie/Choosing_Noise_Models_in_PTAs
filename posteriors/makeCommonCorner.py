#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:53:00 2019

@author: dreardon

Runs basic white, red, and DM noise model for all pulsars in datadir
"""

import os
import json
import sys
import numpy as np
from enterprise_extensions import model_utils
import matplotlib.pyplot as plt
import corner
#from scipy.signal import correlate

# def acor(arr):
#     arr -= np.mean(arr)
#     auto_correlation = correlate(arr, arr, mode='full')
#     auto_correlation = auto_correlation[auto_correlation.size//2:]
#     auto_correlation /= auto_correlation[0]
#     indices = np.where(auto_correlation<0.5)[0]
#     if len(indices)>0:
#         if indices[0] == 0:
#             indices[0] = 1
#         return indices[0]
#     else:
#         return 1


def make_noise_files(psrname, chain, pars, outdir='noiseFiles/'):
    x = {}
    for ct, par in enumerate(pars):
        x[par] = np.median(chain[:, ct])

    os.system('mkdir -p {}'.format(outdir))
    with open(outdir + '/{}_noise.json'.format(psrname), 'w') as fout:
        json.dump(x, fout, sort_keys=True, indent=4, separators=(',', ': '))


def make_noise_files_maxlike(psrname, chain, pars, outdir='noiseFiles_maxlike/'):
    x = {}
    for ct, par in enumerate(pars):
        ind = np.argmax(chain[:, -4])  # find index of maximum posterior
        x[par] = chain[ind, ct] 

    os.system('mkdir -p {}'.format(outdir))
    with open(outdir + '/{}_noise.json'.format(psrname), 'w') as fout:
        json.dump(x, fout, sort_keys=True, indent=4, separators=(',', ': '))

def make_chain_total(pars, chain_total, outdir = 'noiseFiles/chains/', noise_model=''):
    
    os.system('mkdir -p {}'.format(outdir))
    outname = outdir + '{}_chain.npy'.format(noise_model.replace('bugfix_', ''))
    parname = outdir + '{}_pars.npy'.format(noise_model.replace('bugfix_', ''))
    np.save(parname, pars)
    print('saving chain')
    np.save(outname, chain_total)
    return True

def plot_corner(chain, pars, cpars, noise_model='', outdir='noiseFiles/', ylim=(None, None)):

    chain = chain[:, :-4]
    if cpars is not None:

        all_indices = []

        for cp in cpars:
            indices = [cp in p for p in pars]
            all_indices.append(indices)
        indices = np.any(all_indices, axis=0)

        corner_pars = pars[indices]
        corner_file_label = noise_model+'_' + '_'.join(cpars)
    else:
        # plot all pars
        indices = np.array([True for p in pars])
        corner_file_label = noise_model

    chain_corner = chain[:, indices]
    fig = corner.corner(chain_corner, bins=30, labels=corner_pars,
                        quantiles=(0.16, 0.5, 0.84), show_titles=True)
    for ax in fig.axes:
        xlab = ax.get_xlabel()
        ylab = ax.get_ylabel()
        ti = ax.get_title()
        ax.set_title(ti, fontsize=9)
        ax.set_xlabel(xlab, fontsize=9)
        ax.set_ylabel(ylab, fontsize=9)
        #if "log10_A" in ylab:
        #    ax.set_ylim(ylim)
        #elif "log10_A" in xlab:
        #    ax.set_xlim(ylim)
    os.system('mkdir -p {}'.format(outdir))
    figsavename = outdir + '/' + corner_file_label + '.png'
    print("saving figure")
    print(figsavename)
    plt.savefig(figsavename, dpi=300, bbox_inches='tight')
    print("done saving figure")


#Run with 3 arguments
dir = sys.argv[1]
noise_model = sys.argv[2]#'commonNoise_pl_nocorr_fixgam_freespec_monopole'#commonNoise_pl_nocorr_fixgam_pl_hdnoauto_fixgam'
#cpars = None
cpars = ['gw']
#If more than 3 parameters are used parse the third argument in a list of parameters
#either comma separated of space-separated
if len(sys.argv) > 3:
    if ',' in sys.argv[3]:
        cpars = [p.strip() for p in sys.argv[3].split(',')]
    else:
        cpars = [p.strip() for p in sys.argv[3].split(' ')]

#datadir = os.path.abspath("./data/" + str(dir))
datadir = os.path.abspath(str(dir) + '/chains')
#datadir = os.path.abspath(str(dir))
print('printing datadir')
print(datadir)

max_nchain = 1  # number of chains to read

#if _chain.npy already exists load the parameters
if os.path.exists(datadir + '/' +  noise_model + '/_chain.npy'):
    chain_total = np.load(datadir + '/' + noise_model +'/_chain.npy')
    pars = np.load(datadir + '/' + noise_model + '/_pars.npy')
    max_nchain = 1

first_chain = True
for i in range(1, max_nchain+1):
    outdir = os.path.abspath(str(dir))
    print('printing outdir')
    print(outdir)
    chainfile = outdir + '/' + noise_model +  '/chain_1.txt'
    print('printing chainfile')
    print(chainfile)

    if not os.path.exists(chainfile):
        continue
    if os.path.getsize(chainfile) == 0:
        continue
    
    # print('loading {}...'.format(chainfile))
    # try:
    #     try:
    #         chain_i = np.loadtxt(chainfile).squeeze()
    #     except FileNotFoundError:
    #         print("chain file missing. continuing")
    #         continue
    #     print('loaded {}'.format(chainfile))
    # except ValueError as e:
    #     print(e)
    #     continue
    # try:
    #     pars = np.loadtxt('/pars.txt', dtype=np.unicode_)
    # except FileNotFoundError:
    #     print("pars file not found - attempting to use old pars file")

    chain_i = np.loadtxt(chainfile).squeeze()
    pars = np.loadtxt(outdir + '/' + noise_model + '/pars.txt', dtype=np.unicode_)
    pp = model_utils.PostProcessing(chain_i, pars)
    # pp.plot_trace()
    # plt.savefig(chainfile.replace('.txt', '.png'))
    # plt.close()

    # Burn larger of 25% or first 25000
    burn = int(max([0.25*chain_i.shape[0], 25000]))
    #burn = 150
    chain = chain_i[burn:, :]
    print('here')
    print("length of chain: ", len(chain))
    if chain.size == 0:
        continue
    if first_chain:
        chain_total = chain
        first_chain = False
    else:
        chain_total = np.concatenate((chain_total, chain)).squeeze()
    #print('chain')
    #print(chain)
    #pp = model_utils.PostProcessing(chain_i, pars)
    #pp.plot_trace()
    #plt.savefig(chainfile.replace('.txt', '.png'))
    #plt.close()
    #print('saved {}'.format(chainfile.replace('.txt', '.png')))

    # Burn larger of 25% or first 25000
    #    burn = int(max([0.25*chain_i.shape[0], 25000]))
    #nburnmin = 1500
    #burn = int(max([0.25*chain_i.shape[0], nburnmin]))
    #print(burn)
    #gw_log10A_ind = np.squeeze(np.argwhere(np.logical_and(['gw' in p for p in pars], ['log10_A' in p for p in pars])))

    #thin = 0
    #3mx = 0

    #good_pars = 0
    #inds = np.argwhere([('nmodel' not in p and 'gw' in p and 'bin' not in p) for p in pars])
    # for indi in inds:
    #     try:
    #         thin_i = acor(chain_i[burn:, indi])
    #     except ValueError:
    #         continue
    #     thin += thin_i
    #     good_pars += 1
    #     #print("Autocorrelation on {} is {}".format(pars[ind], thin_i))
    # if good_pars == 0: continue
    # thin /= good_pars
    # thin = int(thin)    
    # #thin = 50
    # print('autocorrelation length = {}'.format(thin))
    # chain = chain_i[burn::thin, :]


    #gw_log10A_ind = np.squeeze(np.argwhere(np.logical_and(['gw' in p for p in pars], ['log10_A' in p for p in pars])))
    
    #good_inds = chain_i[:, gw_log10A_ind] < -13.6
    #chain_i = chain_i[good_inds]
    #chain = chain_i[burn:, :]
    #if chain.size == 0:
    #    print('fewer than {} samples. Continuing'.format(nburnmin))
    #    continue
    #print(chain.shape)
    
    #if cpars is not None:
    #    plot_corner(chain, pars, cpars, outdir=outdir)
        #print('Made corner plot for {} using chain {}'.format(psrname, str(i)))

    #if first_chain:
    #    chain_total = chain
    #    first_chain = False
    #else:
    #    if chain.shape[1] != chain_total.shape[1]:
    #        continue
    #    chain_total = np.concatenate((chain_total, chain)).squeeze()

# Now, save noise files
np.save(outdir + '/' + noise_model + "/chain.npy", chain_total)
#cpars = None
try:
    #make_noise_files(pars1, chain_total, pars, outdir=datadir+'/noiseFiles/')
    #make_noise_files_maxlike(psrname, chain_total, pars, outdir=datadir+'/noiseFiles_maxlike/')
    print(chain_total.shape)
    #print('Wrote noise files for {}'.format(psrname))
    cpars = ['gw']
    #make_chain_total(pars, chain_total, outdir = datadir + '/corners/', noise_model=noise_model)
    if cpars is not None:
        plot_corner(chain_total, pars, cpars, noise_model=noise_model, outdir=outdir)#, ylim=(None, -13.6))
        print('Made corner plot')
except Exception as e:
    print(e)
    #print('No valid chains for {}'.format(psrname))
