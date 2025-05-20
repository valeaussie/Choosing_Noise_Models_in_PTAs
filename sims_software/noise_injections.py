#############################################
## Created on Fri August 30 13:29:00 2024
##
## Author: Valentina Di Marco
##
## Adds noise sources to simulated TOAs
##############################################

import libstempo as T
import libstempo.toasim as LT
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt
import math
import os
import glob
import json
import random

### Units
day = 24 * 3600 # Days in seconds
year = 365.25 * day # Year in seconds

AU = sc.astronomical_unit
c = sc.speed_of_light
pc = sc.parsec
AU_light_sec = AU/c
AU_pc = AU/pc

datadir = '/fred/oz005/users/vdimarco/Choosing_Noise_Models_in_PTAs/sims_software/WN_PTASim/output/real_0'

## Pulsar names
psrnames = ['J0030+0451',
            'J0125-2327',
            'J0437-4715',
            'J0613-0200',
            'J0614-3329',
            'J0711-6830',
            'J0900-3144',
            'J1017-7156',
            'J1022+1001',
            'J1024-0719',
            'J1045-4509',
            'J1125-6014',
            'J1446-4701',
            'J1545-4550',
            'J1600-3053',
            'J1603-7202',
            'J1643-1224',
            'J1713+0747',
            'J1730-2304',
            'J1744-1134',
            'J1832-0836',
            'J1857+0943',
            'J1902-5105',
            'J1909-3744',
            'J1933-6211',
            'J1939+2134',
            'J2124-3358',
            'J2129-5721',
            'J2145-0750',
            'J2241-5236']

########################
### Functions definition
########################

### Disperiosn mesure noise
def add_dm(psr, A, gamma, idx=-2, components=30, seed=None):
    """Add dispersion measure noise variations with P(f) = A^2 / (12 pi^2) (f year)^-gamma,
    using `components` Fourier bases.
    Optionally take a pseudorandom-number-generator seed."""    

    if seed is not None:
        np.random.seed(seed)

    t = psr.toas()
    fref = 1400
    v = (psr.freqs / fref)**idx

    minx, maxx = np.min(t), np.max(t)
    x = (t - minx) / (maxx - minx)
    T = (day / year) * (maxx - minx)

    size = 2 * components
    F = np.zeros((psr.nobs, size), "d")
    f = np.zeros(size, "d")

    for i in range(components):
        F[:, 2 * i] = np.cos(2 * math.pi * (i + 1) * x)
        F[:, 2 * i + 1] = np.sin(2 * math.pi * (i + 1) * x)

        f[2 * i] = f[2 * i + 1] = (i + 1) / T

    norm = A**2 * year**2 / (12 * math.pi**2 * T) 
    prior = norm * f ** (-gamma)

    y = np.sqrt(prior) * np.random.randn(size)
    psr.stoas[:] += (1.0 / day) * v * np.dot(F, y)


### Chromatic noise
def add_ch(psr, A, gamma, idx=-4, components=30, seed=None):
    """Add chromatic noise variations with P(f) = A^2 / (12 pi^2) (f year)^-gamma,
    using `components` Fourier bases.
    Optionally take a pseudorandom-number-generator seed."""    

    if seed is not None:
        np.random.seed(seed)

    t = psr.toas()
    fref = 1400
    v = (psr.freqs / fref)**idx

    minx, maxx = np.min(t), np.max(t)
    x = (t - minx) / (maxx - minx)
    T = (day / year) * (maxx - minx)

    size = 2 * components
    F = np.zeros((psr.nobs, size), "d")
    f = np.zeros(size, "d")

    for i in range(components):
        F[:, 2 * i] = np.cos(2 * math.pi * (i + 1) * x)
        F[:, 2 * i + 1] = np.sin(2 * math.pi * (i + 1) * x)

        f[2 * i] = f[2 * i + 1] = (i + 1) / T

    norm = A**2 * year**2 / (12 * math.pi**2 * T) 
    prior = norm * f ** (-gamma)

    y = np.sqrt(prior) * np.random.randn(size)
    psr.stoas[:] += (1.0 / day) * v * np.dot(F, y)

### Solar wind
def add_sw(psr, A, gamma, components=30, seed=None):
    """Add Solar wind variations with P(f) = A^2 / (12 pi^2) (f year)^-gamma,
    using `components` Fourier bases. 
    The electron density at 1 AU distance (nearth) is 1.0
    Optionally take a pseudorandom-number-generator seed."""    

    if seed is not None:
        np.random.seed(seed)

    t = psr.toas()

    theta ,R_earth = theta_impact(psr.earth_ssb, psr.sun_ssb, psr.psrPos)
    dt_sol_wind = (4.148808e3 /(psr.freqs**2))*dm_solar(1.0, theta, R_earth)

    minx, maxx = np.min(t), np.max(t)
    x = (t - minx) / (maxx - minx)
    T = (day / year) * (maxx - minx)

    size = 2 * components
    F = np.zeros((psr.nobs, size), "d")
    f = np.zeros(size, "d")

    for i in range(components):
        F[:, 2 * i] = np.cos(2 * math.pi * (i + 1) * x)
        F[:, 2 * i + 1] = np.sin(2 * math.pi * (i + 1) * x)

        f[2 * i] = f[2 * i + 1] = (i + 1) / T

    norm = A**2 * year**2 / (12 * math.pi**2 * T)
    prior = norm * f ** (-gamma)

    y = np.sqrt(prior) * np.random.randn(size)
    psr.stoas[:] += (1.0 / day) * dt_sol_wind * np.dot(F, y)

def theta_impact(earthssb,sunssb,pos_t):
    """Computes the solar impact angle"""
    earth = earthssb[:,:3]
    sun = sunssb[:,:3]
    pulsar = pos_t[:,:3]
    earthsun = earth - sun
    R_earth = np.sqrt(np.einsum('ij,ij->i', earthsun, earthsun))
    Re_cos_theta_impact = np.einsum('ij,ij->i', earthsun, pulsar)
    theta = np.arccos(-Re_cos_theta_impact / R_earth)
    return theta, R_earth

def _dm_solar_close(n_earth, r_earth):
    return (n_earth * AU_light_sec * AU_pc / r_earth)

def _dm_solar(n_earth, theta, r_earth):
    return ((np.pi - theta) *
            (n_earth * AU_light_sec * AU_pc
             / (r_earth * np.sin(theta))))

def dm_solar(n_earth, theta, r_earth):
    """
    Calculates Dispersion measure due to 1/r^2 solar wind density model.
    ::param :n_earth Solar wind proto/electron density at Earth (1/cm^3)
    ::param :theta: angle between sun and line-of-sight to pulsar (rad)
    ::param :r_earth :distance from Earth to Sun in (light seconds).
    See You et al. 2007 for more details.
    """
    return np.where(np.pi - theta >= 1e-5,
                    _dm_solar(n_earth, theta, r_earth),
                    _dm_solar_close(n_earth, r_earth))


## Jumps
## Poisson arrival times rate in jumps per year
def next_arrival(rate):
    rate = rate/365 # a jump a year
    return -math.log(1.0 - random.random()) / rate

## Generates jumps amplitudes. Max is given in microseconds
def jumps_amplitudes(num_jumps, max_size_ns):
    max_size = max_size_ns/2
    jumps = np.zeros(num_jumps)
    for j in range(num_jumps):
        size = round(random.uniform(0, max_size), 3)
        sign = random.choice([-1, 1])
        jumps[j] = size * sign
    return jumps

## Calculate the maximum time span of the simulated TOAs
def calculate_max_time(psrnames, datadir):
    """Calculate the maximum time span of the simulated TOAs
    by taking the maximum of the time span of the pulsars."""
    times = []
    for psrname in psrnames:
        parfile = (datadir + '/' + psrname + '.par')
        timfile = (datadir + '/' + psrname + '.tim')
        psr = T.tempopulsar(parfile=parfile,timfile=timfile)
        start = int((psr['START'].val) // 1)
        end = int((psr['FINISH'].val) // 1)
        tot_days = end - start
        times.append(tot_days)
    max_time = max(times)
    return max_time, start

## Create jumps
def make_jumps(max_time, start, max_size_ns, rate, seed=None):
    """create jumps taking one simulated TOAs (for the time).
    Optionally take a pseudorandom-number-generator seed."""

    if seed is not None:
        np.random.seed(seed)

    n = 0
    mjds = []
    while n < max_time:
        arrival = next_arrival(rate)
        n += arrival
        mjds.append(start+n)
    mjds.pop()
    num_jumps = len(mjds)
    jumps_ampls = jumps_amplitudes(num_jumps, max_size_ns)*1e-9
    print(len(mjds))
    return mjds, jumps_ampls

## Add jumps
def add_jumps(mjds, jumps, psr, freqs=None):
    """Add jumps to the simulated TOAs.
    Optionally take a pseudorandom-number-generator seed."""
    toas = psr.stoas
    indexes = []
    for i in mjds:
        # ff the value is in the array
        if i in toas:
            index = np.argwhere(toas == i)[:][:]
        else:
        # Find the largest value in the array that is less than the given value
            closest_val = toas[toas < i].max()
            if freqs == None:
                index = np.argwhere(toas >= closest_val) 
            else :
                freqmin, freqmax = freqs
                index = np.argwhere((toas >= closest_val) * (psr.freqs >= freqmin) * (psr.freqs <= freqmax))[:][:]
        indexes.append(index)
    for i, idx in enumerate(indexes):
        toas[idx] += jumps[i]/86400
    return toas


##################################
### Create new .par and .tim files
##################################

dm_gammas = {
    "J0030+0451":3.4,
    "J0613-0200":5.9,
    "J0900-3144":3.4,
    "J1024-0719":3.1,
    "J1446-4701":3.1,
    "J1603-7202":2.9,
    "J1730-2304":2.3,
    "J1902-5105":3.4,
    "J1939+2134":6.2,
    "J2145-0750":4.2,
    "J0125-2327":3.3,
    "J0614-3329":3.2,
    "J1017-7156":3.2,
    "J1045-4509":1.4,
    "J1545-4550":3.3,
    "J1643-1224":0.6,
    "J1832-0836":3.2,
    "J1909-3744":4.3,
    "J2124-3358":4.7,
    "J2241-5236":3,
    "J0437-4715":3.3,
    "J0711-6830":1.2,
    "J1022+1001":3.2,
    "J1125-6014":3.9,
    "J1600-3053":3.4,
    "J1713+0747":2.9,
    "J1744-1134":2.3,
    "J1857+0943":4.9,
    "J1933-6211":3.4,
    "J2129-5721":3.4
}

dm_log10As = {
    "J0030+0451":-16.3,
    "J0613-0200":-15.4,
    "J0900-3144":-16.5,
    "J1024-0719":-17.1,
    "J1446-4701":-17,
    "J1603-7202":-17.2,
    "J1730-2304":-16,
    "J1902-5105":-16.1,
    "J1939+2134":-14.6,
    "J2145-0750":-14.5,
    "J0125-2327":-16.9,
    "J0614-3329":-16.4,
    "J1017-7156":-16.1,
    "J1045-4509":-14.5,
    "J1545-4550":-16.6,
    "J1643-1224":-12.7,
    "J1832-0836":-16.9,
    "J1909-3744":-14,
    "J2124-3358":-14.9,
    "J2241-5236":-14.5,
    "J0437-4715":-14.3,
    "J0711-6830":-13.1,
    "J1022+1001":-16.6,
    "J1125-6014":-14.2,
    "J1600-3053":-15.9,
    "J1713+0747":-17.5,
    "J1744-1134":-15.9,
    "J1857+0943":-14.7,
    "J1933-6211":-16.3,
    "J2129-5721":-16.6
}

tn_log10As = {
    "J0030+0451":-16.8,
    "J0125-2327":-13.4,
    "J0437-4715":-13.48,
    "J0613-0200":-13.6,
    "J0614-3329":-13.8,
    "J0711-6830":-14.1,
    "J0900-3144":-12.7,
    "J1017-7156":-12.89,
    "J1022+1001":-13.8,
    "J1024-0719":-13.9,
    "J1045-4509":-12.38,
    "J1125-6014":-13.2,
    "J1446-4701":-16.9,
    "J1545-4550":-13.7,
    "J1600-3053":-13.2,
    "J1603-7202":-13.2,
    "J1643-1224":-12.9,
    "J1713+0747":-13.9,
    "J1730-2304":-13.5,
    "J1744-1134":-14.2,
    "J1832-0836":-13.5,
    "J1857+0943":-13.3,
    "J1902-5105":-13.1,
    "J1909-3744":-13.66,
    "J1933-6211":-16.9,
    "J1939+2134":-12.91,
    "J2124-3358":-17.3,
    "J2129-5721":-13.7,
    "J2145-0750":-13.5,
    "J2241-5236":-14
}

tn_gammas = {
    "J0030+0451":3.4,
    "J0125-2327":3.2,
    "J0437-4715":2.5,
    "J0613-0200":2.4,
    "J0614-3329":5,
    "J0711-6830":3.2,
    "J0900-3144":2,
    "J1017-7156":2.3,
    "J1022+1001":2.5,
    "J1024-0719":3.5,
    "J1045-4509":2.9,
    "J1125-6014":3.6,
    "J1446-4701":2.7,
    "J1545-4550":4.4,
    "J1600-3053":2.3,
    "J1603-7202":2.3,
    "J1643-1224":2.3,
    "J1713+0747":0.3,
    "J1730-2304":2.5,
    "J1744-1134":3.2,
    "J1832-0836":4.5,
    "J1857+0943":2.4,
    "J1902-5105":1.4,
    "J1909-3744":2,
    "J1933-6211":3.32,
    "J1939+2134":2.8,
    "J2124-3358":3.02,
    "J2129-5721":3.1,
    "J2145-0750":1.8,
    "J2241-5236":2.7
}

max_time, start = calculate_max_time(psrnames, datadir)
mjds, jumps = make_jumps(max_time, start, 200, 10, seed=None)

for i in range(100):
    savedir = f'/fred/oz005/users/vdimarco/Choosing_Noise_Models_in_PTAs/sims_sanity_check/sim_{i}/'
    print("saving in ", savedir)
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    for psrname in psrnames:
        parfile = (datadir + '/' + psrname + '.par')
        timfile = (datadir + '/' + psrname + '.tim')
        psr = T.tempopulsar(parfile=parfile,timfile=timfile)  ##### here timfile has simulated TOAs

        #### Add red noise
        tn_log10_A = tn_log10As[psrname]
        tn_gamma = tn_gammas[psrname]
        LT.add_rednoise(psr, 10**tn_log10_A, tn_gamma)

        # ### Add chromatic noise in all pulsars
        ch_log10_A = -13.5
        ch_gamma = 2.5
        add_ch(psr,10**(ch_log10_A),ch_gamma)

        #### Add solar wind
        # sw_gamma = 1.6
        # sw_log10_A = -6.0
        # add_sw(psr,10**(sw_log10_A),sw_gamma)

        ### Add dispersion measure noise
        dm_log10_A = dm_log10As[psrname]
        dm_gamma = dm_gammas[psrname]
        add_dm(psr, 10**(dm_log10_A), dm_gamma)

        #### Add jumps
        # new_sto = add_jumps(mjds, jumps, psr, freqs=(3000,3200))
        # for i in range(new_sto.size):
        #     psr.stoas[i] = new_sto[i]

        ### Save par and tim files
        psr.savepar(f"{savedir}/{psrname}.par")
        psr.savetim(f"{savedir}/{psrname}.tim")


