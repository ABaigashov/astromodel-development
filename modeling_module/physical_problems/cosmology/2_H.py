# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 20:46:40 2018

@author: Artyom1982
"""

import math
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec

from Graphic import plotChain

from scipy import *
from scipy.integrate import *

def hubble(Omega0,epsilon,x):
    res=np.sqrt((1-Omega0-epsilon)*(1+x)**3+Omega0+epsilon*(1+x)**6)
    return res

def likelihood(y):
    A=0
    B=0
    C=0
    Omega1=y[:,0]
    epsilon1=y[:,1]
    data = np.loadtxt("data/HUBBLE.txt", delimiter='\t', dtype=np.float)

    for i in range(21):
        A=A+(data[i,1])**2/(data[i,2]**2)
        B=B+data[i,1]*hubble(Omega1,epsilon1,data[i,0])/(data[i,2]**2)
        C=C+(hubble(Omega1,epsilon1,data[i,0]))**2/(data[i,2]**2)

    res=np.exp(-0.5*(A-B**2/C-13.62))
    return res

param_names = ['y1', 'y2']
x_initial = np.array([0.69,0.01])
n_params = len(param_names)
from colossus.utils import mcmc
walkers = mcmc.initWalkers(x_initial, nwalkers=200, random_seed = 156)
chain_thin, chain_full, _ = mcmc.runChain(likelihood, walkers)

mcmc.analyzeChain(chain_full, param_names = param_names,percentiles=[68.27, 95.45])


plotChain(chain_full, param_names)
plt.show()
