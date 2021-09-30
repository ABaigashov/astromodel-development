# -*- coding: utf-8 -*-

import math
import numpy as np

import warnings
warnings.filterwarnings("ignore")

from scipy import *
from scipy.integrate import *

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec










def chi_square(mu,mu0,err):
    A0,B0,C0 = 0,0,0
    for i in range(1,len(mu0)):
        A = A0 + (mu[i] - mu0[i])**2/err[i]**2
        B = B0 + (mu[i] - mu0[i])/err[i]**2
        C = C0 + 1/err[i]**2
        C0 = C
        B0 = B
        A0 = A
    return  A-B**2/C, 10**((B/C+42.384)/5)

def chi_square_2(H,H0,err):
    A0,B0,C01 = 0,0,0
    for i in range(len(H0)):
        A = A0 + H0[i]**2/err[i]**2
        B = B0 + H0[i]*H[i]/err[i]**2
        C1 = C01 + H[i]**2/err[i]**2
        C01 = C1
        B0 = B
        A0 = A
    return  A-B**2/C1, B/C1


def definition_opt(Omega_min, chi_max, a1, a2, a3):
    Omega_d=Omega_min
    w0 = a1
    c0 = a3
    zeta0 = a2
    #w0=y[:,1]
    chi_opt = chi_max

    data = np.loadtxt("SNE580.txt", delimiter='\t', dtype=np.float)
    z0 = data[:,0]
    mu0 = data[:,1]
    err = data[:,2]

    for k in range(600):
        Omega = Omega_d +0.001*k
        w = w0
        res0 = mu_versus_z(z0,Omega,w,zeta0,c0)[1]
        res2=chi_square(res0,mu0,err)
        res1 = res2[0]
        H_opt = res2[1]
        if res1<chi_opt:
            chi_opt=res1
            Omega_opt = Omega

    return chi_opt, Omega_opt, H_opt

def definition_opt_2(Omega_min, chi_max, a1, a2, a3):
    Omega_d=Omega_min
    w0 = a1
    zeta0 = a2
    c0 = a3
    #w0=y[:,1]
    chi_opt = chi_max

    data_2 = np.loadtxt("data/HUBBLE3.txt", delimiter=' ', dtype=np.float)
    z1 = data_2[:,0]
    H0 = data_2[:,1]
    err1 = data_2[:,2]

    for k in range(600):
        Omega = Omega_d +0.001*k
        w = w0
        res0 = hubble_versus_z(z1,Omega,w,zeta0,c0)[1]
        res2=chi_square_2(res0,H0,err1)
        res1 = res2[0]
        H_opt = res2[1]
        if res1<chi_opt:
            chi_opt=res1
            Omega_opt = Omega

    return chi_opt, Omega_opt, H_opt

def definition_int(Omega_opt, chi_opt, delta1,a1,a2,a3):
    Omega_d=Omega_opt
    w0 = a1
    zeta0 = a2
    c0 = a3
    #w0=y[:,1]
    chi = chi_opt

    data = np.loadtxt("data/SNE580.txt", delimiter='\t', dtype=np.float)
    z0 = data[:,0]
    mu0 = data[:,1]
    err = data[:,2]

    for k in range(600):
        Omega = Omega_d + 0.001*k
        w = w0
        res0 = mu_versus_z(z0,Omega,w,zeta0,c0)[1]
        res1=chi_square(res0,mu0,err)[0]
        if res1>chi+delta1:
            Omega_max = Omega-0.001
            Omega1 = Omega+0.001
            res0 = mu_versus_z(z0,Omega1,w,zeta0,c0)[1]
            res1=chi_square(res0,mu0,err)[0]
            if res1>chi+delta1:
                Omega2 = Omega1+0.001
                res0 = mu_versus_z(z0,Omega2,w,zeta0,c0)[1]
                res1=chi_square(res0,mu0,err)[0]
                if res1>chi+delta1:
                    break

    for k in range(600):
        Omega = Omega_opt - 0.001*k
        w = w0
        res0 = mu_versus_z(z0,Omega,w,zeta0,c0)[1]
        res1=chi_square(res0,mu0,err)[0]
        if res1>chi+delta1:
            Omega_min = Omega-0.001
            Omega1 = Omega-0.001
            res0 = mu_versus_z(z0,Omega1,w,zeta0,c0)[1]
            res1=chi_square(res0,mu0,err)[0]
            if res1>chi+delta1:
                Omega2 = Omega1-0.001
                res0 = mu_versus_z(z0,Omega2,w,zeta0,c0)[1]
                res1=chi_square(res0,mu0,err)[0]
                if res1>chi+delta1:
                    break

    return Omega_max, Omega_min

def definition_int_2(Omega_opt, chi_opt, delta1,a1,a2,a3):
    Omega_d=Omega_opt
    w0 = a1
    zeta0 = a2
    c0 = a3

    #w0=y[:,1]
    chi = chi_opt

    data_2 = np.loadtxt("data/HUBBLE3.txt", delimiter=' ', dtype=np.float)
    z1 = data_2[:,0]
    H0 = data_2[:,1]
    err1 = data_2[:,2]

    for k in range(600):
        Omega = Omega_d +0.001*k
        w = w0
        res0 = hubble_versus_z(z1,Omega,w,zeta0,c0)[1]
        res1=chi_square_2(res0,H0,err1)[0]
        if res1>chi+delta1:
            Omega_max = Omega-0.001
            Omega1 = Omega+0.001
            res0 = hubble_versus_z(z1,Omega1,w,zeta0,c0)[1]
            res1=chi_square_2(res0,H0,err1)[0]
            if res1>chi+delta1:
                Omega2 = Omega1+0.001
                res0 = hubble_versus_z(z1,Omega2,w,zeta0,c0)[1]
                res1=chi_square_2(res0,H0,err1)[0]
                if res1>chi+delta1:
                    break

    for k in range(600):
        Omega = Omega_opt - 0.001*k
        w = w0
        res0 = hubble_versus_z(z1,Omega,w,zeta0,c0)[1]
        res1=chi_square_2(res0,H0,err1)[0]
        if res1>chi+delta1:
            Omega_min = Omega-0.001
            Omega1 = Omega-0.001
            res0 = hubble_versus_z(z1,Omega1,w,zeta0,c0)[1]
            res1=chi_square_2(res0,H0,err1)[0]
            if res1>chi+delta1:
                Omega2 = Omega1-0.001
                res0 = hubble_versus_z(z1,Omega2,w,zeta0,c0)[1]
                res1=chi_square_2(res0,H0,err1)[0]
                if res1>chi+delta1:
                    break

    return Omega_max, Omega_min

w=-1
zeta0=0
c0=0


#w=-1
#zeta0 = 0
#c0=0

#w=-1
#zeta0 = 0.0
#c0=0.0

OmegaD0 = 0.65

result = hubble_versus_z(z1,OmegaD0,w,zeta0,c0)[1]

res1=chi_square_2(result,H0,err1)[0]

print(res1)

#w=-1
#zeta0 = 0
#c0=0

#result = definition_opt_2(0.65, 100, w, zeta0, c0)
##
#print(result)
#
##
#result = definition_int_2(result[1],result[0],1, w,zeta0,c0)
###
#print(result)


def likelihood(y):
    Omega_d=y[:,0]
    w0 = -1
    #w0=y[:,1]
    res = np.zeros(len(Omega_d))

    data = np.loadtxt("data/SNE580.txt", delimiter='\t', dtype=np.float)
    z0 = data[:,0]
    mu0 = data[:,1]
    err = data[:,2]

    for k in range(len(Omega_d)):
        Omega = Omega_d[k]
        #w = w0[k]
        w = w0
        res0 = mu_versus_z(z0,Omega,w,0,0)[1]
        rrr=chi_square(res0,mu0,err)[0]
        result=np.exp(-0.5*rrr+0.5*552.8887)
        res[k] = result


    return res

def likelihood2(y):
    Omega_d=y[:,0]
    w0 = -1
    #w0=y[:,1]
    res = np.zeros(len(Omega_d))

    data_2 = np.loadtxt("data/HUBBLE3.txt", delimiter=' ', dtype=np.float)
    z1 = data_2[:,0]
    H0 = data_2[:,1]
    err1 = data_2[:,2]

    for k in range(len(Omega_d)):
        Omega = Omega_d[k]
        #w = w0[k]
        w = w0
        res0 = hubble_versus_z(z1,Omega,w,0,0)[1]
        rrr=chi_square_2(res0,H0,err1)[0]
        result=np.exp(-0.5*rrr+0.5*19.282992)
        res[k] = result


    return res


param_names = ['y1', 'y2']

param_names = ['y1']

x_initial = np.array([0.724,-1])

x_initial = np.array([0.72])

n_params = len(param_names)

from colossus.utils import mcmc
walkers = mcmc.initWalkers(x_initial, initial_step=0.1, nwalkers=6, random_seed = 156)

chain_thin, chain_full, _ = mcmc.runChain(likelihood2, walkers)

mcmc.analyzeChain(chain_full, param_names = param_names,percentiles=[68.27, 95.45, 99.73])
