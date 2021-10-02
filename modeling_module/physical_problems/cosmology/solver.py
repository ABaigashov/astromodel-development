import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.tri as mtri
from sympy import symbols, sympify
import matplotlib.pyplot as plt
from scipy import *
from scipy.integrate import *


# from mpl_toolkits.mplot3d import *
# from matplotlib import cm
# import matplotlib.animation as animation
############################
#####УРАВНЕНИЕ ПУАССОНА#####
############################
path = 'modeling_module/physical_problems/cosmology/'

def Friedmann_eqs(ICS,a,eq):
    t, H, rho_m, rho_r, rho_d = ICS
    rho_s = symbols('rho_s')
    p = eq.subs(rho_s, rho_d)
    dt_da = 1/(H*a)
    dH_da = -(1/(H*a))*(1/2)*(rho_d+rho_m+(4/3)*rho_r+p)

    drho_m_da = -(3/a)*rho_m
    drho_r_da = -(4/a)*rho_r
    drho_d_da = -(3/a)*(rho_d+p)

    return dt_da, dH_da, drho_m_da, drho_r_da, drho_d_da

def Friedmann_eqs_2(ICS,a,eq):
    a, H, rho_m, rho_r, rho_d = ICS
    rho_s = symbols('rho_s')
    p = eq.subs(rho_s, rho_d)
    da_dt = H*a
    dH_dt = -(1/2)*(rho_d+rho_m+(4/3)*rho_r+p)

    drho_m_dt = -(3*H)*rho_m
    drho_r_dt = -(4*H)*rho_r
    drho_d_dt = -(3*H)*(rho_d+p)

    return da_dt, dH_dt, drho_m_dt, drho_r_dt, drho_d_dt

class Task_maker():

    def __init__(self, config):
        self.config = config

        self.omega_m = float(self.config.omega_m)
        self.omega_r = float(self.config.omega_r)
        self.z_max = float(self.config.z_max)
        self.t_max = float(self.config.t_max)
        self.H_0 = float(self.config.H_0)
        self.dark_components = []
        self.omega_d = []
        self.name_SNE = self.config.name_SNE
        self.name_Hubble = self.config.name_Hubble
        self.row_SNE = [int(self.config.row_1),int(self.config.row_2),int(self.config.row_3)]
        self.row_Hubble = [int(self.config.row_4),int(self.config.row_5),int(self.config.row_6)]

        if self.config.dark_components:
            for compts in self.config.dark_components:
                if compts.equation_d:
                    self.dark_components.append((compts.equation_d))
                if compts.omega_d:
                    self.omega_d.append(float(compts.omega_d))

class Cosmology_data():

    def __init__(self, name, rows):
        self.name = name
        self.rows = rows
        self.z0 = []
        self.parameter = []
        self.err = []

    def Data_loader(self):
        name = path + "data/" + self.name
        data = np.loadtxt(name, delimiter='\t', dtype=np.float)
        self.z0 = data[:,self.rows[0]-1]
        self.parameter = data[:,self.rows[1]-1]
        self.err = data[:,self.rows[2]-1]

class Cosmology_calculus():

    def __init__(self, Task, Data_1, Data_2):
        self.omega_m = Task.omega_m
        self.omega_r = Task.omega_r
        self.omega_d = Task.omega_d[0]
        self.redshifts_1 = Data_1.z0
        self.redshifts_2 = Data_2.z0
        self.EOS = sympify(Task.dark_components[0])
        self.H_0 = Task.H_0
        self.mu_m = []
        self.H_m = []
        self.mu_i = []
        self.H_i = []
        self.z_i = []
        self.DA = []
        self.mu_o = Data_1.parameter
        self.H_o = Data_2.parameter
        self.Zmax = Task.z_max
        self.Tmax = Task.t_max/13.6
        self.err_mu = Data_1.err
        self.err_H = Data_2.err
        self.scale_factor = []
        self.HUBBLE = []

    def mu_versus_z(self):
        DL = np.zeros(len(self.redshifts_1))
        dist = np.zeros(len(DL))
        for i in range(len(DL)):
            dist0=0
            a = 1/(self.redshifts_1[i]+1)
            N = int((1-a)//0.01)
            if N<2:
                N=2
            scale = np.linspace(1,a,N)
            ICS = 0, 1, 3*self.omega_m, 3*self.omega_r, 3*self.omega_d
            sol = odeint(Friedmann_eqs, ICS, scale, args = (self.EOS,))
            for j in range(1,N):
                #dist[i]=dist0 - (scale[j+1]-scale[j])/(sol[j,1]*scale[j]**2)
                dist[i]=dist0 - 2*(sol[j,0]-sol[j-1,0])/(scale[j]+scale[j-1])
                dist0 = dist[i]
            DL[i] = (1+self.redshifts_1[i])*dist[i]
            self.mu_m.append(5*np.log10(DL[i]))

    def mu_diagram(self):
        DL = np.zeros(200)
        dist = np.zeros(len(DL))
        for i in range(len(DL)):
            dist0=0
            z = 0.0001 + self.Zmax*i/200
            a = 1/(z+1)
            N = int((1-a)//0.01)
            if N<2:
                N=2
            scale = np.linspace(1,a,N)
            ICS = 0, 1, 3*self.omega_m, 3*self.omega_r, 3*self.omega_d
            sol = odeint(Friedmann_eqs, ICS, scale, args = (self.EOS,))
            for j in range(1,N):
                #dist[i]=dist0 - (scale[j+1]-scale[j])/(sol[j,1]*scale[j]**2)
                dist[i]=dist0 - 2*(sol[j,0]-sol[j-1,0])/(scale[j]+scale[j-1])
                dist0 = dist[i]
            self.H_i.append(self.H_0*sol[N-1,1])
            DL[i] = (1+z)*dist[i]
            self.DA.append(DL[i]/(1+z)**2)
            self.z_i.append(z)
            self.mu_i.append(5*np.log10(DL[i])-5*np.log10(self.H_0/100)+43.16-0.713)

    def integration(self):
        for i in range(201):
            t = 0.000 + self.Tmax*i/200
            times = np.linspace(0,t,2)
            ICS = 1, 1, 3*self.omega_m, 3*self.omega_r, 3*self.omega_d
            sol = odeint(Friedmann_eqs_2, ICS, times, args = (self.EOS,))
            self.scale_factor.append(sol[1,0])
            self.HUBBLE.append(sol[1,1]*self.H_0)

    def hubble_versus_z(self):
        z = np.zeros(len(self.redshifts_2))
        for i in range(len(z)):
            a = 1/(self.redshifts_2[i]+1)
            N=2
            scale = np.linspace(1,a,N)
            ICS = 0, 1, 3*self.omega_m, 3*self.omega_r, 3*self.omega_d
            sol = odeint(Friedmann_eqs, ICS, scale, args = (self.EOS,))
            self.H_m.append(sol[N-1,1])

    def visualization(self):
        fig = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
        ax = fig.add_subplot(111)
        fig2 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
        ax2 = fig2.add_subplot(111)
        fig3 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
        ax3 = fig3.add_subplot(111)

        ax.plot(self.z_i,self.mu_i)
        ax.errorbar(self.redshifts_1,self.mu_o, yerr=self.err_mu, fmt=".")

        ax.minorticks_on()
        ax.set_title('SNe_Ia_magnitude_vs_redshift_diagram')
        ax.set_xlabel('z')
        ax.set_ylabel('apparent_magnitude')
        ax.grid(which='major',linewidth = 2)
        ax.grid(which='minor')
        fig.savefig(path + 'results/SNe_Ia')

        ax2.plot(self.z_i,self.H_i)
        ax2.errorbar(self.redshifts_2,self.H_o, yerr=self.err_H,fmt=".")
        ax2.minorticks_on()
        ax2.set_title('Hubble_parameter_vs_redshift')
        ax2.set_xlabel('z')
        ax2.set_ylabel('H, km/s/Mpc')
        ax2.grid(which='major',linewidth = 2)
        ax2.grid(which='minor')
        fig2.savefig(path + 'results/H(z)')

        ax3.plot(self.scale_factor,self.HUBBLE)
        ax3.minorticks_on()
        ax3.set_title('Angular_diameter_distance')
        ax3.set_xlabel('z')
        ax3.set_ylabel('angular_diameter_distance')
        ax3.grid(which='major',linewidth = 2)
        ax3.grid(which='minor')
        fig3.savefig(path + 'results/DA(z)')
