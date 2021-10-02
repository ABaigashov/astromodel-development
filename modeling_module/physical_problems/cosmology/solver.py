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

        self.name_SNE = config.name_SNE
        self.name_Hubble = config.name_Hubble
        self.plot_diagram_1 = config.task_1
        self.plot_diagram_2 = config.task_2
        self.plot_diagram_3 = config.task_3
        self.plot_diagram_4 = config.task_4
        self.plot_diagram_5 = config.task_5
        self.plot_diagram_6 = config.task_6
        self.plot_diagram_7 = config.task_7
        self.plot_diagram_8 = config.task_8

        self.dark_components = []
        self.omega_d = []
        self.row_SNE = [int(config.row_1),int(config.row_2),int(config.row_3)]
        self.row_Hubble = [int(config.row_4),int(config.row_5),int(config.row_6)]

        if config.dark_components:
            for compts in config.dark_components:
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

    def __init__(self, config, Task, Data_1, Data_2):
        self.omega_m = float(config.omega_m)
        self.omega_r = float(config.omega_r)
        self.z_max = float(config.z_max)
        self.t_max = float(config.t_max)/13.6*float(config.H_0)/70
        self.omega_d = Task.omega_d[0]
        self.redshifts_1 = Data_1.z0
        self.redshifts_2 = Data_2.z0
        self.EOS = sympify(Task.dark_components[0])
        self.H_0 = float(config.H_0)

        self.mu_m = []
        self.H_m = []
        self.mu_i = []
        self.H_i = []
        self.z_i = []
        self.t_i = []
        self.T_i = []
        self.DA = []
        self.scale_factor = []
        self.HUBBLE = []
        self.Omega_d = []
        self.Omega_m = []
        self.OMEGA_d = []
        self.OMEGA_m = []

        self.mu_o = Data_1.parameter
        self.H_o = Data_2.parameter
        self.err_mu = Data_1.err
        self.err_H = Data_2.err

    def mu_diagram(self):
        DL = np.zeros(200)
        dist = np.zeros(len(DL))
        for i in range(len(DL)):
            dist0=0
            z = 0.0001 + self.z_max*i/200
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
            DL[i] = (1+z)*dist[i]
            self.H_i.append(self.H_0*sol[N-1,1])
            self.t_i.append(sol[N-1,0]*13.6*self.H_0/70)
            self.DA.append(DL[i]/(1+z)**2)
            self.Omega_d.append(sol[N-1,2]/(3*sol[N-1,1]**2))
            self.Omega_m.append(sol[N-1,4]/(3*sol[N-1,1]**2))
            self.z_i.append(z)
            self.mu_i.append(5*np.log10(DL[i])-5*np.log10(self.H_0/100)+43.16-0.713)

    def integration(self):
        for i in range(201):
            t = 0.000 + self.t_max*i/200
            times = np.linspace(0,t,2)
            ICS = 1, 1, 3*self.omega_m,3*self.omega_r, 3*self.omega_d
            sol = odeint(Friedmann_eqs_2, ICS, times, args = (self.EOS,))
            self.scale_factor.append(sol[1,0])
            self.T_i.append(t*13.6*70/self.H_0)
            self.OMEGA_d.append(sol[1,4]/(3*sol[1,1]**2))
            self.OMEGA_m.append(sol[1,2]/(3*sol[1,1]**2))
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

    def visualization(self, Task):

        if Task.plot_diagram_1 == True:

            fig = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax = fig.add_subplot(111)
            ax.plot(self.z_i,self.mu_i)
            ax.errorbar(self.redshifts_1,self.mu_o, yerr=self.err_mu, fmt=".")
            ax.minorticks_on()
            ax.set_title('SNe_Ia_magnitude_vs_redshift_diagram')
            ax.set_xlabel('z')
            ax.set_ylabel('apparent_magnitude')
            ax.grid(which='major',linewidth = 2)
            ax.grid(which='minor')
            fig.savefig(path + 'results/SNe_Ia')

        if Task.plot_diagram_2 == True:

            fig2 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax2 = fig2.add_subplot(111)
            ax2.plot(self.z_i,self.H_i)
            ax2.errorbar(self.redshifts_2,self.H_o, yerr=self.err_H,fmt=".")
            ax2.minorticks_on()
            ax2.set_title('Hubble_parameter_vs_redshift')
            ax2.set_xlabel('z')
            ax2.set_ylabel('H, km/s/Mpc')
            ax2.grid(which='major',linewidth = 2)
            ax2.grid(which='minor')
            fig2.savefig(path + 'results/H(z)')

        if Task.plot_diagram_3 == True:

            fig3 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax3 = fig3.add_subplot(111)
            ax3.plot(self.z_i,self.t_i)
            ax3.minorticks_on()
            ax3.set_title('time_vs_redshift')
            ax3.set_xlabel('z')
            ax3.set_ylabel('t, Gyr')
            ax3.grid(which='major',linewidth = 2)
            ax3.grid(which='minor')
            fig3.savefig(path + 'results/Backlook_time')

        if Task.plot_diagram_4 == True:

            fig5 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax5 = fig5.add_subplot(111)
            ax5.plot(self.T_i,self.scale_factor)
            ax5.minorticks_on()
            ax5.set_title('scale_factor_in_future')
            ax5.set_xlabel('t, Gyr')
            ax5.set_ylabel('a')
            ax5.grid(which='major',linewidth = 2)
            ax5.grid(which='minor')
            fig5.savefig(path + 'results/a(t)')

        if Task.plot_diagram_5 == True:

            fig4 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax4 = fig4.add_subplot(111)
            ax4.plot(self.T_i,self.HUBBLE)
            ax4.minorticks_on()
            ax4.set_title('Hubble_parameter_in_future')
            ax4.set_xlabel('t, Gyr')
            ax4.set_ylabel('H, km/s/Mpc')
            ax4.grid(which='major',linewidth = 2)
            ax4.grid(which='minor')
            fig4.savefig(path + 'results/H(t)')

        if Task.plot_diagram_6 == True:

            fig5 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax5 = fig5.add_subplot(111)
            ax5.plot(self.z_i,self.Omega_d)
            ax5.plot(self.z_i,self.Omega_m)
            ax5.minorticks_on()
            ax5.set_title('Fraction_of_densities_in_past')
            ax5.set_xlabel('z')
            ax5.set_ylabel('Omega_d, Omega_m')
            ax5.grid(which='major',linewidth = 2)
            ax5.grid(which='minor')
            fig5.savefig(path + 'results/Omega_in_past')

        if Task.plot_diagram_7 == True:

            fig5 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax5 = fig5.add_subplot(111)
            ax5.plot(self.T_i,self.OMEGA_d)
            ax5.plot(self.T_i,self.OMEGA_m)
            ax5.minorticks_on()
            ax5.set_title('Fraction_of_densities_in_past')
            ax5.set_xlabel('t, Gyr')
            ax5.set_ylabel('Omega')
            ax5.grid(which='major',linewidth = 2)
            ax5.grid(which='minor')
            fig5.savefig(path + 'results/Omega_in_future')

        if Task.plot_diagram_7 == True:

            fig4 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax4 = fig4.add_subplot(111)
            ax4.plot(self.z_i,self.DA)
            ax4.minorticks_on()
            ax4.set_title('Angular_diameter_distance')
            ax4.set_xlabel('z')
            ax4.set_ylabel('angular_diameter_distance')
            ax4.grid(which='major',linewidth = 2)
            ax4.grid(which='minor')
            fig4.savefig(path + 'results/DA(z)')
