import numpy as np
from sympy import symbols, sympify
import matplotlib.pyplot as plt
from scipy import *
from scipy.integrate import *
import os
import shutil
import zipfile
from urllib import request

# from mpl_toolkits.mplot3d import *
# from matplotlib import cm
# import matplotlib.animation as animation
############################
#####УРАВНЕНИЕ ПУАССОНА#####
############################
path = 'modeling_module/physical_problems/cosmology/'

def Friedmann_eqs(ICS,a,eq):
    t, H, rho_m, rho_r, rho_d, rho_k = ICS
    rho_s = symbols('rho_s')
    p = eq.subs(rho_s, rho_d)
    dt_da = 1/(H*a)
    dH_da = -(1/(H*a))*(1/2)*(rho_d+rho_m+(4/3)*rho_r+p+(2/3)*rho_k)

    drho_m_da = -(3/a)*rho_m
    drho_r_da = -(4/a)*rho_r
    drho_d_da = -(3/a)*(rho_d+p)
    drho_k_da = -(2/a)*rho_k

    return dt_da, dH_da, drho_m_da, drho_r_da, drho_d_da, drho_k_da

def Friedmann_eqs_2(ICS,a,eq):
    a, H, rho_m, rho_r, rho_d, rho_k = ICS
    rho_s = symbols('rho_s')
    p = eq.subs(rho_s, rho_d)
    da_dt = H*a
    dH_dt = -(1/2)*(rho_d+rho_m+(4/3)*rho_r+p+2*rho_k/3)

    drho_m_dt = -(3*H)*rho_m
    drho_r_dt = -(4*H)*rho_r
    drho_d_dt = -(3*H)*(rho_d+p)
    drho_k_dt = -(2*H)*rho_k

    return da_dt, dH_dt, drho_m_dt, drho_r_dt, drho_d_dt, drho_k_dt

class Cosmology_data():

    def __init__(self, name, rows):
        self.name = name
        self.rows = rows
        self.z0 = []
        self.parameter = []
        self.err = []

    def Data_loader(self):
        remote_url = self.name
        name_file = path + "data/" + "1.txt"
        request.urlretrieve(remote_url, name_file)
        data = np.loadtxt(name_file, delimiter='\t', dtype=np.float)
        self.z0 = data[:,self.rows[0]-1]
        self.parameter = data[:,self.rows[1]-1]
        self.err = data[:,self.rows[2]-1]

class Cosmology_calculus():

    def __init__(self, config, model, Data_1, Data_2):

        self.omega_m = float(model.omega_m)
        self.omega_r = float(model.omega_r)
        self.omega_d = float(model.omega_d)
        self.omega_k = 1 - self.omega_m - self.omega_r - self.omega_d
        self.EOS = sympify(model.equation_d)
        self.title_of_model = model.title_of_model

        self.z_max = float(config.z_max)
        self.t_max = float(config.t_max)/13.6*float(config.H_0)/70
        self.H_0 = float(config.H_0)
        self.task_1 = config.task_1

        self.redshifts_1 = Data_1.z0
        self.redshifts_2 = Data_2.z0

        self.mu_m = []
        self.H_m = []

        self.mu_i = []
        self.H_i = []
        self.t_i = []
        self.z_i = []

        self.DA = []
        self.Omega_d = []
        self.Omega_m = []

        self.T_i = []
        self.scale_factor = []
        self.HUBBLE = []
        self.OMEGA_d = []
        self.OMEGA_m = []

        self.mu_o = Data_1.parameter
        self.H_o = Data_2.parameter
        self.err_mu = Data_1.err
        self.err_H = Data_2.err

        self.chi_square_mu = 10000
        self.chi_square_H = 10000
        self.H_opt = 72

    def mu_diagram(self):
        DL = np.zeros(200)
        DM = np.zeros(200)
        dist = np.zeros(len(DL))
        for i in range(len(DL)):
            dist0=0
            z = 0.005 + self.z_max*i/200
            a = 1/(z+1)
            N = int((1-a)//0.01)
            if N<2:
                N=2
            scale = np.linspace(1,a,N)
            ICS = 0, 1, 3*self.omega_m, 3*self.omega_r, 3*self.omega_d, 3*self.omega_k
            sol = odeint(Friedmann_eqs, ICS, scale, args = (self.EOS,))
            if self.task_1 == True:
                for j in range(1,N):
                    dist[i]=dist0 - 2*(sol[j,0]-sol[j-1,0])/(scale[j]+scale[j-1])
                    dist0 = dist[i]
                if self.omega_k>0:
                    DM[i] = np.sinh(np.sqrt(self.omega_k)*dist[i])/np.sqrt(self.omega_k)
                if self.omega_k<0:
                    DM[i] = np.sin(np.sqrt(-self.omega_k)*dist[i])/np.sqrt(-self.omega_k)
                if self.omega_k==0:
                    DM[i] = dist[i]
                DL[i] = (1+z)*DM[i]
                self.DA.append(DL[i]*(13.6/3.2616)*(70/self.H_0)/(1+z)**2)
                self.mu_i.append(5*np.log10(DL[i])-5*np.log10(self.H_0/100)+43.16-0.713)
            self.H_i.append(self.H_0*sol[N-1,1])
            self.t_i.append(sol[N-1,0]*13.6*70/self.H_0)
            self.Omega_d.append(sol[N-1,2]/(3*sol[N-1,1]**2))
            self.Omega_m.append(sol[N-1,4]/(3*sol[N-1,1]**2))
            self.z_i.append(z)


    def integration(self):
        for i in range(201):
            t = 0.000 + self.t_max*i/200
            times = np.linspace(0,t,2)
            ICS = 1, 1, 3*self.omega_m,3*self.omega_r, 3*self.omega_d, 3*self.omega_k
            sol = odeint(Friedmann_eqs_2, ICS, times, args = (self.EOS,))
            self.scale_factor.append(sol[1,0])
            self.T_i.append(t*13.6*70/self.H_0)
            self.OMEGA_d.append(sol[1,2]/(3*sol[1,1]**2))
            self.OMEGA_m.append(sol[1,4]/(3*sol[1,1]**2))
            self.HUBBLE.append(sol[1,1]*self.H_0)

    def hubble_versus_z(self):
        z = np.zeros(len(self.redshifts_2))
        for i in range(len(z)):
            a = 1/(self.redshifts_2[i]+1)
            N=2
            scale = np.linspace(1,a,N)
            ICS = 0, 1, 3*self.omega_m, 3*self.omega_r, 3*self.omega_d, 3*self.omega_k
            sol = odeint(Friedmann_eqs, ICS, scale, args = (self.EOS,))
            self.H_m.append(sol[N-1,1])

    def mu_versus_z(self):
        DL = np.zeros(len(self.redshifts_1))
        DM = np.zeros(len(self.redshifts_1))
        dist = np.zeros(len(DL))
        for i in range(len(DL)):
            dist0=0
            a = 1/(self.redshifts_1[i]+1)
            N = int((1-a)//0.01)
            if N<2:
                N=2
            scale = np.linspace(1,a,N)
            ICS = 0, 1, 3*self.omega_m, 3*self.omega_r, 3*self.omega_d, 3*self.omega_k
            sol = odeint(Friedmann_eqs, ICS, scale, args = (self.EOS,))
            for j in range(1,N):
                #dist[i]=dist0 - (scale[j+1]-scale[j])/(sol[j,1]*scale[j]**2)
                dist[i]=dist0 - 2*(sol[j,0]-sol[j-1,0])/(scale[j]+scale[j-1])
                dist0 = dist[i]
            if self.omega_k>0:
                DM[i] = np.sinh(np.sqrt(self.omega_k)*dist[i])/np.sqrt(self.omega_k)
            if self.omega_k<0:
                DM[i] = np.sin(np.sqrt(-self.omega_k)*dist[i])/np.sqrt(-self.omega_k)
            if self.omega_k==0:
                DM[i] = dist[i]
            DL[i] = (1+self.redshifts_1[i])*DM[i]
            self.mu_m.append(5*np.log10(DL[i])-5*np.log10(self.H_0/100)+43.16-0.713)

    def chi_square_magnitude(self):
        A0,B0,C0 = 0,0,0
        for i in range(0,len(self.mu_m)):
            A = A0 + (self.mu_m[i] - self.mu_o[i])**2/self.err_mu[i]**2
            B = B0 + (self.mu_m[i] - self.mu_o[i])/self.err_mu[i]**2
            C = C0 + 1/self.err_mu[i]**2
            C0 = C
            B0 = B
            A0 = A
        self.chi_square_mu = A-B**2/C
        self.H_opt = 70*10**((B/C+42.384)/5)/299752458

    def chi_square_hubble(self):
        A0,B0,C01 = 0,0,0
        for i in range(len(self.H_o)):
            A = A0 + self.H_o[i]**2/self.err_H[i]**2
            B = B0 + self.H_o[i]*self.H_m[i]/self.err_H[i]**2
            C1 = C01 + self.H_m[i]**2/self.err_H[i]**2
            C01 = C1
            B0 = B
            A0 = A
        self.chi_square_H = A-B**2/C1
        self.H_opt = B/C1

class Model_Var():

    def __init__(self, omega_d, omega_m, omega_r, equation_d, title):

        self.omega_d = omega_d
        self.omega_m = omega_m
        self.omega_r = omega_r
        self.equation_d = equation_d
        self.title_of_model = title



class Visualization:

    def __init__(self, models, output, job):

        self.models = models
        self.output = output
        self.job = job
        os.mkdir(self.output)

    def graphics(self, config):

        if config.task_1 == True:

            legends = []

            fig = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax = fig.add_subplot(111)
            for model in self.models:
                ax.plot(model.z_i,model.mu_i)
                legends.append(model.title_of_model)
            ax.errorbar(model.redshifts_1,model.mu_o, yerr=model.err_mu, fmt=".")
            legends.append('observational_data')
            ax.legend(legends)
            ax.minorticks_on()
            ax.set_title('SNe_Ia_magnitude_vs_redshift')
            ax.set_xlabel('z')
            ax.set_ylabel('M')
            ax.grid(which='major',linewidth = 2)
            ax.grid(which='minor')
            fig.savefig(f'{self.output}/SNe_Ia.png')
            self.job.progress += 1 / 8

        if config.task_2 == True:

            legends = []

            fig2 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax2 = fig2.add_subplot(111)
            for model in self.models:
                ax2.plot(model.z_i,model.H_i)
                legends.append(model.title_of_model)
            ax2.errorbar(model.redshifts_2,model.H_o,yerr=model.err_H,fmt=".")
            legends.append("observational data")
            ax2.legend(legends)
            ax2.minorticks_on()
            ax2.set_title('Hubble_parameter_vs_redshift')
            ax2.set_xlabel('z')
            ax2.set_ylabel('H, km/s/Mpc')
            ax2.grid(which='major',linewidth = 2)
            ax2.grid(which='minor')
            fig2.savefig(f'{self.output}/H(z)')
            self.job.progress += 1 / 8
        #
        if config.task_3 == True:
            legends = []
            fig3 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax3 = fig3.add_subplot(111)
            for model in self.models:
                ax3.plot(model.z_i,model.t_i)
                legends.append(model.title_of_model)
            ax3.legend(legends)
            ax3.minorticks_on()
            ax3.set_title('time_vs_redshift')
            ax3.set_xlabel('z')
            ax3.set_ylabel('t, Gyr')
            ax3.grid(which='major',linewidth = 2)
            ax3.grid(which='minor')
            fig3.savefig(f'{self.output}/Backlook_time')
            self.job.progress += 1 / 8
        #
        if config.task_4 == True:

            legends = []

            fig5 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax5 = fig5.add_subplot(111)
            for model in self.models:
                ax5.plot(model.T_i,model.scale_factor)
                legends.append(model.title_of_model)
            ax5.legend(legends)
            ax5.minorticks_on()
            ax5.set_title('scale_factor_in_future')
            ax5.set_xlabel('t, Gyr')
            ax5.set_ylabel('a')
            ax5.grid(which='major',linewidth = 2)
            ax5.grid(which='minor')
            fig5.savefig(f'{self.output}/a(t)')
            self.job.progress += 1 / 8
        #
        if config.task_5 == True:

            legends = []

            fig4 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax4 = fig4.add_subplot(111)
            for model in self.models:
                ax4.plot(model.T_i,model.HUBBLE)
                legends.append(model.title_of_model)
            ax4.legend(legends)
            ax4.minorticks_on()
            ax4.set_title('Hubble_parameter_in_future')
            ax4.set_xlabel('t, Gyr')
            ax4.set_ylabel('H, km/s/Mpc')
            ax4.grid(which='major',linewidth = 2)
            ax4.grid(which='minor')
            fig4.savefig(f'{self.output}/H(t)')
            self.job.progress += 1 / 8
        #
        if config.task_6 == True:
            legends = []
            fig5 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax5 = fig5.add_subplot(111)
            for model in self.models:
                ax5.plot(model.z_i,model.Omega_d)
                legends.append(model.title_of_model + ', fraction_of_dark_energy')
                ax5.plot(model.z_i,model.Omega_m)
                legends.append(model.title_of_model + ', fraction_of_matter')
            ax5.legend(legends)
            ax5.minorticks_on()
            ax5.set_title('Fraction_of_densities_in_past')
            ax5.set_xlabel('z')
            ax5.set_ylabel('Omega_d, Omega_m')
            ax5.grid(which='major',linewidth = 2)
            ax5.grid(which='minor')
            fig5.savefig(f'{self.output}/Omega_in_past')
            self.job.progress += 1 / 8
        #
        if config.task_7 == True:

            legends = []

            fig5 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax5 = fig5.add_subplot(111)
            for model in self.models:
                ax5.plot(model.T_i,model.OMEGA_d)
                legends.append(model.title_of_model + ', fraction_of_dark_energy')
                ax5.plot(model.T_i,model.OMEGA_m)
                legends.append(model.title_of_model + ', fraction_of_matter')
            ax5.legend(legends)
            ax5.minorticks_on()
            ax5.set_title('Fraction_of_densities_in_future')
            ax5.set_xlabel('t, Gyr')
            ax5.set_ylabel('Omega')
            ax5.grid(which='major',linewidth = 2)
            ax5.grid(which='minor')
            fig5.savefig(f'{self.output}/Omega_in_future')
            self.job.progress += 1 / 8
        #
        if config.task_8 == True:

            legends = []

            fig4 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
            ax4 = fig4.add_subplot(111)
            for model in self.models:
                ax4.plot(model.z_i,model.DA)
                legends.append(model.title_of_model)
            ax4.legend(legends)
            ax4.minorticks_on()
            ax4.set_title('Angular_diameter_distance')
            ax4.set_xlabel('z')
            ax4.set_ylabel('D_A, Gpc')
            ax4.grid(which='major',linewidth = 2)
            ax4.grid(which='minor')
            fig4.savefig(f'{self.output}/DA(z)')
            self.job.progress += 1 / 8

        fantasy_zip = zipfile.ZipFile(f'{self.output}/archive.zip', 'w')

        for folder, subfolders, files in os.walk(f'{self.output}'):

            for file in files:
                if file.endswith('.png'):
                    fantasy_zip.write(os.path.join(folder, file), os.path.relpath(os.path.join(folder,file), f'{self.output}'), compress_type = zipfile.ZIP_DEFLATED)

        fantasy_zip.close()

        return f'{self.output}/archive.zip'
