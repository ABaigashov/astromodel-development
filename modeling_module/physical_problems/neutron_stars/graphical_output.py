import matplotlib.pyplot as plt
import numpy as np

from calculus import EOS, Star
from collector import Result_maker
from sympy import sympify

import os
import shutil
import zipfile

G = 6.67408 * 10**(-11)
sun_mass = 1.989*pow(10,30)
c = 299792458
#path = 'modeling_module/physical_problems/neutron_stars/'

class Model_of_representation():

	def __init__(self, output, task1, task2):
		self.output = output
		os.mkdir(self.output)
		self.task1 = task1
		self.task2 = task2
		self.fantasy_zip = zipfile.ZipFile(f'{self.output}/archive.zip', 'w')
		units_of_representation = self.task1[0][0]
		if units_of_representation == 'MeV/fm3':
		    k0 = 1.60219*10**(-6)*10**(39)
		    self.C1 = 10*c**2*sun_mass / (G*sun_mass/c**2)**3/k0
		if units_of_representation == 'g/cm3':
		    self.C1 = 0.001*sun_mass / (G*sun_mass/c**2)**3
		if units_of_representation == 'erg/cm3':
		    self.C1 = 10*c**2*sun_mass / (G*sun_mass/c**2)**3

	def work(self):
		if self.task1[1]!=["None"]:
			fig = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
			ax = fig.add_subplot(111)
			fig2 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
			ax2 = fig2.add_subplot(111)
			fig3 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
			ax3 = fig3.add_subplot(111)
			legend=[]
			for i in self.task2:
				result=EOS(name_eos=i[0], name_file=i[1], form=i[6], rho_row=i[2], p_row=i[3],
	                    	units_density=i[4], units_pressure=i[5])
				result2=EOS.EOS_maker(result)
				summa=Result_maker()

				j=0
				for s in self.task1[1]:
					r0, dr0, rho0 = 0.0001, 0.01, float(s)/self.C1
					M0 = (4/3)*np.pi*r0**3*rho0
	                #
					result3=Star(r=r0, dr=dr0, rho=rho0,
	                        	eos_p=result2[1], eos_rho=result2[0], M=M0,
	                        	units_of_representation=self.task1[0][0],
	                        	rho_profile=[], p_profile=[], mass_profile=[],
	                        	nu_profile=[],radial=[])

					for l in range(10000):
						Star.update_values(result3)
						Star.profile(result3)
						message=Star.stop_integration(result3)
						if message=='end_of_star':
							break

					Result_maker.append_conf(summa, result3)

					r1=summa.get_radial()[j]
					m1=summa.get_mass_profile()[j]
					rho1=summa.get_rho_profile()[j]
					p1=summa.get_p_profile()[j]
					ax.plot(r1,m1)
					ax2.plot(r1,rho1)
					ax3.plot(r1,p1)
					j=j+1
					legend.append(result.id + ', ' + str(s) + ' ' + self.task1[0][0])

			ax.minorticks_on()
			ax.set_title('mass profile')
			ax.legend(legend)
			ax.set_xlabel('r, km')
			ax.set_ylabel('m(r)')
			ax.grid(which='major',linewidth = 2)
			ax.grid(which='minor')
			fig.savefig(f'{self.output}/mass_profile.png')

			ax2.minorticks_on()
			ax2.set_title('density profile')
			ax2.legend(legend)
			ax2.set_xlabel('r, km')
			ax2.set_ylabel('rho(r), ' + self.task1[0][0])
			ax2.grid(which='major',linewidth = 2)
			ax2.grid(which='minor')
			fig2.savefig(f'{self.output}/density_profile')

			ax3.minorticks_on()
			ax3.set_title('pressure profile')
			ax3.legend(legend)
			ax3.set_xlabel('r, km')
			ax3.set_ylabel('p, ' + self.task1[0][0])
			ax3.grid(which='major',linewidth = 2)
			ax3.grid(which='minor')
			fig3.savefig(f'{self.output}/pressure_profile')


			if self.task1[0][1]=="None":
				for folder, subfolders, files in os.walk(f'{self.output}'):
					for file in files:
						if file.endswith('.png') or file.endswith('.txt'):
							self.fantasy_zip.write(os.path.join(folder, file), os.path.relpath(os.path.join(folder,file), f'{self.output}'), compress_type = zipfile.ZIP_DEFLATED)

				self.fantasy_zip.close()
				return f'{self.output}/archive.zip'

	def work_2(self):

		if self.task1[0][1]!="None":
			fig = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
			ax = fig.add_subplot(111)
			fig2 = plt.figure(figsize=(8,8), facecolor='pink', frameon=True)
			ax2 = fig2.add_subplot(111)
			legend=[]
			rho_min = float(self.task1[0][1])
			rho_max = float(self.task1[0][2])
			number_of_points = int(self.task1[0][3])

			for i in self.task2:
				print("Производятся вычисления для уравнения состояния")
				print(i[0])
				result=EOS(name_eos=i[0], name_file=i[1], form=i[6], rho_row=i[2], p_row=i[3],
	                    	units_density=i[4], units_pressure=i[5])
				result2=EOS.EOS_maker(result)
				summa=Result_maker()

				for k in range(number_of_points+1):
					log_rho0 = np.log(rho_min)+k*np.log(rho_max/rho_min)/number_of_points
					r0, dr0, rho0 = 0.0001, 0.001, np.exp(log_rho0)/self.C1
					M0 = (4/3)*np.pi*r0**3*rho0

					result3=Star(r=r0, dr=dr0, rho=rho0,
			                     eos_p=result2[1], eos_rho=result2[0], M=M0,
			                     units_of_representation=self.task1[0][0],
			                     rho_profile=[], p_profile=[], mass_profile=[],
			                     nu_profile=[],radial=[])

					for l in range(100000):
						Star.update_values(result3)
						Star.profile(result3)
						message=Star.stop_integration(result3)
						if message=='end_of_star':
							break

					Result_maker.append_conf(summa, result3)
				M=summa.get_masses()
				R=summa.get_radii()
				RHO=summa.get_central_densities()
				ax.plot(R,M)
				ax2.plot(RHO,M)
				legend.append(result.id)

				Name=f'{self.output}/'+i[0]+'_M_R'+'.txt'
				G = open(Name, 'w')
				G.write('density')
				G.write(' ')
				G.write('radius')
				G.write(' ')
				G.write('mass')
				G.write('\n')
				for j in range(len(R)):
					G.write(str(RHO[j]))
					G.write(' ')
					G.write(str(R[j]))
					G.write(' ')
					G.write(str(M[j]))
					G.write('\n')
				G.close()

			ax.minorticks_on()
			ax.set_title('Mass-radius relation')
			ax.set_xlabel('R, km')
			ax.set_ylabel('M')
			ax.legend(legend)
			ax.grid(which='major',linewidth = 2)
			ax.grid(which='minor')
			fig.savefig(f'{self.output}/mass_radius_relation')

			ax2.minorticks_on()
			ax2.set_title('Mass-central density relation')
			ax2.set_xlabel('rho0, ' + self.task1[0][0])
			ax2.set_ylabel('M')
			ax2.legend(legend)
			ax2.grid(which='major',linewidth = 2)
			ax2.grid(which='minor')
			fig2.savefig(f'{self.output}/mass_density_relation')

			for folder, subfolders, files in os.walk(f'{self.output}'):
				for file in files:
					if file.endswith('.png') or file.endswith('.txt'):
						self.fantasy_zip.write(os.path.join(folder, file), os.path.relpath(os.path.join(folder,file), f'{self.output}'), compress_type = zipfile.ZIP_DEFLATED)

			self.fantasy_zip.close()

			return f'{self.output}/archive.zip'
