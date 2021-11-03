""" Модуль физических объектов и полей

Типы физических объектов и полей:

    "Point" - материальная точка, объект размеры которого малы посравнению с
    расстояниями до других объектов, с которыми происходит физическое
    взаимодействи.
"""

import numpy as np
from scipy.integrate import odeint
from sympy import exp, log
from sympy import symbols
from sympy import sympify

G = 6.67408 * 10**(-11)
sun_mass = 1.989*pow(10,30)
c = 299792458

def EOS_analytical(expression, p0):
    expr = sympify(expression)
    p = symbols('p')
    sub = [(p,p0)]

    return float(expr.subs(sub))

def EOS_analytical_2(expression, rho0):
    expr = sympify(expression)
    rho = symbols('rho')
    sub = [(rho,rho0)]

    return float(expr.subs(sub))

def EOS_approximator(density, pressure, p):
    for i in range (len(density)):
        if p<pressure[i]:
            de = density[i]/density[i-1]
            dp = pressure[i]/pressure[i-1]
            k = np.log(de)/np.log(dp)
            break
    log_rho=np.log(density[i-1])+k*np.log(p/pressure[i-1])
    return np.exp(log_rho), k

def EOS_approximator_2(density, pressure, e):
    for i in range (len(density)):
        if e<density[i]:
            de = density[i]/density[i-1]
            dp = pressure[i]/pressure[i-1]
            k = np.log(dp)/np.log(de)
            break
    log_p=np.log(pressure[i-1])+k*np.log(e/density[i-1])
    return np.exp(log_p)

def TOV_func(s, r, k, rho1, p1):
    M0, nu0, p0 = s
    dMdr = 4*np.pi*r**2*rho1*(p0/p1)**k
    dnudr = 2*(4*np.pi*r**3*p0+M0)/(r*(r-2*M0))
    dpdr = -(4*np.pi*r**3*p0+M0)*(rho1*(p0/p1)**k+p0)/(r*(r-2*M0))

    return dMdr, dnudr, dpdr

def TOV_func_2(s, r, expression):
    M0, nu0, p0 = s
    p = symbols('p')
    expr = sympify(expression)
    sub = [(p,p0)]
    rho0 = float(expr.subs(sub))
    dMdr = 4*np.pi*r**2*rho0
    dnudr = 2*(4*np.pi*r**3*p0+M0)/(r*(r-2*M0))
    dpdr = -(4*np.pi*r**3*p0+M0)*(rho0+p0)/(r*(r-2*M0))

    return dMdr, dnudr, dpdr

class EOS:

    def __init__(self,
                name_file, form,
                rho_row, p_row,
                units_density, units_pressure):
        path = 'modeling_module/physical_problems/ns/'
        if form=="table":
            self.name_file = path + "EOS/" + name_file + ".txt"
        else:
            self.name_file = 'analytical'
        self.id = name_file
        if form=="table":
            self.pressure = []
            self.density = []
        else:
            self.pressure = p_row
            self.density = rho_row

        self.units_pressure = units_pressure
        self.units_density = units_density

        if form=="table":
            self.p_row = int(p_row)
            self.rho_row = int(rho_row)

        if self.units_density == 'MeV/fm3':
            k0 = 1.60219*10**(-6)*10**(39)
            self.K1 = 10*c**2*sun_mass / (G*sun_mass/c**2)**3/k0
        if self.units_density == 'g/cm3':
            self.K1 = 0.001*sun_mass / (G*sun_mass/c**2)**3
        if self.units_density == 'erg/cm3':
            self.K1 = 10*c**2*sun_mass / (G*sun_mass/c**2)**3

        if self.units_pressure == 'MeV/fm3':
            k0 = 1.60219*10**(-6)*10**(39)
            self.K2 = 10*c**2*sun_mass / (G*sun_mass/c**2)**3/k0
        if self.units_pressure == 'dyne/cm2':
            self.K2 = 10*c**2*sun_mass / (G*sun_mass/c**2)**3
        if self.units_pressure == 'erg/cm3':
            self.K2 = 10*c**2*sun_mass / (G*sun_mass/c**2)**3

    def EOS_maker(self):
        if  self.name_file!="analytical":
            with open(self.name_file,'r') as file:
                lines = file.read().split('\n')
            file.close()

            for i in range(len(lines)-1):
                rho = float((lines[i].split())[self.rho_row-1])/self.K1
                p = float((lines[i].split())[self.p_row-1])/self.K2
                self.density.append(rho)
                self.pressure.append(p)
        else:
            P, Rho, p, p0, rho, rho0 = symbols('P Rho p p0 rho rho0')
            p0 = sympify(self.pressure)/self.K2
            rho0 = sympify(self.density)/self.K1

            sub_rho = [(Rho,self.K1*rho)]
            sub_p = [(P,self.K2*p)]

            rho0 = rho0.subs(sub_p)
            self.density = str(rho0)
            self.pressure = str(p0.subs(sub_rho))
            print(self.density)

        return self.density, self.pressure

class Star:
    """ Класс, определяющий массивный сферический объект

    """
    def __init__(self,
                 r,dr,
                 rho,
                 eos_p,eos_rho,units_of_representation,
                 M=0,p=0,nu=0,
                 rho_profile=[],p_profile=[],
                 mass_profile=[],nu_profile=[],
                 radial=[]):

        self.r = r
        self.dr = dr
        self.rho = rho
        self.nu = nu
        self.M = M
        self.eos_rho = eos_rho
        self.eos_p = eos_p
        if type(self.eos_rho)!=str:
            self.p = EOS_approximator_2(eos_rho, eos_p, self.rho)
        else:
            self.p = EOS_analytical_2(eos_p, self.rho)

        self.mass_profile = mass_profile
        self.radial = radial
        self.p_profile = p_profile
        self.rho_profile = rho_profile
        self.nu_profile = nu_profile
        self.units_of_representation = units_of_representation

        if self.units_of_representation == 'MeV/fm3':
            k0 = 1.60219*10**(-6)*10**(39)
            self.C1 = 10*c**2*sun_mass / (G*sun_mass/c**2)**3/k0
        if self.units_of_representation == 'g/cm3':
            self.C1 = 0.001*sun_mass / (G*sun_mass/c**2)**3
        if self.units_of_representation == 'erg/cm3':
            self.C1 = 10*c**2*sun_mass / (G*sun_mass/c**2)**3

    def update_values(self):

        """ Возвращает новые координаты точки, определенные
            физическим взаимодействием.
        """
        s0 = self.M, self.nu, self.p
        rad = []

        if type(self.eos_rho)!=str:
            k = EOS_approximator(self.eos_rho, self.eos_p, self.p)[1]
        rad.append(self.r)
        rad.append(self.r+self.dr)
        if type(self.eos_rho)!=str:
            solution = odeint(TOV_func,s0,rad,args=(k,self.rho,self.p))
        else:
            solution = odeint(TOV_func_2,s0,rad,args=(self.eos_rho,))

        self.M = solution[1,0]
        self.nu = solution[1,1]
        self.p = solution[1,2]
        if type(self.eos_rho)!=str:
            self.rho = EOS_approximator(self.eos_rho, self.eos_p, solution[1,2])[0]
        else:
            self.rho = EOS_analytical(self.eos_rho, self.p)
        self.r = self.r+self.dr

        return self.M, self.nu, self.p, self.rho, self.r

    def profile(self):

        r0 = 0.001 * G * sun_mass / c**2

        self.nu_profile.append(self.nu)
        self.rho_profile.append(self.rho * self.C1)
        self.p_profile.append(self.p * self.C1)
        self.mass_profile.append(self.M)
        self.radial.append(self.r * r0)

    def stop_integration(self):
        if self.rho<1e-7 and self.dr>=0.0002:
            self.dr=self.dr/10
        message=' '
        if self.rho<1e-12:
            message='end_of_star'
        return message
