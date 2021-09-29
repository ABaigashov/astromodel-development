import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.tri as mtri
import sympy as sym
from sympy import symbols
from sympy import sympify


# from mpl_toolkits.mplot3d import *
# from matplotlib import cm
# import matplotlib.animation as animation
############################
#####УРАВНЕНИЕ ПУАССОНА#####
############################
path = 'modeling_module/physical_problems/cosmology/'

class Task_maker():

    def __init__(self, config):
        self.config = config

        self.omega_m = self.config.omega_m
        self.omega_r = self.config.omega_r
        self.z_max = self.config.z_max
        self.t_max = self.config.t_max
        self.Hubble = self.config.H_0
        self.dark_components = []

        if self.config.dark_components:
            for compts in self.config.dark_components:
                if compts.equation_d:
                    self.dark_components.append(compts.equation_d)
                if compts.omega_d:
                    self.omega_d.append(compts.omega_d)
