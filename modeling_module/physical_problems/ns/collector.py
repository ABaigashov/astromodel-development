""" Модуль физических объектов и полей

Типы физических объектов и полей:

    "Point" - материальная точка, объект размеры которого малы посравнению с
    расстояниями до других объектов, с которыми происходит физическое
    взаимодействи.
"""

import numpy as np
from calculus import Star

G = 6.67408 * 10**(-11)
sun_mass = 1.989*pow(10,30)
c = 299792458

class Result_maker:

    def __init__(self):
        self.star_configurations = []

    def append_conf(self, Star):
        self.star_configurations.append(Star)

    def get_mass_profile(self):
        return [p.mass_profile for p in self.star_configurations]

    def get_rho_profile(self):
        return [p.rho_profile for p in self.star_configurations]

    def get_p_profile(self):
        return [p.p_profile for p in self.star_configurations]

    def get_masses(self):
        return [p.M for p in self.star_configurations]

    def get_central_densities(self):
        return [p.rho_profile[0] for p in self.star_configurations]

    def get_radii(self):
        r0 = 0.001 * G * sun_mass / c**2
        return [p.r*r0 for p in self.star_configurations]

    def get_radial(self):
        return [p.radial for p in self.star_configurations]
