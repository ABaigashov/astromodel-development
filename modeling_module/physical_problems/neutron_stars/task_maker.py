#""" Модуль для парсинга данных, полученных из json файла

import numpy as np


# ------------------ General scene parameters ------------------
class Task_maker:

    def __init__(self, config):
        self.config = config

		# keep incomming parameters inside 'self'

        self.EoS_descriptor = []
        self.task_descriptor = []

    def EoS_append(self):

        if self.config.EoS_parameters:
            EoS_parameters_1 = self.config.EoS_parameters
            for EoS_parameters in EoS_parameters_1:
                if EoS_parameters.eos_name:
                    name = EoS_parameters.eos_name
                if EoS_parameters.eos_file:
                    eos_file = EoS_parameters.eos_file
                if EoS_parameters.analytical:
                    if EoS_parameters.analytical=="yes":
                        form = "analytical"
                    else:
                        form = "table"
                else:
                    form = "table"
                if EoS_parameters.density_row:
                    density_row = EoS_parameters.density_row
                else:
                    density_row = 1
                if EoS_parameters.pressure_row:
                    pressure_row = EoS_parameters.pressure_row
                else:
                    pressure_row = 2
                if EoS_parameters.units_density:
                    units_density = EoS_parameters.units_density
                else:
                    units_density = "g/cm3"
                if EoS_parameters.units_pressure:
                    units_pressure = EoS_parameters.units_pressure
                else:
                    units_pressure = "dyne/cm2"

                EoS_descriptor = [name, eos_file, density_row, pressure_row, units_density,
                                 units_pressure, form]
                self.EoS_descriptor.append(EoS_descriptor)

    def Task_append(self):
        if self.config.units_of_representation:
            units_of_representation = self.config.units_of_representation
        else:
            units_of_representation = "MeV/fm3"
        if self.config.rho_min:
            rho_min = self.config.rho_min
        else:
            rho_min = "None"
        if self.config.rho_max:
            rho_max = self.config.rho_max
        else:
            rho_max = "None"
        if self.config.number_of_points:
            number_of_points = self.config.number_of_points
        else:
            number_of_points = 50
        Task_descriptor = [units_of_representation, rho_min, rho_max, number_of_points]
        self.task_descriptor.append(Task_descriptor)
        rho_i = "None"
        if self.config.rho_i:
            rho_i = []
            # values = self.config.rho_i
            # for i in values:
            #     rho_i.append(i)
            rho_i.append(self.config.rho_i)
        self.task_descriptor.append(rho_i)
