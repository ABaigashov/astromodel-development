# Discription ...

import random
import numpy as np
from physics import constants

class GeneralGenerators:
    def __init__(self, config):
        self.config = config

    def general_parameters(self, generators):

        point_parametrs = []

        for generator in generators:
            if generator.equal_parameters:
                for i in range(generator.particals_number):
                    point_object = {
            			"id": f"{generator.name}_id_{i}",
            			"name": f"{generator.name}_name_{i}",
                        "index": i + 10**6,
                        "charge": generator.charge_scale * random.random(),
                        "delay": generator.delay_scale * random.random(),
                        "color": generator.color,
                        "mass": generator.mass_scale * random.random(),
                        "radius": generator.radius_scale * random.random(),
                        "trajectory": generator.trajectory,
                        "K": generator.K,
                        "destroy": generator.destroy
                        }
                    point_parametrs.append(point_object)
            else:
                for i in range(generator.particals_number):
                    point_object = {
            			"id": f"{generator.name}_id_{i}",
            			"name": f"{generator.name}_name_{i}",
                        "index": i + 10**6,
                        "charge": generator.charge_scale,
                        "delay": generator.delay_scale,
                        "color": generator.color,
                        "mass": generator.mass_scale,
                        "radius": generator.radius_scale,
                        "radius": generator.radius_scale,
                        "trajectory": generator.trajectory,
                        "K": generator.K,
                        "destroy": generator.destroy
                        }
                    point_parametrs.append(point_object)

        return point_parametrs


class RandomGenerators:

    def __init__(self, config):
        self.config = config
        self.generator = GeneralGenerators(self.config)

    def output_points(self):

        point_parametrs = self.generator.general_parameters(self.config.random_generators)

        counter = 0
        for generator in self.config.random_generators:
            for j in range(generator.particals_number):
                point_parametrs[j+counter]["coords"] = [
                    random.uniform(-1, 1) * generator.coordinate_scale,
                    random.uniform(-1, 1) * generator.coordinate_scale,
                    0]
                point_parametrs[j+counter]["speed"] = [
                    random.random() * generator.velocity_scale,
                    random.random() * generator.velocity_scale,
                    0]
            counter += j + 1

        return point_parametrs


class RingsGenerators(GeneralGenerators):

    def __init__(self, config):
        self.config = config
        self.generator = GeneralGenerators(self.config)

    def output_points(self):

        point_parametrs = self.generator.general_parameters(self.config.ring_generator)

        counter = 0
        for ring in self.config.ring_generator:
            coordinate, velocity  = self.distribution_func(ring)
            for j in range(ring.particals_number):
                point_parametrs[j+counter]["coords"] = coordinate[j]
                point_parametrs[j+counter]["speed"] = velocity[j]
            counter += j + 1

        return point_parametrs

    def distribution_func(self, ring):

        # Создание таблицы с параметрами координат и скоростей объектов
        coordinate = np.ndarray(shape=(ring.particals_number, 3))
        velocity = np.ndarray(shape=(ring.particals_number, 3))

        # Создание цикла для определения координат и скоростей к заданному объекту
        for i in range(0, ring.particals_number, 1):

            # Определение ула для объектов
            alpha = 2 * np.pi / ring.particals_number * i

            # Определение координат для первого кольца
            x = ring.coords_center[0] + ring.ring_radius * np.sin(alpha)
            y = ring.coords_center[1] + ring.ring_radius * np.cos(alpha)
            z = ring.coords_center[2]

            # Подставление координат в таблицу
            coordinate[i, 0] = x
            coordinate[i, 1] = y
            coordinate[i, 2] = z

            # Определение проекций скоростей в разных четвертях тригонометрической окружности
            Vx = ring.speed_center[0] - ring.point_velocity * np.cos(alpha)
            Vy = ring.speed_center[1] + ring.point_velocity * np.sin(alpha)
            Vz = ring.speed_center[2]

            # Подставление проекций скоростей  в  таблицу
            velocity[i, 0] = Vx
            velocity[i, 1] = Vy
            velocity[i, 2] = Vz

        return coordinate, velocity


class ElipsesGenerators(GeneralGenerators):
    def __init__(self, config):
        self.config = config
        self.generator = GeneralGenerators(self.config)

    def output_points(self):

        point_parametrs = self.generator.general_parameters(self.config.ellipse_generators)

        counter = 0
        for ellipse in self.config.ellipse_generators:
            coordinate, velocity  = self.distribution_func(ellipse)
            for j in range(ellipse.particals_number):
                point_parametrs[j+counter]["coords"] = coordinate[j]
                point_parametrs[j+counter]["speed"] = velocity[j]
            counter += j + 1

        return point_parametrs

    def distribution_func(self, ellipse):

        # Создание таблицы с параметрами координат и скоростей объектов
        coordinate = np.ndarray(shape=(ellipse.particals_number+1, 3))
        velocity = np.ndarray(shape=(ellipse.particals_number+1, 3))

        b = ellipse.major_axis * np.sqrt(1-ellipse.exentricity**2) #малая полуось эллипса
        p = b**2 / ellipse.major_axis # фокальный параметр орбиты

        # Создание цикла для определения координат и скоростей к заданному объекту
        for i in range(0, ellipse.particals_number+1, 1):
            # Определение ула для объектов
            alpha = 2 * np.pi / ellipse.particals_number * i
            # Определение координат для первого кольца
            x = ellipse.coords_center[0] + ellipse.focus * p * np.cos(alpha+ellipse.apside_angle)/(1+ellipse.exentricity*np.cos(alpha))
            y = ellipse.coords_center[1] + ellipse.focus * p * np.sin(alpha+ellipse.apside_angle)/(1+ellipse.exentricity*np.cos(alpha))
            z = ellipse.coords_center[2]

            # Подставление координат в таблицу
            coordinate[i, 0] = x
            coordinate[i, 1] = y
            coordinate[i, 2] = z

            rad = np.sqrt((x-ellipse.coords_center[0])**2 + (y-ellipse.coords_center[1])**2)
            V = np.sqrt(constants['G'] * ellipse.mass_center * (2/rad - 1/ellipse.major_axis))

            # Определение проекций скоростей

            # x-я и y-я компоненты скорости на эллиптической орбите при gamma=0
            Vx0 = - V * np.sin(alpha) / np.sqrt(1+2*ellipse.exentricity*np.cos(alpha)+ellipse.exentricity**2)
            Vy0 = V * (ellipse.exentricity + np.cos(alpha)) / np.sqrt(1+2*ellipse.exentricity*np.cos(alpha)+ellipse.exentricity**2)
            Vz0 = 0

            # При повороте на угол gamma часовой стрелке скорости находятся путем
            # Преобразования поворота
            Vx = ellipse.speed_center[0] + Vx0 * np.cos(ellipse.apside_angle) - Vy0 * np.sin(ellipse.apside_angle)
            Vy = ellipse.speed_center[1] + Vx0 * np.sin(ellipse.apside_angle) + Vy0 * np.cos(ellipse.apside_angle)
            Vz = ellipse.speed_center[2]

            # Подставление проекций скоростей  в  таблицу
            velocity[i, 0] = Vx
            velocity[i, 1] = Vy

        return coordinate, velocity
