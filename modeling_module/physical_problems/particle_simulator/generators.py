# Discription ...

import random
import numpy as np

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
                            "radius": generator.radius_scale * random.random()
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
                            "radius": generator.radius_scale
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
                point_parametrs[j+counter]["coords"] = [random.uniform(-1, 1) * generator.coordinate_scale,
                                                        random.uniform(-1, 1) * generator.coordinate_scale,
                                                        0]
                point_parametrs[j+counter]["speed"] = [random.random() * generator.velocity_scale,
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
        for generator in self.config.ring_generator:
            print(generator)
            coordinate, velocity  = self.rings_generator(generator)
            for j in range(generator.particals_number):
                point_parametrs[j+counter]["coords"] = coordinate[j]
                point_parametrs[j+counter]["speed"] = velocity[j]
            counter += j + 1

        return point_parametrs

    def rings_generator(self, **kw):

        """ Функция распределения частиц в кольце объекта, с параметрами:
            xc - координата "х" центра кольца
            yc - координата "у" центра кольца
            vxc - "х" компонента скорости центрального объекта
            vyc - "у" компонента скорости центрального объекта
            r - радиус кольца, относительно центра
            N - количество точек в кольце
            V - скорость движения частиц в кольце, направление скорости касательно
                к радиусу, проведенному к точке кольца из центра (хс, ус)
        """

        # Создание таблицы с параметрами координат и скоростей объектов
        coordinate = np.ndarray(shape=(kw['N'], 3))
        velocity = np.ndarray(shape=(kw['N'], 3))

        # Создание цикла для определения координат и скоростей к заданному объекту
        for i in range(0, kw['N'], 1):

            # Определение ула для объектов
            alpha = 2 * np.pi / kw['N'] * i

            # Определение координат для первого кольца
            x = kw['xc'] + kw['r'] * np.sin(alpha)
            y = kw['yc'] + kw['r'] * np.cos(alpha)
            z = 0

            # Подставление координат в таблицу
            coordinate[i, 0] = x
            coordinate[i, 1] = y
            coordinate[i, 2] = z

            # Определение проекций скоростей в разных четвертях тригонометрической окружности
            Vx = kw['vxc'] - kw['V'] * np.cos(alpha)
            Vy = kw['vyc'] + kw['V'] * np.sin(alpha)
            Vz = 0

            # Подставление проекций скоростей  в  таблицу
            velocity[i, 0] = Vx
            velocity[i, 1] = Vy
            velocity[i, 2] = Vz

        return coordinate, velocity


class ElipsesGenerators(GeneralGenerators):
    pass
