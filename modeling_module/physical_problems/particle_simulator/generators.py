# Discription ...

import random


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
    pass


class ElipsesGenerators(GeneralGenerators):
    pass
