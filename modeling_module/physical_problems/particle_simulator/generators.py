# Discription ...

import random

class GeneratorsBase:
    def __init__(self, config):
        self.config = config

    def output_points(self):

        point_parametrs = []

        # for i in range(self.config.particals_number):
        #     point_object = {
        # 			"id": f"id_{i}",
        # 			"name": f"name_{i}",
        #             "index": i + 10**5,
        #             "charge": self.config.charge * random.random(),
        #             "delay": self.config.delay * random.random(),
        #             "coords": [random.uniform(-1, 1) * self.config.COORDINATE_SIZE,
        #                        random.uniform(-1, 1) * self.config.COORDINATE_SIZE,
        #                        0],
        #             "color":[int(random.random()*255),
        #                      int(random.random()*255),
        #                      int(random.random()*255)],
        #             "speed": [random.random(), random.random(), 0],
        #             "mass": self.config.mass * random.random(),
        #             "radius": self.config.radius * random.random()
        #             }
        #     point_parametrs.append(point_object)
        #
        # return point_parametrs


        ae = 149e+9
        PARTICALS_NUMBER = 20
        COORDINATE_SIZE = ae
        VELOCITY_SIZE = 30000
        MASS_SIZE = 1e+25
        RADIUS_SIZE = 0.05 * ae

        for i in range(PARTICALS_NUMBER):
            point_object = {
        			"id": f"id_{i}",
        			"name": f"name_{i}",
                    "index": i + 10**5,
                    "charge": 0,
                    "delay": 0,
                    "coords": [random.uniform(-1, 1) * COORDINATE_SIZE,
                               random.uniform(-1, 1) * COORDINATE_SIZE,
                               0],
                    "color":[int(random.random()*255), int(random.random()*255), int(random.random()*255)],
                    "speed": [random.random() * VELOCITY_SIZE, random.random() * VELOCITY_SIZE, 0],
                    "mass": MASS_SIZE * random.random(),
                    "radius": RADIUS_SIZE * random.random()
                    }
            point_parametrs.append(point_object)

        return point_parametrs
