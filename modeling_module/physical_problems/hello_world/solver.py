import numpy as np
from unpacking import Parser
import json
import time
import logging.config
import logging
from logging_config import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)
logger = logging.getLogger('my_logger')

class Solver:
    def __init__(self):
        super().__init__()
        self.logger = logger
        self.logger.debug('Solving')
        try:
            self.f = open('init_file.json', )
            self.data = json.load(self.f)
        except Exception as e:
            self.logger.error(f"Error in json file: {e} ")
            raise
        self.parser = Parser(self.data, logger)
        self.parser.objects_parse()
        self.objects = self.parser.objects
        self.N = len(self.objects)
        all_time = self.parser.years * self.parser.seconds_in_year
        self.t = np.linspace(0, all_time, self.parser.frames)
        self.dt = self.t[1]
        self.G = self.parser.G
        self.k = self.parser.k
        self.wall = self.parser.wall
        self.frames = self.parser.frames
        self.seconds = time.time()

    def get_dvx_dt(self, a, b):

        ax = 0.0
        try:
            ax += (-self.G *
                   b.mass * (
                           a.x - b.x) /
                   ((a.x - b.x) ** 2 + (a.y - b.y) ** 2) ** 1.5)
            ax += (self.k *
                   a.electric_charge * b.electric_charge / a.mass * (
                           a.x - b.x) /
                   ((a.x - b.x) ** 2 + (a.y - b.y) ** 2) ** 1.5)
        except ZeroDivisionError as e:
            self.logger.error(e)
        return float(ax)

    def get_dvy_dt(self, a, b):
        ay = 0.0
        try:

            ay += (-self.G *
                   b.mass * (
                           a.y - b.y) /
                   ((a.x - b.x) ** 2 + (a.y - b.y) ** 2) ** 1.5)
            ay += (self.k *
                   a.electric_charge * b.electric_charge / a.mass * (
                           a.y - b.y) /
                   ((a.x - b.x) ** 2 + (a.y - b.y) ** 2) ** 1.5)
        except ZeroDivisionError as e:
            self.logger.error(e)
        return float(ay)

    def hit(self, a, b):
        r = (a.x - b.x) ** 2 + (a.y - b.y) ** 2
        if r <= (100 * 2) ** 2: # !!!!!! radius_earth
            a.vx = (2 * b.mass * b.vx + a.vx * (a.mass - b.mass)) / (a.mass + b.mass)
            a.vy = (2 * b.mass * b.vy + a.vy * (a.mass - b.mass)) / (a.mass + b.mass)
        else:
            a.vx = a.vx
            b.vy = b.vy
        return a.vx, b.vy

    def wall_hit(self, a):
        if a.x <= - self.wall:
            a.vx = -a.vx
            a.x = - self.wall
        if a.x >= self.wall:
            a.vx = -a.vx
            a.x = self.wall
        if a.y <= -2 * self.wall:
            a.vy = -a.vy
            a.y = - 2 * self.wall
        if a.y >= 2 * self.wall:
            a.vy = -a.vy
            a.y = 2 * self.wall
        else:
            a.vx = a.vx
            a.vy = a.vy
        return a.vx, a.vy

    def calc_object(self, object):
        for object_n in self.objects:
            if object == object_n:
                continue
            object.vx += self.dt * self.get_dvx_dt(object, object_n)
            object.vy += self.dt * self.get_dvy_dt(object, object_n)
            self.wall_hit(object)
            self.hit(object, object_n)
        self.euler(object)

    def euler(self, object):
        object.x += self.dt * object.vx
        object.y += self.dt * object.vy

    def output_coords_func(self):
        x_coords = np.zeros((len(self.t), self.N))
        y_coords = np.zeros((len(self.t), self.N))
        time_now = time.time()
        for j in range(len(self.t)):
            for i in range(len(self.objects)):
                self.calc_object(self.objects[i])
                x_coords[j, i] = self.objects[i].x
                y_coords[j, i] = self.objects[i].y

        self.logger.info(f"Time for calculating coords: {time.time() - time_now}")
        return x_coords, y_coords
