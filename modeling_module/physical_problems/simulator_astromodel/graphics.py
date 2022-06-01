from solver import Solver
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from logger import Logger


class Graphics(Solver, Logger):
    def __init__(self):
        super().__init__()
        self.logger.debug('Start drawing')
        output_coords = self.output_coords_func()
        self.x_point = output_coords[0]
        self.y_point = output_coords[1]
        self.edge = 2*3*10**11
        self.x_line = np.zeros((self.N, len(self.t)))
        self.y_line = np.zeros((self.N, len(self.t)))
        self.points = []
        self.lines = []
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlim(-self.edge, self.edge)
        self.ax.set_ylim(-self.edge, self.edge)
        self.creator()
    def creator(self):
        for j in range(self.N):
            for i in range(len(self.t)):
                self.x_line[j][i] = self.x_point[i][j]
                self.y_line[j][i] = self.y_point[i][j]

        for j in range(self.N):
            point, = plt.plot([], [], 'o', ms=5)
            self.points.append(point)
            line, = plt.plot([], [], '-')
            self.lines.append(line)

    def animate(self, i):
        for j in range(self.N):
            self.points[j].set_data(self.x_point[i][j], self.y_point[i][j])
            self.lines[j].set_data(self.x_line[j][:i], self.y_line[j][:i])
        return self.points, self.lines

    def run(self):
        plt.axis('equal')
        self.logger.info("DONE")
        self.logger.info(f"                                 \n"
                          "  ____  _   _ ____  _____ ____   \n"
                          " / ___|| | | |  _ \| ____|  _ \  \n"
                          " \___ \| | | | |_) |  _| | |_) | \n"
                          "  ___) | |_| |  __/| |___|  _ <  \n"
                          " |____/ \___/|_|   |_____|_| \_\  \n ")
        ani = FuncAnimation(self.fig, self.animate, frames=self.frames, interval=50)
		ani.save(self.output + '.gif')

		# Return path to file
		return self.output + '.gif'
