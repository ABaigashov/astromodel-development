import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation


class BallModel:
    def __init__(self, config, output, job):
        self.config = config
        self.output = output
        self.job = job

    def start(self):
        fig, ax = plt.subplots()

        for b in self.config.balls:
            self.ball, = plt.plot([], [], self.config.b.marker, color=self.config.b.ball_color,
                                  label=self.config.b.label, lw=2)

        edge = self.config.edge
        plt.axis('equal')
        ax.set_xlim(-edge, edge)
        ax.set_ylim(-edge, edge)

        ani = FuncAnimation(fig,
                            self.animate,
                            frames=100,
                            interval=30
                            )

        return ani.save(f"{self.output}.gif")

    def circle_move(self, R, vx0, vy0, time):
        x0 = vx0 * time
        y0 = vy0 * time
        alpha = np.arange(0, 2 * np.pi, 0.1)
        x = x0 + R * np.cos(alpha)
        y = y0 + R * np.sin(alpha)
        return x, y

    def animate(self, i):
        self.ball.set_data(self.circle_move(R=self.config.b.radius, vx0=self.config.b.vx0, vy0=self.config.b.vy0, time=i))
