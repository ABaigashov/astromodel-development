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
        print(self.config.marker, self.config.ball_color, self.config.ball_color)
        ball, = plt.plot([], [], self.config.marker, color=self.config.ball_color, label=self.config.label, lw=2)

        def circle_move(R, vx0, vy0, time):
            x0 = vx0 * time
            y0 = vy0 * time
            alpha = np.arange(0, 2*np.pi, 0.1)
            x = x0 + R*np.cos(alpha)
            y = y0 + R*np.sin(alpha)
            return x, y

        edge = self.config.edge
        plt.axis('equal')
        ax.set_xlim(-edge, edge)
        ax.set_ylim(-edge, edge)

        def animate(i):
            ball.set_data(circle_move(R=self.config.radius, vx0=self.config.vx0, vy0=self.config.vy0, time=i))
        
        ani = FuncAnimation(fig,
                            animate,
                            frames=100,
                            interval=30
                        )
                        
        return ani.save(f"{self.output}.gif")
        