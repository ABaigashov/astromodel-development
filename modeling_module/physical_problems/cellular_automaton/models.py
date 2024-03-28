import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
from moviepy.editor import VideoFileClip
from matplotlib import colors
from PIL import Image


class Model():

	def __init__(self, config, output, job):

		self.config = config
		self.output = output
		self.job = job

		############################### GENERAP OPTIONS ###############################
		if type(self.config.frames_interval) is not int:
			self.interval = 300
		elif self.config.frames_interval < 0 or self.config.frames_interval > 1000:
			self.interval = 300
		else:
			self.interval = self.config.frames_interval
		
		if type(self.config.number_of_frames) is not int:
			self.frames = 100
		elif self.config.number_of_frames < 0 or self.config.number_of_frames > 500:
			self.frames = 100
		else:
			self.frames = self.config.number_of_frames

		if type(self.config.number_of_cells) is not int:
			self.number_of_cells = 10
		elif self.config.number_of_cells < 0 or self.config.number_of_cells > 100:
			self.number_of_cells = 10
		else:
			self.number_of_cells = self.config.number_of_cells

		############################### OBJECTS OPTIONS ###############################
		if self.config.model == 0:
			self.FULL = 255
		else:
			self.FULL = 1
			self.FIRE = 2
		self.EMPTY = 0
		
		self.grid = np.zeros((self.number_of_cells, self.number_of_cells))

		for cellular in self.config.single_cellulars:
			if type(cellular.x_coord) is not int or type(cellular.y_coord) is not int:
				pass
			elif cellular.x_coord < 0 or cellular.y_coord < 0:
				pass
			else:
				self.grid[cellular.x_coord, cellular.y_coord] = self.FULL

		for square in self.config.squares:
			if type(square.square_coord_x) is not int or type(square.square_coord_y) is not int:
				pass
			elif square.square_coord_x < 0 or square.square_coord_y < 0:
				pass
			else:
				self.grid[square.square_coord_x, square.square_coord_y] = self.FULL
				self.grid[square.square_coord_x+1, square.square_coord_y] = self.FULL
				self.grid[square.square_coord_x, square.square_coord_y+1] = self.FULL
				self.grid[square.square_coord_x+1, square.square_coord_y+1] = self.FULL

		for line in self.config.lines:
			if type(line.line_coord_x) is not int or type(line.line_coord_y) is not int:
				pass
			elif line.line_coord_x < 0 or line.line_coord_y < 0:
				pass
			else:
					self.grid[line.line_coord_x, line.line_coord_y] = self.FULL
					self.grid[line.line_coord_x, line.line_coord_y] = self.FULL
					self.grid[line.line_coord_x, line.line_coord_y] = self.FULL
		
		for fire in self.config.fires:
			if type(fire.fire_coord_x) is not int or type(fire.fire_coord_y) is not int:
				self.fire_x, self.fire_y = 0, 0
				break
			elif fire.fire_coord_x < 0 or fire.fire_coord_y < 0:
				self.fire_x, self.fire_y = 0, 0
				break
			else:
				self.fire_x, self.fire_y = fire.fire_coord_x, fire.fire_coord_y
				break
		else:
			self.fire_x, self.fire_y = 0, 0


		############################### GRAPHICAL OPTIONS ###############################
		self.output_flag = self.config.output_flag

		if self.config.model == 0:
			self.fig, self.ax = plt.subplots()
			self.mat = self.ax.matshow(self.grid, cmap='Greys')
		else:
			self.fig = plt.figure()
			self.ax = self.fig.add_subplot(111)

			self.neighbourhood = ((-1,-1), (-1,0), (-1,1), (0,-1), (0, 1), (1,-1), (1,0), (1,1))

			self.colors_list = [(0.2,0,0), (0,0.5,0), (1,0,0), 'orange']
			self.cmap = colors.ListedColormap(self.colors_list)
			self.bounds = [0, 1, 2, 3]
			self.norm = colors.BoundaryNorm(self.bounds, self.cmap.N)

			self.im = self.ax.imshow(self.grid, cmap=self.cmap, norm=self.norm)
			self.im.set_data(self.grid)
        	
		
	def game_of_life(self, frame_number):
		self.job.progress = frame_number / 100

		new_grid = self.grid.copy()

		for i in range(self.number_of_cells):
			for j in range(self.number_of_cells):
				total = (self.grid[i, (j-1)%self.number_of_cells] + self.grid[i, (j+1)%self.number_of_cells] + 
						 self.grid[(i-1)%self.number_of_cells, j] + self.grid[(i+1)%self.number_of_cells, j] + 
						 self.grid[(i-1)%self.number_of_cells, (j-1)%self.number_of_cells] + self.grid[(i-1)%self.number_of_cells, (j+1)%self.number_of_cells] + 
						 self.grid[(i+1)%self.number_of_cells, (j-1)%self.number_of_cells] + self.grid[(i+1)%self.number_of_cells, (j+1)%self.number_of_cells]) / 255
				
				if self.grid[i, j] == self.FULL:
					if total < 2 or total > 3:
						new_grid[i, j] = self.EMPTY
				else:
					if total == 3:
						new_grid[i, j] = self.FULL

		self.mat.set_data(new_grid.transpose())
		self.ax.set_title(f'Шаг: {frame_number+1}')
		self.grid = new_grid
	

	def forest_fire(self, frame_number):
		self.job.progress = frame_number / 100

		self.im.set_data(self.grid)
		self.ax.set_title(f'Шаг: {frame_number+1}')

		new_grid = np.zeros((self.number_of_cells, self.number_of_cells))

		for ix in range(1, self.number_of_cells-1):
			for iy in range(1, self.number_of_cells-1):
				if self.grid[iy, ix] == self.FULL:
					new_grid[iy, ix] = self.FULL
					for dx, dy in self.neighbourhood:
						if abs(dx) == abs(dy) and np.random.random() < 0.573:
							continue
						if self.grid[iy+dy, ix+dx] == self.FIRE:
							new_grid[self.fire_x, self.fire_y] = self.EMPTY
							new_grid[iy, ix] = self.FIRE
							break
					else:
						new_grid[self.fire_x, self.fire_y] = self.FIRE
		self.grid = new_grid


	def output_animation(self):
		if self.config.model == 0:
			for i in range(1, self.number_of_cells+1):
				self.ax.text(-0.1, i, f'{i}', va='center', ha='center')
				self.ax.text(i, -0.1, f'{i}', va='center', ha='center')

			plt.xticks(np.arange(0.5, self.number_of_cells+0.5, 1), [])
			plt.yticks(np.arange(0.5, self.number_of_cells+0.5, 1), [])
			plt.xlim(0.5, self.number_of_cells+0.5)
			plt.ylim(0.5, self.number_of_cells+0.5)
			plt.grid()

			if self.number_of_cells <= 10:
				scale = 5
			elif self.number_of_cells <= 40:
				scale = 10
			elif self.number_of_cells >= 40 and self.number_of_cells <= 80:
				scale = 20
			else:
				scale = 30

			self.fig.set_figheight(scale)
			self.fig.set_figwidth(scale)

			if self.output_flag == 0:
				ani = animation.FuncAnimation(self.fig, self.game_of_life, interval=self.interval, save_count=self.frames)
				ani.save(self.output + '.gif')
			else:
				plt.savefig(self.output + '.png')
				img = [Image.open(self.output + '.png')]
				img[0].save(fp=self.output + '.gif', format='GIF', append_images=img, save_all=True, duration=200, loop=0)
		else:
			for i in range(self.number_of_cells):
				if self.number_of_cells < 50:
					alpha = 1
					lw = 0.5
				else:
					alpha = 0.5
					lw = 0.1
				plt.plot([0.5+i, 0.5+i], [-0.5, self.number_of_cells], color='white', lw=lw, alpha=alpha)
				plt.plot([-0.5, self.number_of_cells], [0.5+i, 0.5+i], color='white', lw=lw, alpha=alpha)

				self.ax.set_axis_off()

			if self.output_flag == 0:
				ani = animation.FuncAnimation(self.fig, self.forest_fire, interval=self.interval, frames=self.frames)
				ani.save(self.output + '.gif', writer='imagemagick')
			else:
				plt.savefig(self.output + '.png')
				img = [Image.open(self.output + '.png')]
				img[0].save(fp=self.output + '.gif', format='GIF', append_images=img, save_all=True, duration=200, loop=0)
		
		clip = VideoFileClip(self.output + '.gif')
		clip.write_videofile(self.output + '.mp4')

		return self.output + '.mp4'


	def render(self):
		return self.output_animation()


def _model_from_config(config, output, job):
	model = Model
	return model(config, output, job)