
# This is hand-maded file. It contains model representation
# for this problem. You DON'T NEED to create exatcly
# this file. It's YOUR CHOICE to create help files like this

# Import LOCAL python file with some utils
import utils


# Import some required libraries
from moviepy.video.io.bindings import mplfig_to_npimage
from mpl_toolkits.mplot3d.axes3d import Axes3D
from moviepy.editor import VideoClip
import matplotlib.pyplot as plt
import numpy as np
import json


# Crating base model representation
class BaseModel:

	# Some inicialization
	def __init__(self, config, astro_object, output, job):

		# Keep incomming parameters inside 'self'
		self.astro_object = astro_object
		self.config = config
		self.output = output
		self.job = job

		# Loading point objects and fields by
		# functions in 'utils' help file
		utils.load_point_objects(config, astro_object)
		utils.load_fields(config, astro_object)
		utils.load_walls(config, astro_object)

		# Saving the scale of current calculation area
		self.edge, self.scaling, self.label = utils.load_label(self.config.edge)

	# some properties...

	@property
	def colors(self):
		return np.array(self.astro_object.get_colors()) / 255

	@property
	def activity(self):
		return np.array(self.astro_object.get_activity())

	@property
	def coords(self):
		return np.array(self.astro_object.get_coords()) / self.scaling

	@property
	def coords_walls(self):
		return np.array(self.astro_object.get_coords_walls()) / self.scaling

	@property
	def radius(self):
		return np.array(self.astro_object.get_radius()) / self.scaling

	@property
	def trajectory(self):
		return np.array(self.astro_object.get_trajectories())


# Expanding base model
class JsonModel(BaseModel):

	# Render method to render current problem
	# Returns: Full path to rendered file
	def render(self):

		# Set global problem settings
		data = {
			"name": str(self.config.name),
			"step": str(self.config.step),
			"steps_number": str(self.config.steps_number),
			"frames_gap": str(self.config.frames_gap),
			"fps": str(self.config.fps),
			"edge": str(self.edge)
		}

		# Loop throw every step
		for i in range((self.config.steps_number // self.config.frames_gap)):

			# Update parameters in 'astro_object'
			for j in range(self.config.frames_gap):
				self.astro_object.update_dynamic_parametrs(self.config.step, (i * self.config.frames_gap + j) * self.config.step)

			# Write calculated parameters to frame
			data[f"frame_{i}"] = [{
				"coords": list(self.coords),
				"colors": list(self.colors),
				"radius": list(self.radius),
				"trajectory": list(self.trajectory)
			}]

			# Changing value on progressbar
			self.job.progress = i / (self.config.steps_number // self.config.frames_gap)

		# Write all created data to file
		with open(self.output + '.json', 'w') as f:
			json.dump(data, f, indent=4)

		# Return path to file
		# return self.output + '.json'


# Expanding base model with 'matplotlib' render
class PlotModel(BaseModel):

	# Some inicialization
	def __init__(self, *args, **kw):

		# Basic init of parent class
		super().__init__(*args, **kw)

		# Some 'matplotlib.pyplot' stuff
		plt.gcf().set_size_inches(8, 8)
		plt.title(self.config.name)
		plt.minorticks_on()

		# Creating trajectory array
		self.all_trajectory = []

	# Method to render single frame, which calls by
	# 'moviepy.editor.VideoClip.write_videofile' object.
	# Argument :i: frame number
	def get_frame(self, i):

		# Changing value on progressbar
		self.job.progress += 1 / (self.config.steps_number // self.config.frames_gap)

		# Update paremeters in 'astro_object'
		for j in range(self.config.frames_gap):
			self.astro_object.update_dynamic_parametrs(self.config.step, (i * self.config.frames_gap + j) * self.config.step)

		# Count and draw some other things
		self.counting()

		# Returning frame 'numpy.ndarray' representation
		return mplfig_to_npimage(plt.gcf())

	# Render method to render current problem
	# Returns: Full path to rendered file
	def render(self):

		# Creating 'moviepy.editor.VideoClip' object
		# with callback 'self.get_frame' and duration parameter
		animation = VideoClip(self.get_frame, duration=(self.config.steps_number // self.config.frames_gap) / self.config.fps)

		# Writing all created frames to 'animation'
		animation.write_videofile(self.output + '.mp4', fps=self.config.fps)

		# Return path to file
		return self.output + '.mp4'

	# Basic plug
	def counting(self):
		pass


# Creating 2D version of 'PlotModel'
class Plot2DModel(PlotModel):

	# Plug replacement with 2D version
	def counting(self):
		super().counting()

		# Some 'matplotlib.pyplot' stuff
		plt.cla()

		plt.xlabel(self.label)
		plt.ylabel(self.label)

		if self.config.grid:
			plt.grid(which='major', linewidth=1)
			plt.grid(which='minor', linewidth=0.5)

		plt.axis([-self.edge, self.edge] * 2)

		# Basic angular linspace
		angle = np.linspace(0, 2*np.pi, 100)
		#
		for s in range(len(self.coords_walls[0])):
			plt.plot([float(self.coords_walls[0][s,0]),float(self.coords_walls[1][s,0])],
			[float(self.coords_walls[0][s,1]),float(self.coords_walls[1][s,1])],'-k')
		#
		for i in range(self.coords.shape[0]):

			# Drawing a circle of a point object
			plt.plot(
				self.coords[i, 0] + self.activity[i] * self.radius[i] * np.cos(angle),
				self.coords[i, 1] + self.activity[i] * self.radius[i] * np.sin(angle),
				c=self.colors[i]
			)

				# Saving current position to next iterations
			if len(self.all_trajectory) == i:
				self.all_trajectory.append(np.array([self.coords[i]]))
			else:
				self.all_trajectory[i] = np.append(self.all_trajectory[i], np.array([self.coords[i]]), axis=0)
				# Drawing the trajectory of object
				if self.trajectory[i]:
					plt.plot(*self.all_trajectory[i].T,
						'.', ms=1, c=self.colors[i]
					)


# Creating 3D version of 'PlotModel'
class Plot3DModel(PlotModel):

	# Some inicialization
	def __init__(self, *args, **kw):

		# Basic init of parent class
		super().__init__(*args, **kw)

		# Creating 3D version of 'matplotlib.axis'
		self.ax = Axes3D(plt.gcf())

	def counting(self):
		super().counting()

		# Some 'matplotlib.pyplot' stuff
		self.ax.cla()

		self.ax.set_xlabel(self.label)
		self.ax.set_ylabel(self.label)
		self.ax.set_zlabel(self.label)

		self.ax.w_xaxis.gridlines.set_lw(3.0)
		self.ax.w_yaxis.gridlines.set_lw(3.0)
		self.ax.w_zaxis.gridlines.set_lw(3.0)

		self.ax.set_xlim3d(-self.edge, self.edge)
		self.ax.set_ylim3d(-self.edge, self.edge)
		self.ax.set_zlim3d(-self.edge, self.edge)

		# Doubble angular linspace
		theta = np.linspace(0, np.pi, 30)
		alpha = np.linspace(0, 2 * np.pi, 30)

		for i in range(self.coords.shape[0]):

			# Drawing a sphere of a point object
			self.ax.plot(
				(self.coords[i, 0] + self.activity[i] * self.radius[i] * \
					np.cos(alpha).reshape(len(theta), 1) * np.sin(theta)).flatten(),
				(self.coords[i, 1] + self.activity[i] * self.radius[i] * \
					np.sin(alpha).reshape(len(theta), 1) * np.sin(theta)).flatten(),
				(self.coords[i, 2] + self.activity[i] * self.radius[i] * \
					np.cos(theta) * np.ones((len(theta), 1))).flatten(),
				c=self.colors[i]
			)

				# Saving current position to next iterations
			if len(self.all_trajectory) == i:
				self.all_trajectory.append(np.array([self.coords[i]]))
			else:
				self.all_trajectory[i] = np.append(self.all_trajectory[i], np.array([self.coords[i]]), axis=0)
				# Drawing the trajectory of object
				if self.trajectory[i]:
					plt.plot(
						*self.all_trajectory[i].T,
						'.', ms=1, c=self.colors[i]
					)


# Function to get specified model
# Arguments :config: instance of 'Configuration' object
#     :astro_object: instance of 'GlobalInteraction' object
# Returns: specified model class object
def _model_from_config(config, astro_object, output, job):

	# some 'if' constructions....
	if config.output_graphics == 'matplotlib':
		if config.dimensions == 3:
			model = Plot3DModel
		else:
			model = Plot2DModel
	else:
		model = JsonModel

	# Returning specified model
	return model(config, astro_object, output, job)
