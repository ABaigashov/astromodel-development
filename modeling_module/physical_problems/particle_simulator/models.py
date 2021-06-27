
from moviepy.video.io.bindings import mplfig_to_npimage
from mpl_toolkits.mplot3d.axes3d import Axes3D
from progress.bar import ChargingBar as Bar
from moviepy.editor import VideoClip
import matplotlib.pyplot as plt
import os, json, utils
import numpy as np

IS_SERVER = os.getcwd() == '/usr/src'


class SlowBar(Bar):
	suffix = ', remaining: %(remaining_hours)s'

	@property
	def remaining_hours(self):
		if self.eta > 3600:
			time = round(self.eta / 3600, 2)
			ut = 'hours'
		elif (self.eta < 3600) and (self.eta > 60):
			time = round(self.eta / 60, 2)
			ut = 'minutes'
		else:
			time = round(self.eta, 3)
			ut = 'seconds'
		return '{} '.format(time) + '{}'.format(ut)


class BaseModel:

	def __init__(self, astro_object, config, *args, **kw):
		self.astro_object = astro_object
		self.config = config

		utils.load_point_objects(config, astro_object)
		utils.load_fields(config, astro_object)

		self.edge, self.scaling, self.label = utils.load_label(self.config.edge)
		self.bar = SlowBar('Processing', max=self.config.steps_number + 1)

	@property
	def colors(self):
		return np.array(self.astro_object.get_colors()) / 255

	@property
	def activity(self):
		return np.array(self.astro_object.get_activity())

	@property
	def radius(self):
		return np.array(self.astro_object.get_radius())

	@property
	def coords(self):
		return np.array(self.astro_object.get_coords()) / self.scaling

	@property
	def radius_scale(self):
		return np.array(self.astro_object.get_radius()) / self.scaling


class JsonModel(BaseModel):

	def __init__(self, *args, **kw):
		super().__init__(*args, **kw)
		self.render = kw['path_generator']('json', self.render)

	def render(self, path):
		data = {
			"name": str(self.config.name),
			"step": str(self.config.step),
			"steps_number": str(self.config.steps_number),
			"fps": str(self.config.fps),
			"edge": str(self.edge),
			"frames_gap": str(self.config.frames_gap),
			"trajectory": str(self.config.trajectory),
		}
		for i in range(self.config.steps_number):
			self.astro_object.update_dynamic_parametrs(self.config.step, i * self.config.frames_gap * self.config.step)
			data[f"frame_{i}"] = [{
				"coords": list(self.coords),
				"colors": list(self.colors),
				"radius": list(self.radius)
			}]
		with open(path, 'w') as f:
			json.dump(data, f, indent=4)


class PlotModel(BaseModel):

	def __init__(self, *args, **kw):
		super().__init__(*args, **kw)
		plt.gcf().set_size_inches(8, 8)
		plt.title(self.config.name)
		plt.minorticks_on()
		self.all_trajectory = []
		self.render = kw['path_generator']('mp4', self.render)

	def get_frame(self, i):
		self.astro_object.update_dynamic_parametrs(self.config.step, i * self.config.frames_gap * self.config.step)
		self.counting()
		if not IS_SERVER:
			self.bar.next()
		return mplfig_to_npimage(plt.gcf())

	def render(self, path):
		animation = VideoClip(self.get_frame, duration=(self.config.steps_number // self.config.frames_gap) / self.config.fps)
		animation.write_videofile(path, fps=self.config.fps)

	def counting(self):
		pass


class Plot2DModel(PlotModel):

	def counting(self):
		plt.cla()

		plt.xlabel(self.label)
		plt.ylabel(self.label)

		plt.grid(which='major', linewidth=1)
		plt.grid(which='minor', linewidth=0.5)

		plt.axis([-self.edge, self.edge] * 2)
		angle = np.linspace(0, 2*np.pi, 100)

		for i in range(self.coords.shape[0]):

			plt.plot(
				self.coords[i, 0] + self.activity[i] * self.radius_scale[i] * np.cos(angle),
				self.coords[i, 1] + self.activity[i] * self.radius_scale[i] * np.sin(angle),
				c=self.colors[i]
			)

			if self.config.trajectory:

				if len(self.all_trajectory) == i:
					self.all_trajectory.append(np.array([self.coords[i]]))
				else:
					self.all_trajectory[i] = np.append(self.all_trajectory[i], np.array([self.coords[i]]), axis=0)

				plt.plot(
					*self.all_trajectory[i].T,
					'.', ms=1, c=self.colors[i]
				)

class Plot3DModel(PlotModel):

	def __init__(self, *args, **kw):
		super().__init__(*args, **kw)
		self.ax = Axes3D(plt.gcf())

	def counting(self):
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

		theta = np.linspace(0, np.pi, 30)
		alpha = np.linspace(0, 2 * np.pi, 30)

		for i in range(self.coords.shape[0]):

			self.ax.plot(
				(self.coords[i, 0] + self.activity[i] * self.radius_scale[i] * \
					np.cos(alpha).reshape(len(theta), 1) * np.sin(theta)).flatten(),
				(self.coords[i, 1] + self.activity[i] * self.radius_scale[i] * \
					np.sin(alpha).reshape(len(theta), 1) * np.sin(theta)).flatten(),
				(self.coords[i, 2] + self.activity[i] * self.radius_scale[i] * \
					np.cos(theta) * np.ones((len(theta), 1))).flatten(),
				c=self.colors[i]
			)

			if self.config.trajectory:

				if len(self.all_trajectory) == i:
					self.all_trajectory.append(np.array([self.coords[i]]))
				else:
					self.all_trajectory[i] = np.append(self.all_trajectory[i], np.array([self.coords[i]]), axis=0)

				self.ax.plot(
					*self.all_trajectory[i].T,
					'.', ms=1, c=self.colors[i]
				)

def _model_from_config(path_generator, astro_object, config):
	if config.output_graphics == 'matplotlib':
		if config.dimensions == 3:
			model = Plot3DModel
		else:
			model = Plot2DModel
	else:
		model = JsonModel
	return model(astro_object, config, path_generator=path_generator)
