

from sympy.matrices import Matrix
from sympy import symbols, diff
from physics import constants
import numpy as np


class GlobalInteraction:

	def __init__(self, config):
		self.points = []
		self.fields = []
		self.config = config

	def force_calculation(self, point):
		result_force = np.zeros(self.config.dimensions)
		K = self.config.K

		for field in self.fields:

			if self.config.gravity_extended_interaction:
				result_force += gravity_extended_interaction(self.config, point, field)

			if self.config.electricity_extended_interaction:
				result_force += electricity_extended_interaction(self.config, point, field)

		for p in self.points:
			if point.id != p.id and (point.activity and p.activity):
				if p.mass or p.charge:
					distance = np.linalg.norm(point.coords - p.coords, ord=2)

					if self.config.gravity_point_interaction:
						result_force = result_force + gravity_point_interaction(self.config, point, p, distance)

					if self.config.electricity_point_interaction:
						result_force = result_force + electricity_point_interaction(self.config, point, p, distance)

		return result_force

	def update_dynamic_parametrs(self, dt, time):
		self.clean_acceleration()

		K = self.config.K

		for p in self.points:

			if not p.mass and not p.charge and not p.radius:
				p.activity = False

			elif p.delay <= time:
				p.activity = True
				p.force = self.force_calculation(p)

			if p.activity:
				p.prev_pos = p.coords
				p.prev_vel = p.velocity

				if p.mass:
					p.acceleration = p.force / p.mass
				else:
					p.acceleration = p.force

				p.dv1 = p.acceleration
				p.dq1 = p.velocity
				p.coords = p.prev_pos + p.dq1 * dt / 2
				p.velocity = p.prev_vel + p.dv1 * dt / 2

		for p in self.points:
			if p.activity==1:
				p.force = self.force_calculation(p)

				if p.mass != 0:
					p.acceleration = p.force / p.mass
				else:
					p.acceleration = p.force

				p.dv2 = p.acceleration
				p.dq2 = p.velocity
				p.coords = p.prev_pos + p.dq2 * dt / 2
				p.velocity = p.prev_vel + p.dv2 * dt / 2

		for p in self.points:
			if p.activity==1:
				p.force = self.force_calculation(p)

				if p.mass != 0:
					p.acceleration = p.force / p.mass
				else:
					p.acceleration = p.force

				p.dv3 = p.acceleration
				p.dq3 = p.velocity
				p.coords = p.prev_pos + dt * p.dq3
				p.velocity = p.prev_vel + dt * p.dv3

		for p in self.points:
			if p.activity==1:
				p.force = self.force_calculation(p)

				if p.mass != 0:
					p.acceleration = p.force / p.mass
				else:
					p.acceleration = p.force

				p.dv4 = p.acceleration
				p.dq4 = p.velocity

				p.coords = (p.prev_pos + dt * (p.dq1 + 2 * p.dq2 + 2
											 * p.dq3 + p.dq4) / 6)
				p.velocity = (p.prev_vel + dt * (p.dv1 + 2 * p.dv2 + 2
												  * p.dv3 + p.dv4) / 6)

		for p1 in self.points:
			for p2 in self.points:
				if p1.id < p2.id:
					if (p1.activity and p2.activity):

						distance0 = np.linalg.norm(p1.prev_pos - p2.prev_pos, ord=2)
						distance = np.linalg.norm(p1.coords - p2.coords, ord=2)

						if distance <= p1.radius + p2.radius and distance0 > p1.radius + p2.radius:

							v1 = np.linalg.norm(p1.velocity, ord=2)
							v2 = np.linalg.norm(p2.velocity, ord=2)

							V_1 = np.zeros(self.config.dimensions)
							V_2 = np.zeros(self.config.dimensions)
							R_1 = np.zeros(self.config.dimensions)
							R_2 = np.zeros(self.config.dimensions)

							if self.config.dimensions == 3:
								e_x = np.array([0.0, 0.0, 0.0])
								e_y = np.array([0.0, 0.0, 0.0])
								e_z = np.array([0.0, 0.0, 0.0])

								if v1:
									e_x[0] = p1.velocity[0] / v1
									e_x[1] = p1.velocity[1] / v1
									e_x[2] = p1.velocity[2] / v1
								else:
									e_x[0] = 1

								e_z[0] = p2.velocity[1] * e_x[2] - p2.velocity[2] * e_x[1]
								e_z[1] = p2.velocity[2] * e_x[0] - p2.velocity[0] * e_x[2]
								e_z[2] = p2.velocity[0] * e_x[1] - p2.velocity[1] * e_x[0]

								e_z_norm = np.linalg.norm(e_z, ord=2)

								if e_z_norm:
									e_z[0] = e_z[0] / e_z_norm
									e_z[1] = e_z[1] / e_z_norm
									e_z[2] = e_z[2] / e_z_norm
								else:
									if e_x[2] != 1:
										e_z[0] = 0
										e_z[1] = 0
										e_z[2] = 1
									else:
										e_z[0] = 0
										e_z[1] = 1
										e_z[2] = 0

								e_y[0] = e_z[1] * e_x[2] - e_z[2] * e_x[1]
								e_y[1] = e_z[2] * e_x[0] - e_z[0] * e_x[2]
								e_y[2] = e_z[0] * e_x[1] - e_z[1] * e_x[0]

								T = np.array([[e_x[0], e_x[1], e_x[2]],
											  [e_y[0], e_y[1], e_y[2]],
											  [e_z[0], e_z[1], e_z[2]]])

								for i in range(0,3):
									for j in range(0,3):
										V_1[i] = T[i,j] * p1.velocity[j] + V_1[i]
										V_2[i] = T[i,j] * p2.velocity[j] + V_2[i]
										R_1[i] = T[i,j] * p1.coords[j] + R_1[i]
										R_2[i] = T[i,j] * p2.coords[j] + R_2[i]

							elif self.config.dimensions == 2:
								V_1[0] = p1.velocity[0]
								V_1[1] = p1.velocity[1]
								V_2[0] = p2.velocity[0]
								V_2[1] = p2.velocity[1]
								R_1[0] = p1.coords[0]
								R_1[1] = p1.coords[1]
								R_2[0] = p2.coords[0]
								R_2[1] = p2.coords[1]

							if v1:
								theta1 = np.arccos(V_1[0] / v1)
							else:
								theta1 = 0

							if v2:
								theta2 = np.arccos(V_2[0] / v2)
							else:
								theta2 = 0

							if (R_1[1] - R_2[1]) < 0:
								phi = - np.arccos((R_1[0] - R_2[0]) / distance) + 2 * np.pi
							else:
								phi = np.arccos((R_1[0] - R_2[0]) / distance)

							if V_1[1] < 0:
								theta1 = - theta1 + 2 * np.pi

							if V_2[1] < 0:
								theta2 = - theta2 + 2 * np.pi

							if not p1.mass==0 and not p2.mass:
								V_1[0] = v1 * np.cos(theta1 - phi) * (1 - K) \
										* np.cos(phi) / 2\
										+ ((1 + K) * v2 * np.cos(theta2 - phi))\
										* np.cos(phi) / 2\
										+ K * v1 * np.sin(theta1 - phi) * np.cos(phi + np.pi / 2)

								V_1[1] = v1 * np.cos(theta1 - phi) * (1 - K) \
										* np.sin(phi) / 2 \
										+ ((1 + K) * v2 * np.cos(theta2 - phi)) \
										* np.sin(phi) / 2 \
										+ K * v1 * np.sin(theta1 - phi) * np.sin(phi + np.pi / 2)

								V_2[0] = v2 * np.cos(theta2 - phi) * (1 - K) \
										* np.cos(phi) / 2\
										+ ((1 + K) * v1 * np.cos(theta1 - phi)) \
										* np.cos(phi) / 2\
										+ K * v2 * np.sin(theta2 - phi) * np.cos(phi + np.pi / 2)

								V_2[1] = v2 * np.cos(theta2 - phi) * (1 - K) \
										* np.sin(phi) / 2 \
										+ ((1 + K) * v1 * np.cos(theta1 - phi)) \
										* np.sin(phi) / 2\
										+ K * v2 * np.sin(theta2 - phi) * np.sin(phi + np.pi / 2)
							else:
								V_1[0] = v1 * np.cos(theta1 - phi) * (p1.mass - K * p2.mass) \
										* np.cos(phi) / (p1.mass + p2.mass)\
										+ ((1 + K) * p2.mass * v2 * np.cos(theta2 - phi))\
										* np.cos(phi) / (p1.mass + p2.mass)\
										+ K * v1 * np.sin(theta1 - phi) * np.cos(phi + np.pi / 2)

								V_1[1] = v1 * np.cos(theta1 - phi) * (p1.mass - K * p2.mass) \
										* np.sin(phi) / (p1.mass + p2.mass) \
										+ ((1 + K) * p2.mass * v2 * np.cos(theta2 - phi)) \
										* np.sin(phi) / (p1.mass + p2.mass) \
										+ K * v1 * np.sin(theta1 - phi) * np.sin(phi + np.pi / 2)

								V_2[0] = v2 * np.cos(theta2 - phi) * (p2.mass - K * p1.mass) \
										* np.cos(phi) / (p1.mass + p2.mass)\
										+ ((1 + K) * p1.mass * v1 * np.cos(theta1 - phi)) \
										* np.cos(phi) / (p1.mass + p2.mass)\
										+ K * v2 * np.sin(theta2 - phi) * np.cos(phi + np.pi / 2)

								V_2[1] = v2 * np.cos(theta2 - phi) * (p2.mass - K * p1.mass) \
										* np.sin(phi) / (p1.mass + p2.mass) \
										+ ((1 + K) * p1.mass * v1 * np.cos(theta1 - phi)) \
										* np.sin(phi) / (p1.mass + p2.mass)\
										+ K * v2 * np.sin(theta2 - phi) * np.sin(phi + np.pi / 2)

							v_1 = np.zeros(self.config.dimensions)
							v_2 = np.zeros(self.config.dimensions)
							if self.config.dimensions == 3:
								for j in range(0,3):
									for i in range(0,3):
										v_1[j] = T[i,j] * V_1[i] + v_1[j]
										v_2[j] = T[i,j] * V_2[i] + v_2[j]

							elif self.config.dimensions == 2:
								v_1[0] = V_1[0]
								v_1[1] = V_1[1]
								v_2[0] = V_2[0]
								v_2[1] = V_2[1]

							for i in range(self.config.dimensions):
								p1.velocity[i] = v_1[i]
								p2.velocity[i] = v_2[i]

							if not K:
								if p1.mass >= p2.mass:
									p1.mass = p1.mass + p2.mass
									p1.charge = p1.charge + p2.charge
									p1.radius = (p1.radius**3 + p2.radius**3)**(1/3)
									p2.mass = 0
									p2.charge = 0
									p2.radius = 0
								else:
									p2.mass = p1.mass + p2.mass
									p2.charge = p1.charge + p2.charge
									p2.radius = (p1.radius**3 + p2.radius**3)**(1/3)
									p1.mass = 0
									p1.charge = 0
									p1.radius = 0

	def clean_acceleration(self):
		for p in self.points:
			p.clean_acceleration()

	def append(self, *args, **kwargs):
		self.points.append(Point(*args, **kwargs))

	def append_fields(self, *args, **kwargs):
		self.fields.append(Field(*args, **kwargs))

	def get_coords(self):
		return [p.coords for p in self.points]

	def get_activity(self):
		return [p.activity for p in self.points]

	def get_radius(self):
		return [p.radius * self.config.scale_faktor for p in self.points]

	def get_colors(self):
		return [p.color for p in self.points]


def gravity_point_interaction(config, point1, point2, r):
	gravity_force = np.zeros(config.dimensions)

	for i in range(config.dimensions):

		if point1.mass and point2.mass:
			gravity_force[i] = -constants['G'] * point1.mass * point2.mass * (point1.coords[i]-point2.coords[i]) / (r**3)

		elif not point1.mass:
			gravity_force[i] = -constants['G'] * point2.mass * (point1.coords[i]-point2.coords[i]) / (r**3)

		elif not point2.mass:
			gravity_force[i] = -constants['G'] * point1.mass * (point1.coords[i]-point2.coords[i]) / (r**3)

	return gravity_force


def gravity_extended_interaction(config, point1, field):
	components = field.get_fields(point1.coords)
	grav_force = np.zeros(config.dimensions)

	if point1.mass:
		for i in range(config.dimensions):
			grav_force[i] = point1.mass * components[2][i]
	else:
		for i in range(config.dimensions):
			grav_force[i] = components[2][i]

	return grav_force


def electricity_point_interaction(config, point1, point2, r):
	coulomb_force = np.zeros(config.dimensions)

	for i in range(config.dimensions):
		if point1.mass and point2.mass:
			coulomb_force[i] = constants['k'] * point1.charge * point2.charge * (point1.coords[i] - point2.coords[i]) / (r**3)

	return coulomb_force


def electricity_extended_interaction(config, point1, field):
	lorentz_force = np.zeros(config.dimensions)
	components = field.get_fields(point1.coords)

	if point1.mass and config.dimensions == 3:
		lorentz_force[0] = point1.velocity[1] * components[1][2] * point1.charge\
						 - point1.velocity[2] * components[1][1] * point1.charge\
						 + point1.charge * components[0][0]

		lorentz_force[1] = point1.velocity[2] * components[1][0] * point1.charge\
						 - point1.velocity[0] * components[1][2] * point1.charge\
						 + point1.charge * components[0][1]

		lorentz_force[2] = point1.velocity[0] * components[1][1] * point1.charge\
						 - point1.velocity[1] * components[1][0] * point1.charge\
						 + point1.charge * components[0][2]

	if point1.mass and config.dimensions == 2:

		lorentz_force[0] = point1.velocity[1] * components[1][2] * point1.charge\
						 + point1.charge * components[0][0]

		lorentz_force[1] = - point1.velocity[0] * components[1][2] * point1.charge\
						 + point1.charge * components[0][1]

	return lorentz_force



class Field:

	def __init__(self, Em_x=0, Em_y=0,
				 Em_z=0, Hm_x=0, Hm_y=0,
				 Hm_z=0, phi=0):

		self.Em = Matrix([[0,0,0]])
		self.Hm = Matrix([[0,0,0]])
		self.grav = Matrix([[0,0,0]])

		self.Em[0] = Em_x
		self.Em[1] = Em_y
		self.Em[2] = Em_z
		self.Hm[0] = Hm_x
		self.Hm[1] = Hm_y
		self.Hm[2] = Hm_z

		x, y, z = symbols('x y z')
		self.grav[0] = -diff(phi, x)
		self.grav[1] = -diff(phi, y)
		self.grav[2] = -diff(phi, z)

	def get_fields(self, coords):

		x, y, z = symbols('x y z')
		sub = [(x, coords[0]), (y, coords[1])]
		electric_field = self.Em.subs(sub)
		magnetic_field = self.Hm.subs(sub)
		grav_field = self.grav.subs(sub)

		return electric_field, magnetic_field, grav_field

class Point:

	def __init__(self, coords, velocity, id=0,
				 activity=1, delay=0, mass=1.0,
				 charge=1.0, radius=0, color=None):

		self.id = id
		self.coords = coords
		self.velocity = velocity
		self.acceleration = np.zeros(len(coords))
		self.mass = mass
		self.charge = charge
		self.radius = radius
		self.color = color
		self.delay = delay

		if delay == 0:
			self.activity = 1
		else:
			self.activity = 0

		self.force = self.dv1 = self.dv2 = self.dv3 = \
		self.dv4 = self.dq1 = self.dq2 = self.dq3 = \
		self.dq4 = self.prev_pos = self.prev_vel = None

	def update_position(self, dt):
		pass

	def clean_acceleration(self):
		self.acceleration = self.acceleration * 0
		self.dv1 = self.dv2 = self.dv3 = self.dv4 = \
		self.dq1 = self.dq2 = self.dq3 = self.dq4 = 0
