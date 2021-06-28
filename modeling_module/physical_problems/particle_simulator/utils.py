
from sympy import symbols, sympify
from physics import constants
import numpy as np
import random


def load_label(edge):

	if edge < 1:
		scaling = 1e-2
		label = 'cantimeters'

	if edge >= 1 and edge < 1e3:
		scaling = 1
		label = 'meters'

	if edge >= 1e3 and edge < 1e6:
		scaling = 1e3
		label = 'km'

	if edge >= 1e6 and edge < 1e9:
		scaling = 1e6
		label = 'thousands km'

	if edge >= 1e9 and edge < 1e11:
		scaling = 1e9
		label = 'millions km'

	if edge >= 1e11 and edge < 3e15:
		scaling = constants['ae']
		label = 'au'

	if edge >= 3e15 and edge < 3e18:
		scaling = constants['pc']
		label = 'pc'

	if edge >= 3e18 and edge < 3e21:
		scaling = 1e3 * constants['pc']
		label = 'kpc'

	if edge >= 3e21:
		scaling = 1e6 * constants['pc']
		label = 'Mpc'

	return edge / scaling, scaling, label


def load_fields(config, astro_object):

	Emx = Emy = Emz = Hmx = Hmy = Hmz = PhiG = 0

	for field in config.EM_fields:

		Emx += sympify(field.electricity[0])
		Emy += sympify(field.electricity[1])
		Emz += sympify(field.electricity[2])
		Hmx += sympify(field.magnetic[0])
		Hmy += sympify(field.magnetic[1])
		Hmz += sympify(field.magnetic[2])

	for field in config.G_fields:
		PhiG += sympify(field.gravity)

	astro_object.append_fields(Emx, Emy, Emz, Hmx, Hmy, Hmz, PhiG)
	return Emx, Emy, Emz, Hmx, Hmy, Hmz, PhiG


def load_point_objects(config, astro_object):

	for point in config.point_objects:

		coordinates = np.ndarray(shape=(1, config.dimensions))
		velocities = np.ndarray(shape=(1, config.dimensions))

		for i in range(config.dimensions):
			coordinates[0, i] = point.coords[i]
			velocities[0, i] = point.speed[i]

		astro_object.append(
			*coordinates, *velocities, charge=point.charge, delay=point.delay,
			color=point.color, mass=point.mass, radius=point.radius, id=point.id
		)