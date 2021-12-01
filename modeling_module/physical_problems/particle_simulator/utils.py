
# This is hand-maded file. It calls in some files while
# solving the problem. You DON'T NEED to create exatcly
# this file. It's YOUR CHOICE to create help files like this

# Import LOCAL python dict with some physical constants
from physics import constants
import generators
# Import some required libraries
from sympy import sympify
import numpy as np

# Function to scale calculating area
# Argument :edge: size (in meters)
# Returns (scaled, scaling, units)
def load_label(edge):

	# some 'if' constructions....

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
		scaling = constants['au']
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

	# Return created values
	return edge / scaling, scaling, label


# Function to load fields to the astro_object
# Arguments :config: instance of 'Configuration' object
#     :astro_object: instance of 'GlobalInteraction' object
def load_fields(config, astro_object):

	# Creating start empty values
	Emx = Emy = Emz = Hmx = Hmy = Hmz = PhiG = 0

	for field in config.EM_fields:

		# Adding electro-magnetic field parameters
		Emx += sympify(field.electricity[0])
		Emy += sympify(field.electricity[1])
		Emz += sympify(field.electricity[2])
		Hmx += sympify(field.magnetic[0])
		Hmy += sympify(field.magnetic[1])
		Hmz += sympify(field.magnetic[2])

	for field in config.G_fields:

		# Adding gravitational field parameters
		PhiG += sympify(field.gravity)

	# Load by 'append_fields' procedure
	astro_object.append_fields(Emx, Emy, Emz, Hmx, Hmy, Hmz, PhiG)



# Function to load point objects to the astro_object
# Arguments :config: instance of 'Configuration' object
#     :astro_object: instance of 'GlobalInteraction' object
def load_point_objects(config, astro_object):

	for point in config.point_objects:

		# Creating empty arrays with specific dimension
		coordinates = np.ndarray(shape=(1, config.dimensions))
		velocities = np.ndarray(shape=(1, config.dimensions))

		# Filling arrays by incomming values
		for i in range(config.dimensions):
			coordinates[0, i] = point.coords[i]
			velocities[0, i] = point.speed[i]

		if point.trajectory:
			trajectory0 = point.trajectory
		else:
			trajectory0 = config.trajectory

		# Load by 'append' procedure
		astro_object.append(*coordinates,
							*velocities,
							charge=point.charge,
							delay=point.delay,
							color=point.color,
							mass=point.mass,
							radius=point.radius,
							trajectory=trajectory0,
							id=point.id
							)

	if config.random_generators:

		generator_points = generators.RandomGenerators(config)

		for point in generator_points.output_points():

			# Creating empty arrays with specific dimension
			coordinates = np.ndarray(shape=(1, config.dimensions))
			velocities = np.ndarray(shape=(1, config.dimensions))

			# Filling arrays by incomming values
			for i in range(config.dimensions):
				coordinates[0, i] = point['coords'][i]
				velocities[0, i] = point['speed'][i]

			try:
				trajectory0 = point.trajectory
			except:
				trajectory0 = config.trajectory

			# Load by 'append' procedure
			astro_object.append(*coordinates,
								*velocities,
								charge=point['charge'],
								delay=point['delay'],
								color=point['color'],
								mass=point['mass'],
								radius=point['radius'],
								trajectory=trajectory0,
								id=point['id']
								)

# Function to load walls to the astro_object
# Arguments :config: instance of 'Configuration' object
#     :astro_object: instance of 'GlobalInteraction' object
def load_walls(config, astro_object):

	try:
		for wall in config.wall_objects:

			coordinate_1 = np.ndarray(shape=(1, config.dimensions))
			coordinate_2 = np.ndarray(shape=(1, config.dimensions))

			coordinate_1 = wall.coords_1
			coordinate_2 = wall.coords_2

			astro_object.append_wall(coordinate_1, coordinate_2, id=wall.id)
	except:
		walls = "None"
