import os, json
import numpy as np
from dolfin import *
from mshr import *


# import some required libraries
from time import sleep

# some example fucntion of astronaut state
def astronaut_state(astronaut):

	# checking astronaut age
	if astronaut.age < 21:
		return 'слишком молод'
	if astronaut.age > 80:
		return 'слишком стар'

	# checking astronaut weight
	if astronaut.weight < 60:
		return 'недобор веса'
	if astronaut.weight > 120:
		return 'перебор веса'

	# checking astronaut sex
	if astronaut.sex != 'man':
		return 'В космос только мужики!'

	# all is correct
	return 'допущен'
	

class MyLogger:

	# some inicialization
	def __init__(self):

		# creating output buffer
		self.output = ''

	# log method (same as print)
	def log(self, *strings, sep=' ', end='\n'):

		# printing all by default 'print'
		print(*strings, sep=' ', end='\n')

		# saving all stuff to the buffer
		self.output += sep.join(map(str, strings)) + end

	# save method (save and return )
	def save(self, path):

		# opening log file
		with open(path, 'w') as logfile:

			# write all output suff
			logfile.write(self.output)

		# returning path to the file
		return path

class MeshCreator:

	# some inicialization
	def __init__(self, config, output, job):

		# keep incomming parameters inside 'self'
		self.config = config
		self.output = output
		self.job = job

		# creating some local logger object
		self.logger = MyLogger()

		# running some funny methods
		self.build_rocket()
		self.prepare_astronauts()

	# this is some funny example method
	# to build and prepare our rocket
	def build_rocket(self):

		# printing rocket building start message
		self.logger.log('Подготовка ракеты...')

		# fake checking jet states
		for i in range(self.jet_count):
			self.logger.log(f'Сопло №{i} в порядке')

		# printing rocket building end message
		self.logger.log(f'Ракета \"{self.config.rocket_name}\" готова\n')

	# this is some funny example method
	# to prepare astronauts to fly
	def prepare_astronauts(self):

		# printing astronauts preparing start message
		self.logger.log('Подготовка астронафтов...')

		# loop through every astronaut
		for astronaut in self.astronauts:

			# getting astronaut state
			state = astronaut_state(astronaut)

			# if all correct print this, else print error
			if state == 'допущен':
				self.logger.log(f'Астронафт \"{astronaut.name}\" готов')
			else:
				self.logger.log(f'Астронафт \"{astronaut.name}\" не допущен. Причина: \"{state}\"')

		# printing astronauts preparing end message
		self.logger.log('Астронафты готовы!\n')

	# this example method needs to
	# launch our rocket
	def launch_rocket(self):

		# printing rocket launching start message
		self.logger.log('Обратный отсчёт пошёл...')

		# waiting for start
		for i in range(int(self.config.rocket_time), 0, -1):

			# waiting one second
			sleep(1)

			# change progress bar
			self.job.progress += 1 / (1 + int(self.config.rocket_time))

			# if time is near to start
			# and i is not last tick
			if i < 6:

				# printing last 5 seconds countdown
				self.logger.log(i, '...')

		# waiting one more second
		sleep(1)

		# printing rocket launching end message
		self.logger.log('Ракета успешно стартовала!!!')

		# if message must be posted in news
		if self.config.post_in_news:

			# fake posting messages in news
			self.logger.log('Ваш успех был опубликован в новостях!')

		# saving all output to file
		return self.logger.save(self.output + '.txt')

'''
domain = Rectangle(Point(0, 0), Point(1, 1))
domain.set_subdomain(1, Circle(Point(0.5, 0.5), 0.4))
mesh = generate_mesh(domain, 10)

file = File("2_meshes/mesh_subdomain.pvd")
file << mesh


domain = Rectangle(Point(-1, -1), Point(1, 1))
domain1 = Circle(Point(0.5, 0.5), 0.4)
domain2 = Ellipse(Point(0,0), 0.5, 1, segments=32)
domain3 = CSGRotation(domain2, 1.57)
mesh = generate_mesh(domain-domain2-domain3, 10)

file = File("2_meshes/mesh_complex.pvd")
file << mesh

#3D
domain = Cylinder(Point(0,0,0), Point(0,0,10), 5, 5, 124)
domain2 = Cylinder(Point(0,0,0), Point(0,0,10), 2, 2, 124)
mesh = generate_mesh(domain-domain2, 10)
file = File("2_meshes/cylinder.pvd")
file << mesh

domain = Cone(Point(0,0,10), Point(0,0,0), 5, 32)
mesh = generate_mesh(domain, 10)
file = File("2_meshes/cone.pvd")
file << mesh

domain=Ellipsoid(Point(0,0,0), 1, 2, 3, 15)
mesh = generate_mesh(domain, 10)
file = File("2_meshes/ellipsoid.pvd")
file << mesh

# domain = Circle(Point(1, 0.1), 0.4)
# domain2 = CSGScaling(domain, 2)
# domain3 = domain2 - domain
# mesh = generate_mesh(domain3, 10)
# file = File("2_meshes/ring.pvd")
# file << mesh

# Примеры различных сеток в DOLFIN (синтаксис)
#команда plot не работает
mesh = UnitIntervalMesh(10)
print("Plotting a UnitIntervalMesh")
plot(mesh, title="Unit interval")
vtkfile = File('mesh_1.pvd')
vtkfile << mesh

mesh = UnitSquareMesh(10, 10)
print("Plotting a UnitSquareMesh")
plot(mesh, title="Unit square")
vtkfile = File('mesh_2.pvd')
vtkfile << mesh


mesh = UnitSquareMesh(10, 10, "left")
print("Plotting a UnitSquareMesh")
plot(mesh, title="Unit square (left)")
vtkfile = File('mesh_3.pvd')
vtkfile << mesh


mesh = UnitSquareMesh(10, 10, "crossed")
print("Plotting a UnitSquareMesh")
plot(mesh, title="Unit square (crossed)")
vtkfile = File('mesh_4.pvd')
vtkfile << mesh


mesh = UnitSquareMesh(10, 10, "right/left")
print("Plotting a UnitSquareMesh")
plot(mesh, title="Unit square (right/left)")
vtkfile = File('mesh_5.pvd')
vtkfile << mesh

mesh = RectangleMesh(Point(0.0, 0.0), Point(10.0, 4.0), 10, 10)
print("Plotting a RectangleMesh")
plot(mesh, title="Rectangle")
vtkfile = File('mesh_6.pvd')
vtkfile << mesh

mesh = RectangleMesh(Point(-3.0, 2.0), Point(7.0, 6.0), 10, 10, "right/left")
print("Plotting a RectangleMesh")
plot(mesh, title="Rectangle (right/left)")
vtkfile = File('mesh_7.pvd')
vtkfile << mesh

mesh = UnitCubeMesh(10, 10, 10)
print("Plotting a UnitCubeMesh")
vtkfile = File('mesh_8.pvd')
vtkfile << mesh

mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(10.0, 4.0, 2.0), 10, 10, 10)
print("Plotting a BoxMesh")
vtkfile = File('mesh_9.pvd')
vtkfile << mesh
'''
