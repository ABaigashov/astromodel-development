
from progress.bar import ChargingBar as Bar
import os, json
import numpy as np
from dolfin import *
from mshr import *

IS_SERVER = os.getcwd() == '/usr/src'


class MeshCreator:

	def __init__(self, config, *args, **kw):
		self.config = config

	def creator(self):
		pass

	def run_creator(self):
		print(self.config)
		print(self.config.domains)


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
