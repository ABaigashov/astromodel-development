import numpy as np
from dolfin import *
from mshr import *

def output_file(mesh, key, output):
	if key == "xmls":
		file = File(f"{output}/mesh.pvd")
		file << mesh
		return f"{output}/mesh.pvd"
	elif key == "matplotlib":
		file = File(f"{output}/mesh.pvd")
		file << mesh
		return f"{output}/mesh.pvd"
	else:
		file = File(f"{output}/mesh.pvd")
		file << mesh
		return f"{output}/mesh.pvd"


class MeshCreator:

	def __init__(self, config, output, job):
		self.config = config
		self.output = output
		self.job = job

	def domains_parser(self):

		domains_nul = Rectangle(Point(0, 0), Point(1, 1))
		domains_new = Rectangle(Point(0, 0), Point(1, 1))
		domains = domains_nul - domains_new

		if self.config.dimensions == 3:
			domains_nul=Box(Point(0, 0, 0),Point(1, 1, 1))
			domains_new = Box(Point(0, 0, 0),Point(1, 1, 1))
			domains = domains_nul - domains_new

		if self.config.rectangles:
			for rectangle in self.config.rectangles:
				rect = Rectangle(Point(rectangle.rec_x0, rectangle.rec_y0),
								 Point(rectangle.rec_x1, rectangle.rec_y1))
				if rectangle.angle_of_rotation!=0:
					rect = CSGRotation(rect, Point(rectangle.rot_point_x, rectangle.rot_point_y),
					        rectangle.angle_of_rotation*180/np.pi)
				if rectangle.invert:
					domains -= rect
				else:
					domains += rect

		if self.config.circles:
			for circle in self.config.circles:
				circ = Circle(Point(circle.circ_x_centre, circle.circ_y_centre),
								 circle.circ_radius)
				if circle.angle_of_rotation!=0:
					circ = CSGRotation(circ, Point(rectangle.rot_point_x, rectangle.rot_point_y),
					        rectangle.angle_of_rotation*180/np.pi)
				if circle.invert:
					domains -= circ
				else:
					domains += circ

		if self.config.ellipses:
			for ellipse in self.config.ellipses:
				ell = Ellipse(Point(ellipse.ell_x_centre, ellipse.ell_y_centre),
							    ellipse.ell_a,
								ellipse.ell_b)
				if ellipse.angle_of_rotation!=0:
					ell = CSGRotation(ell, Point(ellipse.rot_point_x, ellipse.rot_point_y),
					        rectangle.angle_of_rotation*180/np.pi)
				if ellipse.invert:
					domains -= ell
				else:
					domains += ell

		if self.config.parallelepipeds:
			for paralls in self.config.parallelepipeds:
				box = Box(Point(paralls.parall_x0, paralls.parall_y0, paralls.parall_z0),
							  Point(paralls.parall_x1, paralls.parall_y1, paralls.parall_z1))
				if paralls.invert:
					domains -= box
				else:
					domains += box

		if self.config.ellipsoides:
			for ellipse in self.config.ellipsoides:
				ell = Ellipsoid(Point(ellipse.sphere_x_centre, ellipse.sphere_y_centre, ellipse.sphere_z_centre),
							    ellipse.sphere_a, ellipse.sphere_b, ellipse.sphere_c)
				if ellipse.invert:
					domains -= ell
				else:
					domains += ell

		# if self.config.cones:
		# 	for con in self.config.cones:
		# 		cone = Cone(Point(con.cone_x0, con.cone_y0, con.cone_z0),
		# 					Point(con.cone_x1, con.cone_y1, con.cone_z1),
		# 					con.cone_radius)
		# 		if con.invert:
		# 			domains -= cone
		# 		else:
		# 			domains += cone

		if self.config.cylinders:
			for cyl in self.config.cylinders:
				cylinder = Cylinder(Point(cyl.cylinder_x0, cyl.cylinder_y0, cyl.cylinder_z0),
							Point(cyl.cylinder_x1, cyl.cylinder_y1, cyl.cylinder_z1),
							cyl.cyl_radius0, cyl.cyl_radius1)
				if cyl.invert:
					domains -= cylinder
				else:
					domains += cylinder

		return domains

	def create_mesh(self):
		if self.config.mesh == "domains":
			mesh = generate_mesh(self.domains_parser(), 10)
			return output_file(mesh, self.config.output_file, self.output)
		else:
			print('Poka pusto!!!')


'''



domain = Rectangle(Point(-1, -1), Point(1, 1))
domain1 = Circle(Point(0.5, 0.5), 0.4)
domain2 = Ellipse(Point(0,0), 0.5, 1, segments=32)
domain3 = CSGRotation(domain2, 1.57)
mesh = generate_mesh(domain-domain2-domain3, 10)

file = File("2_meshes/mesh_complex.pvd")
file << mesh

#3D
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
