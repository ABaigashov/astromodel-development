import numpy as np
from dolfin import *
from mshr import *
#import matplotlib.pyplot as plt

def output_file(mesh, key, output):
	if key == "xmls":
		file = File(f"{output}/mesh.xml")
		file << mesh
		return f"{output}/mesh.xml"
	elif key == "matplotlib":
		plot(mesh)
		filename = f"{output}.png"
		plt.savefig(filename)
		plt.close()
		return f"{output}.png"
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

		if self.config.dimensions==2:
			domains_nul = Rectangle(Point(0, 0), Point(1, 1))
			domains_new = Rectangle(Point(0, 0), Point(1, 1))
			domains = domains_nul - domains_new
		else:
			domains_nul = Box(Point(0, 0, 0), Point(1, 1, 1))
			domains_new = Box(Point(0, 0, 0), Point(1, 1, 1))
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
					circ = CSGRotation(circ, Point(rectangle.rot_point_x, circle.rot_point_y),
					        circle.angle_of_rotation*180/np.pi)
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
					        ellipse.angle_of_rotation*180/np.pi)
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

		if self.config.cones:
			for con in self.config.cones:
				cone = Cone(Point(con.cone_x0, con.cone_y0, con.cone_z0),
							Point(con.cone_x1, con.cone_y1, con.cone_z1),
							con.cone_radius, self.config.divisions)
				if con.invert:
					domains -= cone
				else:
					domains += cone

		if self.config.cylinders:
			for cyl in self.config.cylinders:
				cylinder = Cylinder(Point(cyl.cylinder_x0, cyl.cylinder_y0, cyl.cylinder_z0),
							Point(cyl.cylinder_x1, cyl.cylinder_y1, cyl.cylinder_z1),
							cyl.cyl_radius0, cyl.cyl_radius1, self.config.divisions)
				if cyl.invert:
					domains -= cylinder
				else:
					domains += cylinder

		return domains

	def create_mesh(self):

		if self.config.mesh == "domains":
			mesh = generate_mesh(self.domains_parser(), 10)
		else:
			if self.config.unitintervalmesh:
				for interval in self.config.unitintervalmesh:
					mesh = UnitIntervalMesh(inter.number_of_cells_in_an_interval)

			if self.config.unitsquaremesh:
				for square in self.config.unitsquaremesh:
					mesh = UnitSquareMesh(square.number_of_cells_in_x_direction,
						 	square.number_of_cells_in_y_direction, square.direction_of_the_diagonals)

			if self.config.rectanglemesh:
				for rectangle in self.config.rectanglemesh:
					mesh = RectangleMesh(Point(rectangle.parall_x0, rectangle.parall_y0),
							Point(rectangle.parall_x1, rectangle.parall_x1),
							rectangle.number_of_cells_in_x_direction, rectangle.number_of_cells_in_y_direction,
							rectangle.direction_of_the_diagonals)

			if self.config.unitcubemesh:
				for cube in self.config.unitcubemesh:
					mesh = UnitCubeMesh(cube.number_of_cells_in_x_direction,
							cube.number_of_cells_in_y_direction, cube.number_of_cells_in_z_direction)

			if self.config.boxmesh:
				for box in self.config.boxmesh:
					mesh = BoxMesh(Point(box.parall_x0, box.parall_y0, box.parall_z0),
							Point(box.parall_x1, box.parall_y1, box.parall_z1),
							box.number_of_cells_in_x_direction, box.number_of_cells_in_y_direction,
							box.number_of_cells_in_z_direction)

		return output_file(mesh, self.config.output_file, self.output)




'''


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

mesh = UnitSquareMesh(10, 10, "right/left")
print("Plotting a UnitSquareMesh")
plot(mesh, title="Unit square (right/left)")
vtkfile = File('mesh_5.pvd')
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
