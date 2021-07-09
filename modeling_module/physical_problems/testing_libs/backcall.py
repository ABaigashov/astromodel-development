
def value(configuration):
	return int(configuration.fenics) + \
		int(configuration.dolfin) + \
		int(configuration.mshr) + \
		int(configuration.leopart)