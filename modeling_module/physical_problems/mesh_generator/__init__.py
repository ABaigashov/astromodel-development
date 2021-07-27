from mesh_creator import MeshCreator
import os, sys

sys.path.append(os.path.abspath(os.path.join('..', '..')))
sys.path.append(os.path.abspath('modeling_module'))

from configurator import Configurator, path_generator
configuration = Configurator(sys.argv[1])

mesh = MeshCreator(configuration)
mesh.run_creator()

print(sys.argv[3], model.render(sys.argv[2]))
