from interaction_creator import GlobalInteraction
from models import _model_from_config

import os, sys
sys.path.append(os.path.abspath(os.path.join('..', '..')))
sys.path.append(os.path.abspath('modeling_module'))

from configurator import Configurator, path_generator
configuration = Configurator(sys.argv[1])

astro_object = GlobalInteraction(configuration)
model = _model_from_config(path_generator, astro_object, configuration)

print(sys.argv[3], model.render(sys.argv[2]))
