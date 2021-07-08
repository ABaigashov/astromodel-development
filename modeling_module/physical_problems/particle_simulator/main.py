from interaction_creator import GlobalInteraction
from models import _model_from_config

import os, sys
sys.path.append(os.path.abspath(os.path.join('..', '..')))
sys.path.append(os.path.abspath('modeling_module'))

from configurator import Configurator
configuration = Configurator(sys.argv[1], sys.argv[2])

astro_object = GlobalInteraction(configuration)
model = _model_from_config(astro_object, configuration)


print(sys.argv[3], model.render())
