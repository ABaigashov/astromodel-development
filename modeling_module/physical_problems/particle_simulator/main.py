
##=================< DEFAULT BLOCK >=================##
##                                                   ##
##     This block is specified for configurating     ##
##      and constructing some spcial parameters      ##
##                                                   ##
##===================================================##

# Import some needed libraries
import os, sys

# Adding 'modeling_module' directory from upper folders
# to get 'Configurator' object from 'configurator.py'
sys.path.append(os.path.abspath('modeling_module'))

# Importing 'Configurator'
from configurator import Configurator

# Activating 'Configurator' to parse incomming parameters
configuration = Configurator._decompress(sys.argv[1], sys.argv[2])



##===============< CALCULATING BLOCK >===============##
##                                                   ##
##   This block contains YOUR calculing parameters   ##
##        and you need to write all stuff here       ##
##                                                   ##
##===================================================##

# Importing 'GlobalInteraction' object from 'interaction_creator.py'
# and '_model_from_config' function from 'models.py'
from interaction_creator import GlobalInteraction
from models import _model_from_config
# [Hint]:
#   'interaction_creator' and 'models' is LOCAL FILES, not a python libraries!!!
#    You need to write it by hands (if you need this), not installing from 'pip'

# Initializing 'GlobalInteraction' object with
# current configuration parameters
astro_object = GlobalInteraction(configuration)

# Getting specific render model for
# current configuration parameters
model = _model_from_config(astro_object, configuration)

# Render the model
result_path = model.render()



##===================< END BLOCK >===================##
##                                                   ##
##      You ALWAYS need to write a 'return-key'      ##
##         and result path. It needs to get          ##
##           file format and split output            ##
##                                                   ##
##===================================================##

# Printing the 'return-key' (here it is 'sys.argv[3]')
# and path to result of your calculation
print(sys.argv[3], result_path)
