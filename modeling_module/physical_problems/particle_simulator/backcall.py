
##==================< INFORMATION >==================##
##                                                   ##
##     This file has special name 'backcall' and     ##
##   shows how many backcalls did your code prints.  ##
##   The 'backcall' is the step on a progress line   ##
##   witch calls like'print(end="\b")'. If problem   ##
##   doesn't need any progress line, you can IGNORE  ##
##  this file and DON'T ADD it to problem directory  ##
##                                                   ##
##===================================================##

# Define special function
# Argument :configuration: is 'Configurator' class object
# Returns :int: count of EXACT amount 'backcall' callings
def value(configuration):
	# In our case 'moviepy.editor.VideoClip' object
	# calls 'render' function 'steps_number / frames_gap' times
	return (configuration.steps_number // configuration.frames_gap)