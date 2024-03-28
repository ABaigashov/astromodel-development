
# This is hand-maded file. It contains model representation
# for this problem. You DON'T NEED to create exatcly
# this file. It's YOUR CHOICE to create help files like this
# Import some required libraries
import numpy as np
import json
import time


# Expanding base model with 'matplotlib' render
class Model():

	# Some inicialization
	def __init__(self, config, output, job):

		# Keep incomming parameters inside 'self'
		self.config = config
		self.output = output
		self.job = job


	def log(self, path):
		for i in np.arange(0, 1, 0.1):
			self.job.progress = i
			time.sleep(1)

		output = 'Здравствуйте, я бацила'
		with open(path, 'w') as logfile:

			# write all output suff
			logfile.write(output)

		# returning path to the file
		return path


	def render(self):
		return self.log(self.output + '.txt')


def _model_from_config(config, output, job):
	model = Model
	return model(config, output, job)