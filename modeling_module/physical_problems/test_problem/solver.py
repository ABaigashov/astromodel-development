

# crating some rocket model representation
class SomeRocketModel:

	# some inicialization
	def __init__(self, config, output, job):

		# keep incomming parameters inside 'self'
		self.config = config
		self.output = output
		self.job = job


	def saver(self):
		with open(self.output, 'w') as logfile:

			# write all output suff
			log = logfile.write(self.output + '.txt')
		# ani.save(self.output + '.gif')
		return log
