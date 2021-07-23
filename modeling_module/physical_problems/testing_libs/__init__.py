
from time import sleep

RESULT = ''
def log(text='', end='\n'):
	global RESULT
	RESULT += str(text) + end
	print(text, end=end)


class Model:

	def init(self, config, output, job):
		self.config = config
		self.output = output
		self.job = job

		self.value = 1 + int(config.fenics) + \
			int(config.dolfin) + \
			int(config.mshr) + \
			int(config.leopart)

	def run(self):
		leopart = mshr = dolfin = fenics = False
		log('Testing started\n')

		if self.config.fenics:
			try:
				log(end='Testing "fenics"... ')
				sleep(5)
				self.job.progress += 1 / self.value
				import fenics
				log('Done')
				fenics = True
			except ImportError:
				log('Fail')
		else:
			fenics = True

		if self.config.dolfin:
			try:
				log(end='Testing "dolfin"... ')
				sleep(5)
				self.job.progress += 1 / self.value
				import dolfin
				log('Done')
				dolfin = True
			except ImportError:
				log('Fail')
		else:
			dolfin = True

		if self.config.mshr:
			try:
				log(end='Testing "mshr"... ')
				sleep(5)
				self.job.progress += 1 / self.value
				import mshr
				log('Done')
				mshr = True
			except ImportError:
				log('Fail')
		else:
			mshr = True

		if self.config.leopart:
			try:
				log(end='Testing "leopart"... ')
				sleep(5)
				self.job.progress += 1 / self.value
				import leopart
				log('Done')
				leopart = True
			except ImportError:
				log('Fail')
		else:
			leopart = True

		if all([leopart, mshr, dolfin, fenics]):
			log('\nTest passed')
		else:
			log('\nTest failed')

		output = self.output + '.log'
		with open(output, 'wb') as f:
			f.write(RESULT.encode())

		return output
