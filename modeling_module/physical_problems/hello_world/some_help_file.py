
# This is hand-maded file. It contains model representation
# for this problem. You DON'T NEED to create exatcly
# this file. It's YOUR CHOICE to create help files like this


# import some required libraries
from time import sleep

# some example fucntion of astronaut state
def astronaut_state(astronaut):

	# checking astronaut age
	if astronaut.ast_age < 21:
		return 'слишком молод'
	if astronaut.ast_age > 80:
		return 'слишком стар'

	# checking astronaut weight
	if astronaut.ast_weight < 60:
		return 'недобор веса'
	if astronaut.ast_weight > 120:
		return 'перебор веса'

	# checking astronaut sex
	if astronaut.ast_sex != 'man':
		return 'В космос только мужики!'

	# all is correct
	return 'допущен'


# crating some example logger class
class MyLogger:

	# some inicialization
	def __init__(self):

		# creating output buffer
		self.output = ''

	# log method (same as print)
	def log(self, *strings, sep=' ', end='\n'):

		# printing all by default 'print'
		print(*strings, sep=' ', end='\n')

		# saving all stuff to the buffer
		self.output += sep.join(map(str, strings)) + end

	# save method (save and return)
	def save(self, path):

		# opening log file
		with open(path, 'w') as logfile:

			# write all output suff
			logfile.write(self.output)

		# returning path to the file
		return path

# crating some rocket model representation
class SomeRocketModel:

	# some inicialization
	def __init__(self, config, output, job):

		# keep incomming parameters inside 'self'
		self.config = config
		self.output = output
		self.job = job

		# creating some local logger object
		self.logger = MyLogger()

		# running some funny methods
		self.build_rocket()
		self.prepare_astronauts()

	# this is some funny example method
	# to build and prepare our rocket
	def build_rocket(self):

		# printing rocket building start message
		self.logger.log('Подготовка ракеты...')

		# fake checking jet states
		for i in range(self.config.jet_count):
			self.logger.log(f'Сопло №{i} в порядке')

		# printing rocket building end message
		self.logger.log(f'Ракета \"{self.config.rocket_name}\" готова\n')

	# this is some funny example method
	# to prepare astronauts to fly
	def prepare_astronauts(self):

		# printing astronauts preparing start message
		self.logger.log('Подготовка астронафтов...')

		# loop through every astronaut
		for astronaut in self.config.astronauts:

			# getting astronaut state
			state = astronaut_state(astronaut)

			# if all correct print this, else print error
			if state == 'допущен':
				self.logger.log(f'Астронафт \"{astronaut.ast_name}\" готов')
			else:
				self.logger.log(f'Астронафт \"{astronaut.ast_name}\" не допущен. Причина: \"{state}\"')

		# printing astronauts preparing end message
		self.logger.log('Астронафты готовы!\n')

	# this example method needs to
	# launch our rocket
	def launch_rocket(self):

		# printing rocket launching start message
		self.logger.log('Обратный отсчёт пошёл...')

		# waiting for start
		for i in range(int(self.config.rocket_time), 0, -1):

			# waiting one second
			sleep(1)

			# change progress bar
			self.job.progress += 1 / (1 + int(self.config.rocket_time))

			# if time is near to start
			# and i is not last tick
			if i < 6:

				# printing last 5 seconds countdown
				self.logger.log(i, '...')

		# waiting one more second
		sleep(1)

		# printing rocket launching end message
		self.logger.log('Ракета успешно стартовала!!!')

		# if message must be posted in news
		if self.config.post_in_news:

			# fake posting messages in news
			self.logger.log('Ваш успех был опубликован в новостях!')

		# saving all output to file
		return self.logger.save(self.output + '.txt')
