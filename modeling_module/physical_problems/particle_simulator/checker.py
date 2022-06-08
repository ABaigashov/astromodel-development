import logging
import numpy as np


logger = logging.getLogger('my_logger')
_log_format = f"%(asctime)s - [%(levelname)s] - %(name)s - (%(filename)s).%(funcName)s(%(lineno)d) - %(message)s"


class Checker:
	def __init__(self, config, astro_object, output, job):
		self.astro_object = astro_object
		self.config = config
		self.output = output
		self.job = job
		self.logger = self.get_logger('my_logger')
		self.logger.debug('Solving')
		if not isinstance(self.config.steps_number, int):
			self.logger.warning(f"Количество шагов должно быть числом ")

		self.a = '2,5'
		self.logger.info(self.config.steps_number)
		self.objects = []
		self.check_comma()

	def check_comma(self):

		self.a = self.a.replace(',', '.')
		self.logger.info(self.a)

	def check_point(self):
		for point in self.config.point_objects:

			# Creating empty arrays with specific dimension
			coordinates = np.ndarray(shape=(1, self.config.dimensions))
			velocities = np.ndarray(shape=(1, self.config.dimensions))

			# Filling arrays by incomming values
			for i in range(self.config.dimensions):
				coordinates[0, i] = point.coords[i]
				self.logger.info(self.config.steps_number)
				velocities[0, i] = point.speed[i]

			charge=point.charge,
			delay=point.delay,
			color=point.color,
			mass=point.mass,
			radius=point.radius,
			trajectory=point.trajectory,
			id=point.id,
			K=point.K,
			destroy=point.destroy

	def get_file_handler(self):
		self.file_handler = logging.FileHandler(self.output + '.txt')
		self.file_handler.setLevel(logging.DEBUG)
		self.file_handler.setFormatter(logging.Formatter(_log_format))

		return self.file_handler

	def get_stream_handler(self):
		stream_handler = logging.StreamHandler()
		stream_handler.setLevel(logging.DEBUG)
		stream_handler.setFormatter(logging.Formatter(_log_format))
		return stream_handler

	def get_logger(self, name):
		my_logger = logging.getLogger(name)
		logger.setLevel(logging.DEBUG)
		logger.addHandler(self.get_file_handler())
		logger.addHandler(self.get_stream_handler())
		return my_logger

	def error_tester(self):
		if self.config.output_graphics == 'json':
			self.logger.info( '\n \n Данное представление реузультата \n находится в разработке. \n Ожидайте :)')
		elif self.config.output_graphics == 'vispy':
			self.logger.info ('\n \n Данное представление реузультата \n находится в разработке. \n Ожидайте :)')
		else:
			self.logger.info('\n \n Вы не указали в каком виде \n хотите получить данные моделирования. \n Пожалуйста, вернитесь в редактирование конфигурационного файла \n и заполните поле "Представление результатов"')

	def render(self):
		return self.output + '.txt'


def _log_from_config(config, astro_object, output, job):

	log = Checker
	# Returning log
	return log(config, astro_object, output, job)