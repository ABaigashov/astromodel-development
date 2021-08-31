

import serial


class Telescope:

	def __init__(self, port, baudrate=9600, timeout=1):
		self.port = port
		self.baudrate = baudrate
		self.timeout = timeout
		self.serial = None

	def __enter__(self):
		self.serial = serial.Serial(self.port, baudrate=self.baudrate, timeout=self.timeout)
		return self

	def __exit__(self, exc_type, exc_value, trace):
		self.serial.close()
		self.serial = None
		return False

	def command(self, command):
		command = command.encode() + b'\r'
		if self.serial is None:
			with self as telescope:
				return self._command(telescope.serial, command)
		return self._command(self.serial, command)

	@staticmethod
	def _command(serial, command):
		serial.write(command)
		return serial.readline()

	@staticmethod
	def _decode_coordinate(coordinate):
		if len(coordinate) == 9:
			return tuple(int(coord, base=16) / 65536 for coord in coordinate.split(','))
		if len(coordinate) == 17:
			return tuple(int(coord[:-2], base=16) / 16777216 for coord in coordinate.split(','))

	@classmethod
	def _ra_dec(cls, coordinate):
		_ra, _dec = coordinate
		ra, dec = cls._ra(_ra), cls._dec(_dec)
		return f'RA( {ra} ); DEC( {dec} )'

	@staticmethod
	def _ra(ra):
		_deg = ra * 360
		deg = int(_deg)
		_minutes = 60 * (_deg - deg)
		minutes = int(_minutes)
		_seconds = 60 * (_minutes - minutes)
		seconds = int(_seconds)
		return f'{deg}° {minutes}’ {seconds}”'

	@staticmethod
	def _dec(dec):
		_hours = dec * 24
		hours = int(_hours)
		_minutes = 60 * (_hours - hours)
		minutes = int(_minutes)
		_seconds = 60 * (_minutes - minutes)
		seconds = int(_seconds)
		digit = int(10 * (_seconds - seconds))
		return f'{hours}h {minutes}m {seconds}.{digit}s'
