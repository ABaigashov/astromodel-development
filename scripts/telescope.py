
from typeguard import typechecked
from datetime import datetime
from serial import Serial
from enum import Enum
from typing import *


class TrackMode(Enum):
	OFF = 0
	ALT_AZ = 1
	EQATORIAL = 2
	PEC = 3


@typechecked
class Telescope:

	models = {
		0: 'EQ6 GOTO Series',
		1: 'HEQ5 GOTO Series',
		2: 'EQ5 GOTO Series',
		3: 'EQ3 GOTO Series',
		4: 'EQ8 GOTO Series',
		5: 'AZ-EQ6 GOTO Series',
		6: 'AZ-EQ5 GOTO Series',
		160: 'AllView GOTO Series'
	}
	models.update({id: 'AZ GOTO Series' for id in range(128, 144)})
	models.update({id: 'DOB GOTO Series' for id in range(144, 160)})

	def __init__(self, port: str, baudrate: int = 9600, timeout: Optional[int] = None) -> None:
		self.port = port
		self.baudrate = baudrate
		self.timeout = timeout
		self.serial = None

	def __enter__(self):
		self.serial = Serial(self.port, baudrate=self.baudrate, timeout=self.timeout)
		return self

	def __exit__(self, exc_type, exc_value, trace):
		self.serial.close()
		self.serial = None
		return False

	def execute(self, command: str) -> str:
		command = command.encode() + b'\r'
		if self.serial is None:
			with self as telescope:
				result = self._execute(telescope.serial, command)
		else:
			result = self._execute(self.serial, command)
		return result.decode()[:-1]

	@staticmethod
	def _execute(serial, command: bytes) -> bytes:
		serial.write(command)
		return serial.readline()

	@staticmethod
	def _decode_coordinate(coordinate: str) -> Tuple[float, float]:
		if len(coordinate) == 9:
			return tuple(int(coord, base=16) / 65536 for coord in coordinate.split(','))
		if len(coordinate) == 17:
			return tuple(int(coord[:-2], base=16) / 16777216 for coord in coordinate.split(','))

	@staticmethod
	def _encode_coordinate(coordinate: Tuple[float, float], precise: bool) -> str:
		if precise:
			return ','.join(hex(int(coord * 16777216))[2:].upper() + '00' for coord in coordinate)
		else:
			return ','.join(hex(int(coord * 65536))[2:].upper() for coord in coordinate)

	@staticmethod
	def _format_clock(value: float) -> str:
		_hours = value * 24
		hours = int(_hours)
		_minutes = 60 * (_hours - hours)
		minutes = int(_minutes)
		_seconds = 60 * (_minutes - minutes)
		seconds = int(_seconds)
		digit = int(10 * (_seconds - seconds))
		return f'{hours}h {minutes}m {seconds}.{digit}s'

	@staticmethod
	def _format_degrees(value: float) -> str:	
		_deg = value * 360
		deg = int(_deg)
		_minutes = 60 * (_deg - deg)
		minutes = int(_minutes)
		_seconds = 60 * (_minutes - minutes)
		seconds = int(_seconds)
		return f'{deg}° {minutes}’ {seconds}”'

	@classmethod
	def format_ra_dec(cls, coordinate: Tuple[float, float]) -> str:
		_ra, _dec = coordinate
		ra = cls._format_clock(_ra)
		dec = cls._format_degrees(_dec)
		return f'RA( {ra} ); DEC( {dec} )'

	@classmethod
	def format_azm_alt(cls, coordinate: Tuple[float, float]) -> str:
		_azm, _alt = coordinate
		azm = cls._format_degrees(_azm)
		alt = cls._format_degrees(_alt)
		return f'AZM( {azm} ); ALT( {alt} )'

	@staticmethod
	def format_datetime(time: List[int]) -> str:
		return datetime.fromisoformat(
			'20{5:02d}-{3:02d}-{4:02d} {0}:{1:02d}:{2:02d}.000{q}:00'.format(
				*time, q=(f'+{time[6]:02d}' if time[6] < 128 else f'-{(256 - time[6]):02d}')
			)
		).strftime('%I:%M:%S%p %B %d, %Y %Z')

	def format_location(self, location: List[int]) -> str:
		return '{4}°{5}’{6}” {8}, {0}°{1}’{2}” {9}'.format(*location, *map(
			lambda x: ['WE'[x[1]], 'SN'[x[1]]][x[0]], enumerate(location[3::4])
		))

	def get_ra_dec(self, precise: bool = True) -> Tuple[float, float]:
		coordinate = self.execute('e' if precise else 'E')
		return self._decode_coordinate(coordinate)

	def get_azm_alt(self, precise: bool = True) -> Tuple[float, float]:
		coordinate = self.execute('z' if precise else 'Z')
		return self._decode_coordinate(coordinate)

	def goto_ra_dec(self, ra: float, dec: float, precise: bool = True) -> str:
		coordinate = self._encode_coordinate((ra, dec), precise=precise)
		return self.execute(('r' if precise else 'R') + coordinate)

	def goto_azm_alt(self, azm: float, alt: float, precise: bool = True) -> str:
		coordinate = self._encode_coordinate((azm, alt), precise=precise)
		return self.execute(('b' if precise else 'B') + coordinate)

	def sync_ra_dec(self, ra: float, dec: float, precise: bool = True) -> str:
		coordinate = self._encode_coordinate((ra, dec), precise=precise)
		return self.execute(('s' if precise else 'S') + coordinate)

	def get_tracking(self) -> TrackMode:
		track_mode = ord(self.execute('t'))
		return TrackMode(track_mode)

	def set_tracking(self, mode: TrackMode) -> str:
		return self.execute('T' + chr(mode.value))

	def slew_ra(self, speed: int) -> str:
		assert speed in range(-9, 10)
		char = chr(abs(speed))
		if speed >= 0:
			command = 'P\x02\x10${char}\x00\x00\x00'
		else:
			command = 'P\x02\x10%{char}\x00\x00\x00'
		return self.execute(command)

	def slew_dec(self, speed: int) -> str:
		assert speed in range(-9, 10)
		char = chr(abs(speed))
		if speed >= 0:
			command = 'P\x02\x11${char}\x00\x00\x00'
		else:
			command = 'P\x02\x11%{char}\x00\x00\x00'
		return self.execute(command)

	def stop_slew(self) -> None:
		self.slew_ra(0)
		self.slew_dec(0)

	def get_location(self) -> List[int]:
		location = list(map(ord, self.execute('w')))
		return location

	def set_location(self, location: List[int]) -> str:
		return self.execute('W' + bytes(location).decode())

	def get_datetime(self) -> List[int]:
		time = list(map(ord, self.execute('h')))
		return time

	def set_datetime(self, time: List[int]) -> str:
		return self.execute('H' + bytes(time).decode())

	def get_version(self) -> str:
		version = map(lambda x: int(x, base=16), self.execute('V'))
		return '{0}{1}.{2}{3}.{4}{5}'.format(*version)

	def get_model(self) -> str:
		return self.models[ord(self.execute('m'))]

	def echo(self, message: str) -> str:
		return self.execute('K' + message)

	def is_alignment_complete(self) -> bool:
		return self.execute('J') == '\x01'

	def is_goto_in_progress(self) -> bool:
		return self.execute('J') == '1'

	def cancel_goto(self) -> str:
		return self.execute('M')

	def get_mount_state(self) -> str:
		return self.execute('p')