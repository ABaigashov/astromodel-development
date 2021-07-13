from uuid import uuid4 as random_hash
from sympy import sympify
import json, os, pickle
import numpy as np


def _to_float(value):
	value = value.replace(',', '.')
	# some other stuff to do
	return float(eval(value))

def _parse_expr(value):
	return sympify(value)

def _to_int(value):
	return int(_to_float(value))


class ABS_RecursiveContainer:

	def __init__(self, connection):
		self.storage = {}
		for key, value in connection.items():
			if isinstance(value, dict):
				value = ABS_RecursiveContainer(value)
			if isinstance(value, list):
				temp = []
				for slot in value:
					if isinstance(slot, dict):
						slot = ABS_RecursiveContainer(slot)
					temp.append(slot)
				value = temp
			self.storage[key] = value

	def __getattr__(self, key):
		if key == 'storage':
			return object.__getattribute__(self, 'storage')
		return self.storage[key]

	def __setattr__(self, key, value):
		if key == 'storage':
			return object.__setattr__(self, 'storage', {})
		raise AttributeError

	def __getitem__(self, key):
		return self.storage[key]

	def __repr__(self):
		return f'ABS_RecursiveContainer({self.storage})'


class Configurator:

	def __init__(self, data):
		self._CFG_RAW_DATA = data
		

	@classmethod
	def _decompress(cls, pickled_data, output):
		self = cls(pickle.loads(bytes.fromhex(pickled_data)))
		self.OUTPUT = os.path.join(output, random_hash().hex[:16])
		return self

	@staticmethod
	def to_type(value, case):
		if case['type'] == 'str':
			return value
		elif case['type'] == 'bool':
			return value.lower() in ['true', '1', 'on', 'y']
		elif case['type'] == 'float':
			return _to_float(value)
		elif case['type'] == 'int':
			return _to_int(value)
		elif case['type'] == 'expr':
			return _parse_expr(value)
		raise TypeError(f'type <{case["type"]}> is not defined')

	@classmethod
	def fill_defaults(cls, name, value, defaults, parameters, config):
		if 'units' in config['CASES'][name]:
			try:
				units = cls.to_type(value['units'], config['CASES'][name])
			except:
				units = min([max(v, v**-1) for v in map(eval, config['UNITS'][config['CASES'][name]['units']].values())])
			if config['CASES'][name]['dementional']:
				stack = []
				for v in value['value']:
					try:
						stack.append(cls.to_type(v, config['CASES'][name]))
					except:
						stack.append(defaults[name])
				return units * np.array(stack)
			try:
				return units * cls.to_type(value['value'], config['CASES'][name])
			except:
				return units * defaults[name]
		try:
			return cls.to_type(value, config['CASES'][name])
		except:
			return defaults[name]


	@classmethod
	def parse_general(cls, defaults, parameters, config):
		raw = {}
		for name in set(parameters['GENERAL']) | set(config['GENERAL']):
			if name in defaults:
				raw[name] = cls.fill_defaults(name, parameters['GENERAL'].get(name), defaults, parameters, config)
			else:
				raw[name] = parameters['GENERAL'][name]
		return raw

	@classmethod
	def parse_objects(cls, defaults, parameters, config):
		raw = {}
		if 'OBJECTS' not in config:
			return raw
		for name in config['OBJECTS']:
			raw[name] = []
		for name, objects in parameters['OBJECTS'].items():
			for objekt in objects:
				raw[name].append({})
				for slot in set(config['OBJECTS'][name]['cases']) | set(objekt.keys()):
					if slot in defaults:
						raw[name][-1][slot] = cls.fill_defaults(slot, objekt.get(slot), defaults, parameters, config)
					else:
						raw[name][-1][slot] = objekt[slot]
		return raw

	@classmethod
	def generate_defaults(cls, parameters, config):
		defaults = {}
		for name, case in config['CASES'].items():
			default = case['default']
			defaults[name] = cls.to_type(default, case)
		return defaults

	@classmethod
	def parse_parameters(cls, parameters):
		problem = parameters['PROBLEM']
		if os.getcwd() == '/usr/src':
			config_path = os.path.join('.', 'modeling_module', 'physical_problems', problem, 'config.json')
		else:
			config_path = 'config.json'
		with open(config_path, 'rb') as f:
			config = json.load(f)
		defaults = cls.generate_defaults(parameters, config)
		general = cls.parse_general(defaults, parameters, config)
		objects = cls.parse_objects(defaults, parameters, config)
		return problem, dict(general=general, objects=objects)

	def __getattr__(self, key):
		if key in {'_CFG_RAW_DATA', 'OUTPUT'}:
			return object.__getattribute__(self, key)
		if key in self._CFG_RAW_DATA['general']:
			return self._CFG_RAW_DATA['general'][key]
		if key in self._CFG_RAW_DATA['objects']:
			value = ABS_RecursiveContainer(self._CFG_RAW_DATA['objects'])[key]
			return value
		raise KeyError(key) from None

	def __getitem__(self, key):
		return self.__getattr__(key)
