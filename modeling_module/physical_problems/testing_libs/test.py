
from traceback import format_exception as error_text

import os, sys
sys.dont_write_bytecode = True
sys.path.append(os.path.abspath(os.path.join('..', '..')))

from jobexecutor import Dispatcher

if len(sys.argv) > 1:
	configuration_name = sys.argv[1]
	if configuration_name.find('.json') < 0:
		configuration_name += '.json'
else:
	configuration_name = 'test.json'

problem_dir = os.path.dirname(__file__)
configuration = os.path.join(problem_dir, 'init_files', configuration_name)
output = os.path.join(problem_dir, 'results')

dispatcher = Dispatcher(configuration, output)
dispatcher.init()

result = dispatcher.run()
if isinstance(result, tuple):
	print(''.join(error_text(*result)))
