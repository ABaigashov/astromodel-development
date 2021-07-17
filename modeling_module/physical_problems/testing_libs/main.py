
import os, sys
sys.path.append(os.path.abspath(os.path.join('..', '..')))
sys.path.append(os.path.abspath('modeling_module'))

from configurator import Configurator
configuration = Configurator._decompress(sys.argv[1], sys.argv[2])


RESULT = ''
def log(text='', end='\n'):
	global RESULT
	RESULT += str(text) + end
	print(text, end=end)

from time import sleep
leopart = mshr = dolfin = fenics = False
log('Testing started\n')

if configuration.fenics:
	try:
		log(end='Testing "fenics"... ')
		sleep(5)
		print(end='\b')
		import fenics
		log('Done')
		fenics = True
	except ImportError:
		log('Fail')

if configuration.dolfin:
	try:
		log(end='Testing "dolfin"... ')
		sleep(5)
		print(end='\b')
		import dolfin
		log('Done')
		dolfin = True
	except ImportError:
		log('Fail')

if configuration.mshr:
	try:
		log(end='Testing "mshr"... ')
		sleep(5)
		print(end='\b')
		import mshr
		log('Done')
		mshr = True
	except ImportError:
		log('Fail')

if configuration.leopart:
	try:
		log(end='Testing "leopart"... ')
		sleep(5)
		print(end='\b')
		import leopart
		log('Done')
		leopart = True
	except ImportError:
		log('Fail')

if all([leopart, mshr, dolfin, fenics]):
	log('\nTest passed')
else:
	log('\nTest failed')


output = configuration.OUTPUT + '.log'
with open(output, 'wb') as f:
	f.write(RESULT.encode())


print(sys.argv[3], output)
