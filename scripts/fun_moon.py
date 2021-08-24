
from astropy.io import fits
import cv2, os, time, math
import numpy as np

path = os.path.abspath('../23.08.2021/Moon')

frames = np.array([
	(
		lambda x: (x[0].data, x.close())[0]
	)(
		fits.open(os.path.join(path, name))
	)
	for name in os.listdir(path)
	if name[-5:] == '.fits'
])

size = (np.array(frames[0].shape) / 2).astype(np.int)

while True:
	a = math.pi * 2 * ((time.time() % 15) / 15)
	x0 = int(100 * math.cos(a))
	x = -2 * x0
	y0 = int(100 * math.sin(a))
	y = -2 * y0

	image = frames[0]
	mask = frames[1]

	diff = np.roll(np.roll(
		image.astype(np.int16) - np.roll(np.roll(mask, x, axis=0), y, axis=1).astype(np.int16), x0, axis=0
	), y0, axis=1).astype(np.uint16)

	cv2.imshow('astromodel', cv2.resize(diff, tuple(size)[::-1], interpolation=cv2.INTER_CUBIC))

	cv2.waitKey(10)

	try:
		cv2.getWindowProperty('astromodel', 0)
	except:
		break