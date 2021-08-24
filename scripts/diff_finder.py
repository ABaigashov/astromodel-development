
from astropy.io import fits
import numpy as np
import cv2, os


# [INFO]:
# | 
# | Use W A S D keys to control mask movement.
# | Use ENTER key to print info about state.
# ! 

# [USE "M13" AND "M51" FOR TESTS]:
# | 
# | BEST FOR "M51": [ΔE =  587813312; Δx =   35; Δy =   13;]
# | BEST FOR "M13": [ΔE =  580304428; Δx =   21; Δy =   -9;]
# ! 

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

x = y = 0

while True:
	image = frames[0]
	mask = frames[10]

	diff = 10 * np.abs(
		image.astype(np.int64) - np.roll(np.roll(mask, y, axis=0), x, axis=1).astype(np.int64)
	)

	sx = diff.shape[0] // 10
	sy = diff.shape[1] // 10

	diff = np.abs(diff - np.mean(diff))
	diff = diff.astype(np.uint16)

	cv2.imshow('astromodel', diff[sx*3:sx*7, sy*3:sy*6])
	
	try:
		key = chr(cv2.waitKey(100))
		if key in 'wsad':
			if key == 'w':
				y -= 1
			elif key == 's':
				y += 1
			elif key == 'a':
				x -= 1
			elif key == 'd':
				x += 1
		elif key == '\r':
			print(f'ΔE = {np.sum(diff):>10}; Δx = {x:>4}; Δy = {y:>4};')
	except:
		pass

	try:
		cv2.getWindowProperty('astromodel', 0)
	except:
		break