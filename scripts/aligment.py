
from telescope import Telescope
import numpy as np
import zwoasi
import cv2

with Telescope('COM9') as telescope:
	print(telescope.format_ra_dec(telescope.get_ra_dec()))

zwoasi.init(library_file=__file__[::-1].split('\\',1)[1][::-1]+'\\ASICamera2.dll')
camera = zwoasi.Camera(0)

camera.set_control_value(zwoasi.ASI_GAIN, 150)
camera.set_control_value(zwoasi.ASI_EXPOSURE, 30000)
camera.set_control_value(zwoasi.ASI_WB_B, 99)
camera.set_control_value(zwoasi.ASI_WB_R, 75)
camera.set_control_value(zwoasi.ASI_GAMMA, 50)
camera.set_control_value(zwoasi.ASI_BRIGHTNESS, 50)
camera.set_control_value(zwoasi.ASI_FLIP, 0)
camera.set_image_type(zwoasi.ASI_IMG_RAW16)


def load_image():
	while True

def get_start_image():
	whbi = camera.get_roi_format()
	if whbi[3] == zwoasi.ASI_IMG_RAW8 or whbi[3] == zwoasi.ASI_IMG_Y8:
	    return np.zeros([whbi[1], whbi[0]], dtype=np.uint8)
	elif whbi[3] == zwoasi.ASI_IMG_RAW16:
	    return np.zeros([whbi[1], whbi[0]], dtype=np.uint16)
	elif whbi[3] == zwoasi.ASI_IMG_RGB24:
	    return np.zeros([whbi[1], whbi[0], 3], dtype=np.uint8)

IMAGE = get_start_image()


while True:
	cv2.imshow('astromodel-development', IMAGE)
	key = cv2.waitKey(0)
	if key == 119:
		telescope.slew_ra(1)
	if key == 97:
		telescope.slew_dec(1)
	if key == 115:
		telescope.slew_ra(-1)
	if key == 100:
		telescope.slew_dec(-1)
	if key == 32:
		telescope.cancel_goto()
	if key == 27:
		break