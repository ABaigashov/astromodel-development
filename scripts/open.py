
from astropy.io import fits
import numpy as np
import cv2, os


path = os.path.abspath(r'.././21.09.2021/Orion/005.fits')
photo = fits.open(path)[0].data
photo = cv2.resize(photo, tuple(map((2).__rfloordiv__, photo.shape[::-1])))
cv2.imshow('astromodel', photo)
cv2.waitKey(0)