
from astropy.io import fits
import numpy as np
import cv2, os


path = os.path.abspath('../23.08.2021/M13/001.fits')
photo = fits.open(path)[0].data

cv2.imshow('astromodel', photo)
cv2.waitKey(0)