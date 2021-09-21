
from telescope import Telescope
# # import numpy as np
# # import zwoasi
# # import cv2

with Telescope('COM5', timeout=1) as telescope:
	print(telescope.format_ra_dec(telescope.get_ra_dec()))
# # dll_path = __file__[::-1].split('\\', 1)[1][::-1] + r'\ASICamera2.dll'
# # zwoasi.init(library_file=dll_path)

# # camera = zwoasi.Camera(0)

# # print(camera.set_id(1.1a))

# # camera.close()

# # # a = np.zeros((600, 600), dtype=np.uint8)
# # # cv2.imshow('test', a)

# # # while True:
# # # 	key = cv2.waitKey(0)
# # # 	if key == 119:
# # # 		telescope.slew_ra(1)
# # # 	if key == 97:
# # # 		telescope.slew_dec(1)
# # # 	if key == 115:
# # # 		telescope.slew_ra(-1)
# # # 	if key == 100:
# # # 		telescope.slew_dec(-1)
# # # 	if key == 32:
# # # 		telescope.cancel_goto()
# # # 	if key == 27:
# # # 		break