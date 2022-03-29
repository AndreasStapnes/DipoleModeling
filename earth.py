import numpy as np
from methods import Dipole

#---------------CONSTANTS--------------------

tilt = 23.7 #degrees
m = 4       #A m^2

#---------------------------------------------

theta = tilt * np.pi/180
direction = np.array([np.sin(theta), 0, np.cos(theta)])


earth = Dipole(direction*m)
