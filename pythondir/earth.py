import numpy as np
from methods import Dipole

#---------------CONSTANTS--------------------

tilt = 23.7     #degrees
m = -8.22e22       #A m^2
azim = 0

#---------------------------------------------

theta = tilt * np.pi/180
direction = np.array([np.sin(theta)*np.cos(azim), np.sin(theta)*np.sin(azim), np.cos(theta)])


earth = Dipole(direction*m)
