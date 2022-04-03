import numpy as np
from methods import Dipole

#---------------CONSTANTS--------------------

tilt = 23.7     #degrees
m = -8.22e22       #A m^2

#---------------------------------------------

theta = tilt * np.pi/180
direction = np.array([0, np.sin(theta), np.cos(theta)])


earth = Dipole(direction*m)
