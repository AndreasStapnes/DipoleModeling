from earth import earth

from methods import ChargedParticle, Vec

import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111,projection="3d")

pos = Vec((-10e6, 0, 0))
vel = Vec((1e6, 0, 3e5))
cp = ChargedParticle(pos, vel, 1e-19, 1.66e-27)
cp.reacting_dipoles.append(earth)
poses = []

time = 10
timestep = 1e-5
for i in range(int(time/timestep)):
    cp % timestep
    poses.append(cp.position)
#print(poses)


ax.plot(*np.array(poses).T)
#ax.scatter(*earth.position, color="red", s=50)

x_limits = ax.get_xlim3d()
y_limits = ax.get_ylim3d()
z_limits = ax.get_zlim3d()

x_range = abs(x_limits[1] - x_limits[0])
x_middle = np.mean(x_limits)
y_range = abs(y_limits[1] - y_limits[0])
y_middle = np.mean(y_limits)
z_range = abs(z_limits[1] - z_limits[0])
z_middle = np.mean(z_limits)

# The plot bounding box is a sphere in the sense of the infinity
# norm, hence I call half the max range the plot radius.
plot_radius = 0.5*max([x_range, y_range, z_range])

ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

fig.show()