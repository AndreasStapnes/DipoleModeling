from methods import Dipole, Vec, spatial_span, cross, project_along
import numpy as np
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from earth import earth


xi = Vec((1,0,0))
eta = Vec((0,0,1))
zeta = cross(xi, eta)

pts = 100

rs = spatial_span(xi, eta, -1, 1, pts, -1, 1, pts)
xs, ys = np.meshgrid(np.linspace(-1,1,pts), np.linspace(-1,1,pts))
rs_sizes = np.apply_along_axis(np.linalg.norm, -1, rs)

import matplotlib.pyplot as plt
shape = np.shape(rs)[:-1]
B_vals = earth.B(rs)
#B_vals = np.apply_along_axis(np.linalg.norm, 2, B_vals)
#B_vals = np.where(rs_sizes > 0.1, np.log(B_vals), -15)
#B_vals = np.clip(B_vals, -1e-5,1e-5)
B_vals_size = np.apply_along_axis(np.linalg.norm, 2, B_vals)
B_vals_proj = project_along(B_vals, along=zeta, axis_tip=xi)
B_vals_size = np.clip(B_vals_size, -1e-5, 1e-5)

fig: Figure
ax: Axes
fig, (ax) = plt.subplots(1,1)
ax.set_facecolor("k")
#cmesh = ax.pcolormesh(xs, ys, B_vals, shading='auto')
#ax.pcolormesh(xs, ys, B_vals_size, cmap="inferno", shading="auto")
B_x, B_y = B_vals_proj.transpose((2,0,1))
quiver = ax.streamplot(np.linspace(-1, 1, pts), np.linspace(-1, 1, pts), B_x, B_y, density=1.4, color=B_vals_size,cmap="inferno")

#cbarax = fig.add_axes([0.55, 0.8, 0.2, 0.05])
#ax.contourf(xs, ys, B_vals_size)
ax.set_aspect("equal")
#plt.colorbar(cmesh, cbarax, ax, orientation="horizontal")
fig.show()
