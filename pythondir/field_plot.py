from methods import Dipole, Vec, spatial_span, cross, project_along
import numpy as np
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from earth import earth

#one wishes to plot the B-field projection in a plane which contains the dipole.
xi = Vec((1,0,0))       # first spanning vector of the plane
eta = Vec((0,0,1))      # second spanning vector of the plane
xi_label = r"$x$[m]"
eta_label = r"$z[m]$"
filename="xz_proj"

zeta = cross(xi, eta)   # normal vector of the plane

pts = 100               # amt of points in either direction of the plane one wishes to plot

view_r = 16e6            # distance from the dipole in the xi and eta-direction one wishes to plot
max_B = 0.64e-4            # the maximal B-field-strength one will allow to plot

rs = spatial_span(xi, eta, -view_r, view_r, pts, -view_r, view_r, pts)
            # meshgrid of the points for which one wants to plot B
xs, ys = np.meshgrid(np.linspace(-view_r,view_r,pts), np.linspace(-view_r,view_r,pts))
            # 2d-points in the xi-eta-basis of the relevant points to plot


import matplotlib.pyplot as plt
shape = np.shape(rs)[:-1]
B_vals = earth.B(rs)    #Calculate the B-field at each r-point of the plane



B_vals_proj = project_along(B_vals, along=zeta, axis_tip=xi)
B_vals_size = np.apply_along_axis(np.linalg.norm, 2, B_vals)
B_vals_size = np.clip(B_vals_size, -max_B, max_B)

fig: Figure
ax: Axes
fig, (ax) = plt.subplots(1,1)
ax.set_facecolor("k")
cmesh = ax.pcolormesh(xs, ys, B_vals_size*1e4, shading='auto', cmap="cividis")
#ax.pcolormesh(xs, ys, B_vals_size, cmap="inferno", shading="auto")
B_x, B_y = B_vals_proj.transpose((2,0,1))

quiver = ax.streamplot(np.linspace(-view_r, view_r, pts), np.linspace(-view_r, view_r, pts), B_x, B_y, density=1.4, color="k")
earth_circle = plt.Circle((0,0),6.371e6, color="k")
cbarax = fig.add_axes([0.85, 0.58, 0.05, 0.3])
cbarax.set_title("B[gauss]")
#ax.contourf(xs, ys, B_vals_size)
ax.add_patch(earth_circle)
ax.set_aspect("equal")
plt.colorbar(cmesh, cbarax, ax, orientation="vertical")
ax.set_ylabel(eta_label)
ax.set_xlabel(xi_label)
fig.show()
fig.savefig(f"{filename}.pdf")
