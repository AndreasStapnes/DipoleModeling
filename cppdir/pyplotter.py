import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
import matplotlib.cm as cm

sublines = 40


with open("out.txt", "r") as file:
    ax: Axes
    vals = [[[float(num) for num in elem.strip("()\n").split(",")] for elem in line.split(") (")] for line in file.readlines()]
    vals = np.array(vals)
    positions, velocities = vals.transpose((1,0,2))

    N = len(vals)
    n = int(N/sublines)
    fig = plt.figure()
    ax = fig.add_subplot(111,projection="3d")

    convergence_fig, convergence_ax = plt.subplots(1,1)

    for i in range(sublines-1):
        relevant_positions = positions[n*i:n*(i+1)]
        ax.plot(*relevant_positions.T, color=cm.plasma(i/sublines))
    ax.plot(*positions[n*(sublines-1):].T, color=cm.plasma((sublines-1)/sublines))
    ax.scatter(*positions[-1], color="k")

    convergence_ax.plot(np.linalg.norm(velocities, axis=1))
    convergence_fig.show();


    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    zlim = ax.get_zlim()
    x_mean = np.mean(xlim); y_mean=np.mean(ylim); z_mean = np.mean(zlim);
    x_diff = xlim[1]-xlim[0]; y_diff = ylim[1]-ylim[0]; z_diff = zlim[1] - zlim[0];
    rad = np.max((x_diff, y_diff, z_diff))/2;
    ax.set_xlim([x_mean-rad, x_mean+rad])
    ax.set_ylim([y_mean-rad, y_mean+rad])
    ax.set_zlim([z_mean-rad, z_mean+rad])
    ax.set_xlabel("x[m]")
    ax.set_ylabel("y[m]")
    ax.set_zlabel("z[m]")
    fig.show()
