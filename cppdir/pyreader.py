import matplotlib.pyplot as plt
import numpy as np

with open("out.txt", "r") as file:
    vals = [[float(num) for num in line.strip("()\n").split(",")] for line in file.readlines()]
    vals = np.array(vals).T
    fig = plt.figure()
    ax = fig.add_subplot(111,projection="3d")
    ax.plot(*vals)
    fig.show()
