import matplotlib.pyplot as plt
import numpy as np


def plot(fig, ax, x, y):
    ax.plot(x, y)


def test(ax, points):
    ax.pcolormesh([[0, 1, 2], [0, 1, 2]])


if __name__ == "__main__":
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = np.linspace(0, 1, 3)
    y = np.linspace(0, 1, 3)
    x, y = np.meshgrid(x, y)

    c = ax.pcolormesh(x, y, [[1, 2], [3, 4]], alpha=1, cmap="winter")
    bar = fig.colorbar(c, orientation="horizontal", label="test")
    #bar.ax.set_ylabel("test", rotation="horizontal")
    plt.show()
