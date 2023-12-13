import numpy
import numpy as np
from toolkit import toolkit

point = toolkit.load_mask("planck_point")
galactic = toolkit.load_mask("planck_galactic")

dif = point.map + 2 * galactic.map

print(np.array((np.sum(dif == 0), np.sum(dif == 1), np.sum(dif == 2), np.sum(dif == 3))) / len(dif))

data = 1 - numpy.real(numpy.loadtxt("../../sim_results/point_actual.csv", delimiter=","))

print(np.mean(data), np.std(data))
