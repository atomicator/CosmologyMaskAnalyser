import numpy as np

a = np.array([True, True, True, False, False])
c = np.bool_(np.array([1, 1, 1, 0, 0]))

b = np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, 5]])

print(b[:, a])
print(b[:, c])
