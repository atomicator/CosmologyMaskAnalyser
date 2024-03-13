import matplotlib.pyplot as plt
import numpy as np
from toolkit import toolkit

toolkit.plt_use_tex()

data = np.load("./bias_data/bias_test_planck_const_only.npy")

plt.hist(data[:, 0], bins=25)
plt.show()
print(data)
