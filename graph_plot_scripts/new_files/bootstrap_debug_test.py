import numpy as np

std = 3.4562e-5

print(round(std, +1 - int(np.floor(np.log10(std)))))
