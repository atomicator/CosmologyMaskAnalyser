import numpy as np
import matplotlib.pyplot as plt
from toolkit import toolkit

toolkit.plt_use_tex()

"""
Results:
Exact: galactic - 1.1977 +/- 0.017
       total    - 2.5142 +/- 0.025
       
N      | G      | T
-------|--------|-------
const  | 0.9922 | 2.3232
binary | 1.0226 | 2.3479
1      | 0.9839 | 2.3110
2      | 0.9867 | 2.3142
4      | 1.0159 | 2.3443
8      | 1.0565 | 2.3900
16     | 1.1091 | 2.4481
32     | 1.1433 | 2.4826
64     | 1.1774 | 2.5100
128    | 1.1934 | 2.5336
256    | 1.1915 | 2.5156
512    | 1.1780 | 2.4723
"""

fig = plt.figure()
ax = fig.add_subplot(111)

results = np.array([[0.9922, 2.3232],
                    [1.0226, 2.3479],
                    [0.9839, 2.3110],
                    [0.9867, 2.3142],
                    [1.0159, 2.3443],
                    [1.0565, 2.3900],
                    [1.1091, 2.4481],
                    [1.1433, 2.4826],
                    [1.1774, 2.5100],
                    [1.1934, 2.5336],
                    [1.1915, 2.5156],
                    [1.1780, 2.4723]])

galactic_mean = 1.1977
galactic_error = 0.017

total_mean = 2.5142
total_error = 0.025

galactic_data = np.abs(results[:, 0] - galactic_mean) / galactic_error
total_data = np.abs(results[:, 1] - total_mean) / total_error

n = np.linspace(-2, len(galactic_data) - 3, len(galactic_data))
x = np.power(2, n)

print(n)
print(x)

ax.plot(x, galactic_data, label="Galactic", color="red")
ax.plot(x, total_data, label="Total", color="blue")

ax.plot(x, np.ones(len(x)) * 3, linestyle="dotted", color="black", label=r"$\sigma = 3$")
ax.plot(x, np.ones(len(x)) * 2, linestyle="dashdot", color="black", label=r"$\sigma = 2$")
ax.plot(x, np.ones(len(x)), linestyle="dashed", color="black", label=r"$\sigma = 1$")

ax.legend()
#ax.set_yscale("log")
ax.set_xscale("log", base=2)

ax.set_xticks(x, ["Const", "Hemi"] + list(np.int_(x)[2:]))

ax.set_xlabel("NSIDE")
ax.set_ylabel(r"$\sigma$")
ax.set_ylim(0, 14)
plt.title("The effects of weighting the masked fraction using a binning algorithm")

plt.savefig("../../graphs/binned_results.png")

plt.show()
