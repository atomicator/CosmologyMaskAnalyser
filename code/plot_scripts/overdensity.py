from toolkit import toolkit
import numpy as np
import matplotlib.pyplot as plt

toolkit.plt_use_tex()

error_bar_colors = ["xkcd:aqua blue", "orange", "xkcd:mint green", "pink", "xkcd:burnt sienna"]
line_colors = ["xkcd:electric blue", "red", "xkcd:grass green", "purple", "xkcd:reddish brown"]

#f_s_values = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
#f_s_values = [0.01, 0.02, 0.04, 0.06, 0.08]
f_s_values = [0.009859934289099422,]

excess = np.linspace(0, 0.4, 5000)

f_s_values = list(f_s_values)
#f_s_values = f_s_values[0:1]

for f_s in f_s_values:
    #overdensity = (excess * f_s) / (f_s * (1 - (excess + 1) * f_s))
    #n = 100 * (((overdensity + 1) * f_s) + (1 - f_s) - 1)
    #n = 100 * (overdensity) * f_s
    #hidden_fraction = 100 * (overdensity * f_s) / (overdensity * f_s + 1)
    f_c = (excess + 1) * f_s
    overdensity = (f_c - f_s) / (f_s * (1 - f_c))
    hidden_fraction = overdensity * f_s
    plt.plot(excess, f_c * 100, color=line_colors[f_s_values.index(f_s)],
             label=r"$f_{s}" + rf"= {str(f_s * 100)[0:4]}\%$")
             #label="Excess correction")
    #plt.plot(excess, np.ones(np.shape(excess)) * f_s * 100, linestyle="dotted", color=error_bar_colors[f_s_values.index(f_s)],
    #         #label=rf"Excess $\times {f_s}$")
    #         label="Uncorrelated")

plt.plot(excess, np.ones(np.shape(excess)) * f_s_values[0] * 100, linestyle="dashed", color="k")
plt.xlim(0, 0.4)
plt.ylim(0, 1.5)
plt.xlabel("Excess")
plt.ylabel(r"$f_{c}$ (\%)")
plt.title("Number count correction as a function of excess for a masked sky fraction of 2\%")
plt.legend()
plt.show()

# Excess = 0.074 +/- 0.015
