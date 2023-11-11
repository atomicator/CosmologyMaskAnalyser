import numpy as np
from toolkit import toolkit, bootstrapping

mask = toolkit.PixellMask("../../data/ACT_mask.fits", hdu=1)
#mask = toolkit.HealpyMask("../../data/planck_galactic_mask.fits")

data = []
full_data = []

initial_points = 1000
depth = 1000
iterations = 1
total_points = initial_points * depth * iterations

for i in range(initial_points):
    x = mask.gen_random_point()
    data.append(x)

for i in range(total_points):
    x = mask.gen_random_point()
    full_data.append(x)

mean, std = bootstrapping.bootstrapping(data, depth, iterations)

exact_unmasked_fraction = mask.calc_exact_unmasked_fraction()

print("Results:")
print(f"Exact results:\nUnmasked fraction: {exact_unmasked_fraction}")
print(f"Original data: \nPoints: {len(data)}; Mean: {np.mean(data)}; Deviation: {np.std(data)}; Error (calculated): "
      f"{np.std(data) / np.sqrt(len(data))}; Error (exact): {np.abs(np.mean(data) - exact_unmasked_fraction)}; Error "
      f"(predicted): {np.sqrt(exact_unmasked_fraction * (1 - exact_unmasked_fraction) / len(data))}")
print(f"Bootstrapped data: \nDepth: {depth}; Iterations:{iterations}; Mean: {mean}; Deviation: {std}; "
      f"Error (calculated): {std / np.sqrt(depth * iterations)}; Error (exact): "
      f"{np.abs(mean - exact_unmasked_fraction)}")
print(f"Full data: \nPoints: {len(full_data)}; Mean: {np.mean(full_data)}; Deviation: {np.std(full_data)}; Error (calculated): "
      f"{np.std(full_data) / np.sqrt(len(full_data))}; Error (exact): "
      f"{np.abs(np.mean(full_data) - exact_unmasked_fraction)}; Error "
      f"(predicted): {np.sqrt(exact_unmasked_fraction * (1 - exact_unmasked_fraction) / len(full_data))}")
print("""Calculated errors are estimates of the error from the data alone, exact errors are the difference between the
calculated and exact values. Predicted errors are the errors predicted from a binomial distribution. Depth refers 
to the amount of times the data was sampled, iterations refers to how many times the data was reset.""")

#mean_estimates = []

#for i in range(int(total_points / initial_points)):
#    mean_estimates.append(np.mean(full_data[i * initial_points:(i+1)*initial_points]))
"""
print(mean_estimates)

print("\n")

mean = np.mean(mean_estimates)
print(mean)
print("\n")
std = np.std(mean_estimates)
error = std / np.sqrt(len(mean_estimates))
print(mean, std, error)

plt.title(rf"Mean = {round(mean, 1 - int(np.floor(np.log10(error))))} $\pm$ {round(error, 1 - int(np.floor(np.log10(error))))}")
plt.hist(mean_estimates, color="red")
plt.xlabel("Mean estimates")
plt.ylabel("Frequency")
plt.show()
"""