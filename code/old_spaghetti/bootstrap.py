import numpy as np

bootstrap_samples = 10  # would be significantly higher for the actual simulations
n = 50000000
data = np.random.randint(0, 2, n)
print(len(data))
mean_estimates = []

for i in range(bootstrap_samples):
    print(f"Bootstrap iteration: {i}")
    new_data = np.random.choice(data, n)
    estimate = np.sum(new_data)
    print(estimate)
    mean_estimates.append(estimate)
