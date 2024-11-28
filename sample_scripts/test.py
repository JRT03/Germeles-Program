import numpy as np
import matplotlib.pyplot as plt

# Define parameters for the Gaussian
mean = 0  # Mean of the distribution
std_dev = 2.5  # Standard deviation

# Generate x values
x = np.linspace(-5, 5, 500)
# Calculate Gaussian function
y = (1 / (std_dev * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std_dev) ** 2)

# Create the plot
plt.figure(figsize=(8, 5))
plt.ylim(0.01,0.17)
plt.plot(x, y, color='black',linewidth=24)
plt.savefig('Gaussian_plume_image',dpi=300)
plt.show()
