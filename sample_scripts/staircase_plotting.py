import matplotlib.pyplot as plt
import numpy as np

# 1. Create some data
x = np.linspace(0, 10, 100)
y_data = [
    np.sin(x),            # Data for the 1st plot
    np.cos(x),            # Data for the 2nd plot
    np.tan(x),            # Data for the 3rd plot
    np.exp(-x),           # Data for the 4th plot
]

titles = ['Sine wave', 'Cosine wave', 'Tangent wave', 'Exponential decay']

# 2. Create a 2x2 grid of subplots
fig, axes = plt.subplots(2, 2, figsize=(10, 8))

# Flatten the axes array to make it easier to iterate
axes = axes.flatten()

# 3. Use a loop to fill each subplot in the next available space
for i in range(len(axes)):
    print(type(axes[i]))
    axes[i].plot(x, y_data[i])      # Plot the corresponding data
    axes[i].set_title(titles[i])    # Set the corresponding title

# 4. Adjust layout to prevent overlap
plt.tight_layout()

# 5. Show the plot
plt.show()
