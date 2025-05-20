import json
import matplotlib.pyplot as plt
import numpy as np

# Path to the JSON file saved earlier
json_path = "/fred/oz005/users/vdimarco/Choosing_Noise_Models_in_PTAs/posteriors/mahalanobis_distances/mean_mahalanobis_distances.json"

# Load the JSON file
with open(json_path, 'r') as f:
    mean_distances = json.load(f)

# Extract distances > 2 (these are the chains that failed)
distances = np.array(list(mean_distances.values()))
filtered_distances = distances[distances > 2.5]

# Convert the values to a list of floats (keys are strings in JSON)
distances = [float(val) for val in mean_distances.values()]

# Plot the distribution
plt.figure(figsize=(8, 5))
plt.hist(filtered_distances, bins=18, edgecolor='black', density=True, alpha=0.7)  # use this instead for just histogram

plt.xlabel('Mean Mahalanobis distance')
plt.title(f'Distribution of mean Mahalanobis distances')
plt.tight_layout()

# Save the plot
outdir = "/fred/oz005/users/vdimarco/Choosing_Noise_Models_in_PTAs/posteriors/mahalanobis_distances"
plt.savefig(f"{outdir}/mean_mahalanobis_distribution_plot.png")
plt.savefig(f"{outdir}/mean_mahalanobis_distribution_plot.pdf")
