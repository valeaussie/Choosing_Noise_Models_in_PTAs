import json
import matplotlib.pyplot as plt
import numpy as np

# Path to the JSON file saved earlier
json_path_1 = "/fred/oz005/users/vdimarco/Choosing_Noise_Models_in_PTAs/posteriors/mahalanobis_distances/mean_mahalanobis_distances_sanity_check_1.json"
json_path_2 = "/fred/oz005/users/vdimarco/Choosing_Noise_Models_in_PTAs/posteriors/mahalanobis_distances/mean_mahalanobis_distances_sanity_check_2.json"

# Load the JSON files
with open(json_path_1, 'r') as f:
    mean_distances_1 = json.load(f)
with open(json_path_2, 'r') as f:
    mean_distances_2 = json.load(f)

# Extract distances
distances_1 = np.array(list(mean_distances_1.values()))
#filtered_distances_1 = distances_1[distances_1 > 2.5]
distances_2 = np.array(list(mean_distances_2.values()))
#filtered_distances_2 = distances_2[distances_2 > 2.5]

# Convert the values to a list of floats (keys are strings in JSON)
distances_1 = [float(val) for val in mean_distances_1.values()]
distances_2 = [float(val) for val in mean_distances_2.values()]

# Plot
plt.figure(figsize=(10, 6))
plt.hist(distances_1, bins=50, alpha=0.6, edgecolor='#d3d3d3', density=False, label='Misspecified', color='steelblue')
plt.hist(distances_2, bins=50, alpha=0.6, edgecolor='#d3d3d3', density=False, label='Correctly specified', color='darkorange')

plt.xlabel('Mahalanobis distance', fontsize=16)
# Larger tick font sizes
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)

# Optionally apply log scale if there are long tails
plt.xscale('log')
plt.tight_layout()

# Save the plot
outdir = "/fred/oz005/users/vdimarco/Choosing_Noise_Models_in_PTAs/posteriors/mahalanobis_distances"
plt.savefig(f"{outdir}/mean_mahalanobis_distribution_plot_sanity_check_paper.pdf")