import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import hdbscan
from sklearn.manifold import TSNE
import os

# set your path
file_path = "aging_data.xlsx"

# input excel
data = pd.read_excel(file_path)

# extract_Mean as the data
expr_data = data.filter(regex="_Mean$")
expr_matrix = expr_data.values
compound_names = data["Name"]

#  Z-score normalized
expr_z = StandardScaler().fit_transform(expr_matrix.T).T

#  HDBSCAN cluster parameter setting
clusterer = hdbscan.HDBSCAN(
    min_cluster_size=30,
    min_samples=2,
    cluster_selection_epsilon=0.8,
    metric='euclidean'
)
labels = clusterer.fit_predict(expr_z)

# add labels
data["HDBSCAN_Cluster"] = labels
data["AveragePeakArea"] = expr_matrix.mean(axis=1)

# output the results
print("Cluster counts:\n", data["HDBSCAN_Cluster"].value_counts())

# save the cluster results in Excel
output_excel = os.path.join(os.path.dirname(file_path), "hdbscan_cluster_result.xlsx")
data.to_excel(output_excel, index=False)
print(f"saved in {output_excel}")

# projection（t-SNE）
projection = TSNE().fit_transform(expr_z)
palette = sns.color_palette("hls", len(np.unique(labels)))

plt.figure(figsize=(10, 6))
for cluster_id in np.unique(labels):
    mask = labels == cluster_id
    color = "gray" if cluster_id == -1 else palette[cluster_id % len(palette)]
    plt.scatter(projection[mask, 0], projection[mask, 1],
                label=f"Cluster {cluster_id}", alpha=0.6, color=color)

plt.title("HDBSCAN Clustering (t-SNE projection)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(os.path.dirname(file_path), "hdbscan_tsne_results.png"))
plt.show()

# generalized plots
time_points = expr_data.columns
data["HDBSCAN_Cluster"] = data["HDBSCAN_Cluster"].astype(str)

for cluster_id in data["HDBSCAN_Cluster"].unique():
    cluster_data = expr_z[data["HDBSCAN_Cluster"] == cluster_id]
    if cluster_data.shape[0] == 0:
        continue

    plt.figure(figsize=(10, 6))
    for row in cluster_data:
        plt.plot(time_points, row, color='gray', alpha=0.3)
    plt.title(f"Z-score Time Trends - Cluster {cluster_id} ({cluster_data.shape[0]} compounds)")
    plt.xticks(rotation=45)
    plt.ylabel("Z-score")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(os.path.dirname(file_path), f"trend_cluster_{cluster_id}.png"))
    plt.show()
