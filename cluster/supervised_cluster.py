import pandas as pd
import os
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
import seaborn as sns

# set font size and style
plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.labelsize"] = 26
plt.rcParams["xtick.labelsize"] = 24
plt.rcParams["ytick.labelsize"] = 24


# input file from the HDBSCAN
file_path = "hdbscan_cluster_result.xlsx"
data = pd.read_excel(file_path)

# extract the mean
expr_data = data.filter(regex="_Mean$")
expr_matrix = expr_data.values
time_points = expr_data.columns
compound_names = data["Name"]

# zscore normalized
expr_z = StandardScaler().fit_transform(expr_matrix.T).T
zscore_df = pd.DataFrame(expr_z, columns=time_points, index=compound_names)

# HDBSCAN_Cluster
mask_cluster1 = data["HDBSCAN_Cluster"] == 1
expr_cluster1 = expr_z[mask_cluster1]

# set the cluster number under supervised
kmeans = KMeans(n_clusters=4, random_state=42)
sub_labels = kmeans.fit_predict(expr_cluster1)

#  RefinedCluster
data["RefinedCluster"] = data["HDBSCAN_Cluster"]
data.loc[mask_cluster1, "RefinedCluster"] = sub_labels + 10


# plot
output_folder = "refined_trend_plots"
os.makedirs(output_folder, exist_ok=True)

refined_clusters = sorted(data["RefinedCluster"].unique())

# set the target compounds in specific condition
highlight_condition = "condition"
highlight_compounds = ["mz_544.20846_RT_25.621"]
is_tn_wet = True

# set the x-axis coordinate
custom_xticks = [0, 2, 4, 6, 8, 12]

for cluster_id in refined_clusters:
    members = data[data["RefinedCluster"] == cluster_id]
    if len(members) == 0:
        continue

    zscores = zscore_df.loc[members["Name"]]


    # plots
    plt.figure(figsize=(10, 6))
    for compound in zscores.index:
        is_highlight = is_tn_wet and compound in highlight_compounds
        plt.plot(custom_xticks, zscores.loc[compound],
                 color='red' if is_highlight else 'gray',
                 linewidth=2 if is_highlight else 1,
                 alpha=1.0 if is_highlight else 0.3)
    plt.xticks(custom_xticks)
    plt.xlim(0,12)
    plt.xlabel("Weeks")
    plt.ylabel("Z-score")
    plt.tight_layout()

    plt.savefig(f"{output_folder}/trend_cluster_{cluster_id}.png")
    plt.close()

# save the new cluster results
output_excel = "refined_cluster_result.xlsx"
data.to_excel(output_excel, index=False)

# using t-SNE RefinedCluster
tsne = TSNE(n_components=2, random_state=42, perplexity=30, learning_rate=200)
tsne_proj = tsne.fit_transform(expr_z)

# set the color
unique_clusters = sorted(data["RefinedCluster"].unique())
palette = sns.color_palette("hls", len([x for x in unique_clusters if x != -1]))
color_map = {}

for i, label in enumerate(unique_clusters):
    if label == -1:
        color_map[label] = (0.5, 0.5, 0.5)
    else:
        color_map[label] = palette[i if -1 not in unique_clusters else i - 1]

# plots
plt.figure(figsize=(12, 8))
for label in unique_clusters:
    mask = data["RefinedCluster"] == label
    plt.scatter(tsne_proj[mask, 0], tsne_proj[mask, 1],
                label=f"Cluster {label}",
                color=color_map[label],
                alpha=0.6, s=30)

plt.title("t-SNE Projection of Refined Clusters (Gray = Noise)")
plt.legend(loc="upper right", bbox_to_anchor=(1.25, 1))
plt.grid(True)
plt.tight_layout()
plt.savefig("refined_cluster_tsne.png")
plt.show()


