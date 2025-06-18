# Temporal trend analysis

Two-stage clustering pipeline for analyzing time-series data from aging experiments.

## What it does

This toolkit performs hierarchical clustering on aging data across multiple time points:
1. **Stage 1**: HDBSCAN clustering to identify major groups
2. **Stage 2**: K-means refinement on specific clusters for better granularity

## Requirements

```bash
pip install pandas numpy matplotlib seaborn scikit-learn hdbscan
```

## Usage

### Step 1: Initial HDBSCAN clustering
```bash
python habscan.py
```

**Input**: `aging_data.xlsx` with columns:
- `Name`: Compound names
- `*_Mean`: Mean peak areas at each time point (e.g., `W0_Mean`, `W2_Mean`, etc.)

**Outputs**:
- `hdbscan_cluster_result.xlsx`: Original data + cluster assignments
- `hdbscan_tsne_results.png`: t-SNE visualization of clusters
- `trend_cluster_*.png`: Z-score trends for each cluster

### Step 2: Refined clustering
```bash
python Superised cluster.py
```

Takes the HDBSCAN results and further subdivides Cluster 1 into 4(set based on previous research) sub-clusters.

**Input**: `hdbscan_cluster_result.xlsx` from Step 1

**Outputs**:
- `refined_cluster_result.xlsx`: Data with refined cluster assignments
- `refined_cluster_tsne.png`: t-SNE plot of refined clusters
- `refined_trend_plots/`: Folder containing trend plots for each refined cluster

## Key Parameters

### habscan.py
- `min_cluster_size=30`: Minimum compounds per cluster
- `min_samples=2`: Core point definition
- `cluster_selection_epsilon=0.8`: Cluster distance cutoff

### Superised cluster.py
- `n_clusters=4`: Number of sub-clusters for Cluster 1
- Custom time points on X-axis : `[0, 1, 2, 4, 8, 12]` weeks
- Highlight compounds (TN-wet only): Benzothiazole, 2-amino-benzothiazole

## Data Processing Pipeline

1. **Normalization**: Z-score normalization across time points
2. **Clustering**: HDBSCAN identifies dense regions 
3. **Refinement**: K-means splits large clusters for better resolution
4. **Visualization**: t-SNE projections and temporal trend plots

## Cluster Labeling

- `-1`: Noise points (HDBSCAN)
- `0, 1, 2...`: Original HDBSCAN clusters
- `10, 11, 12, 13`: Refined sub-clusters from original Cluster 1

## Customization

To adapt for your data:
- Adjust clustering parameters based on dataset size
- Modify time point labels in `custom_xticks`
- Change highlighted compounds for different experimental conditions
- Tune t-SNE parameters (`perplexity`, `learning_rate`) for better visualizations

## Notes

- Gray lines in plots = individual compounds
- Red lines = highlighted compounds (if applicable)
- Cluster -1 represents outliers/noise points
- All plots use Arial font for publication readiness
