# Function to update scatter plot and display highlighted cell count
def update_scatter(total_counts_threshold, pct_mito_threshold, n_genes_by_counts_threshold):
    fig, ax = plt.subplots(figsize=(6, 6))

    # Plot all cells
    ax.scatter(adata.obs['total_counts'], adata.obs['n_genes_by_counts'], c='gray', label='All cells', alpha=0.5, s=2)

    # Highlight cells that meet the threshold criteria
    highlight = (
        (adata.obs['total_counts'] > total_counts_threshold) &
        (adata.obs['pct_counts_mt'] < pct_mito_threshold) &
        (adata.obs['n_genes_by_counts'] > n_genes_by_counts_threshold)
    )
    highlighted_cells = adata.obs[highlight]

    ax.scatter(highlighted_cells['total_counts'],
               highlighted_cells['n_genes_by_counts'], c='red', label='Highlighted cells',alpha =0.2, s=2)

    # Labels and plot formatting
    ax.set_xlabel('Total Counts')
    ax.set_ylabel('n_genes_by_counts')
    ax.set_title(f'Cells with Total Counts > {total_counts_threshold}, % Mito < {pct_mito_threshold}, and n_genes_by_counts > {n_genes_by_counts_threshold}')
    ax.legend()

    # Display the count of highlighted cells
    highlighted_count = len(highlighted_cells)
    print(f'Number of highlighted cells: {highlighted_count}')

    plt.show()
    
# Function to apply Louvain clustering and plot UMAP based on resolution
def update_cluster_plot(resolution,cluster='louvain'):
    # Perform Louvain clustering
    sc.tl.louvain(adata, resolution=resolution)
    # Plot the UMAP with Louvain clusters
    sc.pl.umap(adata, color=cluster)
