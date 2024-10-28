# Single-Cell RNA-seq Analysis Workshop

## Overview

This workshop provides a comprehensive guide to analyzing single-cell RNA sequencing (scRNA-seq) data, covering essential steps from quality control to cell type identification. Participants will use Python libraries such as `scanpy` for data preprocessing, visualization, and analysis.

This workshop is part of **6120 Precision and Genomic Medicine**, taught by **Iwijn De Vlaminck**.

## Submodules Covered

1. **Environment Setup**: Install the necessary Python packages (`scanpy`, `louvain`, `ipywidgets`) to prepare your computational environment for scRNA-seq analysis.
2. **Data Preparation**: Load datasets, create directories, and retrieve gene annotations for proper data organization and accessibility.
3. **Quality Control**: Perform quality checks on the cells by filtering based on mitochondrial content, total counts, and gene counts. Use Scrublet to detect potential doublets and remove them.
4. **Normalization and Filtering**: Normalize the scRNA-seq data to mitigate batch effects, extract highly variable genes, and filter out sex-specific genes to refine the dataset.
5. **Dimensionality Reduction**: Apply Principal Component Analysis (PCA) to reduce data dimensionality and Uniform Manifold Approximation and Projection (UMAP) for visualizing the high-dimensional structure of the data.
6. **Clustering and Differential Expression**: Use Louvain clustering to identify distinct cell populations, and determine differentially expressed marker genes. Interactive widgets allow for adjusting clustering parameters.
7. **Cell Type Annotation**: Annotate cell clusters using known marker genes to assign cell type identities.
8. **Pseudotime Analysis**: Infer pseudotime to explore cellular differentiation trajectories, providing insight into the dynamic processes occurring within the cell populations.

## Getting Started

### Prerequisites

- **Python 3.6+**
- **Jupyter Notebook** or **Google Colab**

### Installing Dependencies

To install the required packages, run the following command:

```python
!pip install scanpy python-igraph louvain ipywidgets pybiomart
```

### Data Preparation

- **Datasets**: Sample datasets are stored in the `/content/workshop_datasets/` directory.
- **Directory Setup**: Create directories to store data at different stages of analysis.

## Mini Analysis Vignette

This vignette demonstrates a typical single-cell RNA-seq workflow, encompassing key steps from data preparation to pseudotime analysis. Follow along with the provided notebook to gain hands-on experience.

### Submodule 1: Quality Control

- **Mitochondrial Content**: Filter low-quality cells based on mitochondrial gene content.
- **Doublet Detection**: Use Scrublet to detect and remove doublets.
- **Interactive Quality Control**: Adjust filtering thresholds using interactive sliders to refine the dataset.

### Submodule 2: Normalization and Filtering

- **Normalization**: Normalize gene expression counts to mitigate technical variation.
- **Variable Gene Selection**: Select the top 2000 highly variable genes for further analysis.

### Submodule 3: Dimensionality Reduction

- **PCA**: Perform PCA to reduce dimensionality and visualize variance.
- **UMAP**: Generate UMAP embeddings to visualize clusters and capture major data structures.

### Submodule 4: Clustering and Marker Genes

- **Louvain Clustering**: Apply Louvain clustering to identify distinct cell subpopulations.
- **Marker Genes**: Identify and visualize cluster-specific marker genes using dot plots.

### Submodule 5: Cell Type Annotation

- Annotate clusters based on known marker gene profiles to assign biological identities to cells.

### Submodule 6: Pseudotime Analysis

- **Pseudotime Inference**: Use pseudotime analysis to infer cellular differentiation trajectories and explore dynamic biological processes.

## Usage

### Running the Notebook

1. **Clone the Repository**:

   ```bash
   git clone https://github.com/your_username/your_repository.git
   cd your_repository
   ```

2. **Run the Jupyter Notebook**:

   ```bash
   jupyter notebook single_cell_workshop.ipynb
   ```

3. **Follow Along**: Each section includes explanations, code, and visualizations to guide you through the analysis.

### Example

- Load and normalize data:

  ```python
  import scanpy as sc
  adata = sc.read_h5ad("path/to/your/dataset.h5ad")
  sc.pp.normalize_total(adata, target_sum=1e4)
  sc.pp.log1p(adata)
  ```

## Additional Resources

- **Scanpy Documentation**: [https://scanpy.readthedocs.io/](https://scanpy.readthedocs.io/)
- **Single-Cell Analysis Guide**: [scanpy_usage](https://github.com/scverse/scanpy_usage)

## License

This workshop and accompanying materials are licensed under the MIT License.

## Contact

- **Name**: Ioannis (Yannis) Ntekas
- **Email**: [your_email@example.com](mailto:your_email@example.com)

# Single-Cell RNA-seq Analysis Workshop

## Overview

This workshop provides a comprehensive guide to analyzing single-cell RNA sequencing (scRNA-seq) data, covering essential steps from quality control to cell type identification. Participants will use Python libraries such as `scanpy` for data preprocessing, visualization, and analysis.

This workshop is part of **6120 Precision and Genomic Medicine**, taught by **Iwijn De Vlaminck**.

## Submodules Covered

1. **Environment Setup**: Install the necessary Python packages (`scanpy`, `louvain`, `ipywidgets`) to prepare your computational environment for scRNA-seq analysis.
2. **Data Preparation**: Load datasets, create directories, and retrieve gene annotations for proper data organization and accessibility.
3. **Quality Control**: Perform quality checks on the cells by filtering based on mitochondrial content, total counts, and gene counts. Use Scrublet to detect potential doublets and remove them.
4. **Normalization and Filtering**: Normalize the scRNA-seq data to mitigate batch effects, extract highly variable genes, and filter out sex-specific genes to refine the dataset.
5. **Dimensionality Reduction**: Apply Principal Component Analysis (PCA) to reduce data dimensionality and Uniform Manifold Approximation and Projection (UMAP) for visualizing the high-dimensional structure of the data.
6. **Clustering and Differential Expression**: Use Louvain clustering to identify distinct cell populations, and determine differentially expressed marker genes. Interactive widgets allow for adjusting clustering parameters.
7. **Cell Type Annotation**: Annotate cell clusters using known marker genes to assign cell type identities.
8. **Pseudotime Analysis**: Infer pseudotime to explore cellular differentiation trajectories, providing insight into the dynamic processes occurring within the cell populations.

## Getting Started

### Prerequisites

- **Python 3.6+**
- **Jupyter Notebook** or **Google Colab**

### Installing Dependencies

To install the required packages, run the following command:

```python
!pip install scanpy python-igraph louvain ipywidgets pybiomart
```

### Data Preparation

- **Datasets**: Sample datasets are stored in the `/content/workshop_datasets/` directory.
- **Directory Setup**: Create directories to store data at different stages of analysis.

## Mini Analysis Vignette

This vignette demonstrates a typical single-cell RNA-seq workflow, encompassing key steps from data preparation to pseudotime analysis. Follow along with the provided notebook to gain hands-on experience.

### Submodule 1: Quality Control

- **Mitochondrial Content**: Filter low-quality cells based on mitochondrial gene content.
- **Doublet Detection**: Use Scrublet to detect and remove doublets.
- **Interactive Quality Control**: Adjust filtering thresholds using interactive sliders to refine the dataset.

### Submodule 2: Normalization and Filtering

- **Normalization**: Normalize gene expression counts to mitigate technical variation.
- **Variable Gene Selection**: Select the top 2000 highly variable genes for further analysis.

### Submodule 3: Dimensionality Reduction

- **PCA**: Perform PCA to reduce dimensionality and visualize variance.
- **UMAP**: Generate UMAP embeddings to visualize clusters and capture major data structures.

### Submodule 4: Clustering and Marker Genes

- **Louvain Clustering**: Apply Louvain clustering to identify distinct cell subpopulations.
- **Marker Genes**: Identify and visualize cluster-specific marker genes using dot plots.

### Submodule 5: Cell Type Annotation

- Annotate clusters based on known marker gene profiles to assign biological identities to cells.

### Submodule 6: Pseudotime Analysis

- **Pseudotime Inference**: Use pseudotime analysis to infer cellular differentiation trajectories and explore dynamic biological processes.

## Usage

### Running the Notebook

1. **Clone the Repository**:

   ```bash
   git clone https://github.com/your_username/your_repository.git
   cd your_repository
   ```

2. **Run the Jupyter Notebook**:

   ```bash
   jupyter notebook single_cell_workshop.ipynb
   ```

3. **Follow Along**: Each section includes explanations, code, and visualizations to guide you through the analysis.

### Example

- Load and normalize data:

  ```python
  import scanpy as sc
  adata = sc.read_h5ad("path/to/your/dataset.h5ad")
  sc.pp.normalize_total(adata, target_sum=1e4)
  sc.pp.log1p(adata)
  ```

## Additional Resources

- **Scanpy Documentation**: [https://scanpy.readthedocs.io/](https://scanpy.readthedocs.io/)
- **Single-Cell Analysis Guide**: [scanpy_usage](https://github.com/scverse/scanpy_usage)

## License

This workshop and accompanying materials are licensed under the MIT License.

## Contact

- **Name**: Ioannis (Yannis) Ntekas
- **Email**: [in68@cornell.com](mailto:in68@cornell.com)


