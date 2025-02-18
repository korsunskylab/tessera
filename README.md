# Tessera

**Tessera** is an algorithm for segmenting single-cell resolution spatial omics data into small multicellular tiles whose edges track with natural tissue boundaries. These tiles can then be used in downstream analysis to label and define tissue regions across samples.

![cartoon](img/cartoon.png)

## System Requirements

### OS Requirements

`tessera` is supperted for *macOS* and *Linux*. The package has been tested on the following systems:
* macOS: Sonoma (14.6.1)
* Linux: CentOS 7 (7.9.2009)

### R Dependencies

`tessera` has been tested on R versions >= 4.3. Please consult the DESCRIPTION file for more details on required R packages, including Rcpp for R and C++ integraion.

## Installation

### Install from GitHub
Open R and run:
```R
devtools::install_github('korsunskylab/tessera')
```
This will install dependencies from CRAN.

### Manually Install Dependencies with Conda/Mamba
In the command line:
```bash
mamba create -n tessera_env
mamba activate tessera_env

## If using macOS with Apple Silicon:
# conda config --env --set subdir osx-64

# Install required dependencies (~1 min)
mamba install -c conda-forge r-essentials r-rcpp r-rcpparmadillo r-bh r-devtools \
r-tidyverse r-matrix r-rlang r-r.utils r-sf r-igraph r-furrr r-future r-data.table \
r-geometry r-mclust r-rspectra r-magrittr r-cowplot r-rhpcblasctl

# Optionally install Seurat and additional packages to run vignettes (~15 sec)
mamba install -c conda-forge r-seurat r-ggthemes r-patchwork r-viridis jupyterlab r-irkernel
```
Next, in an R console:
```R
install.packages("harmony")  # ~1 min
devtools::install_github('korsunskylab/tessera', dependencies = FALSE)  # ~1 min
```

## Quick Start

### Standalone Mode

Check out this [vignette](https://github.com/korsunskylab/tessera/blob/main/vignettes/vignette_basic.ipynb) for a quick start tutorial which
demonstrates using `tessera` in standalone mode on a single sample.

The basic usage is as follows:
```R
res = GetTiles(
    X = meta_data$X,         # Vector of cell spatial coordinates
    Y = meta_data$Y,         # Vector of cell spatial coordinates
    counts = counts,         # Gene-by-cell matrix of transcript counts

    embeddings = embeddings, # (Optional) Cell-by-embedding matrix of pre-computed cell embeddings. If missing, embeddings are calculated using PCA.
    meta_data = meta_data,   # (Optional) Additional cell meta data
    meta_vars_include = meta_vars_include, # (Optional) Cell meta data to include in output

    group.by = 'sample_id',  # (Optional) Name of meta_data column that provides sample IDs. If missing, treated as a single sample.
    ...                      # Additional Tessera algorithm parameters
)
dmt = res$dmt                # Tiling results with mapping from cells to tiles
aggs = res$aggs              # Meta data for aggregated tiles
```

### Seurat Objects (Multi-sample)
Tessera can also be applied directly to a Seurat object containing single cells with spatial coordinates.
The `GetTiles` function can use cell embeddings that have already been pre-computed (and integrated, if there are multiple samples).
By default, the output is a pair of Seurat objects: 1) a single-cell object updated with tile assignments for each cell, and 2) a Seurat object
where each entry represents an individual Tessera tile.
```R
future::plan(future::multicore)
res = GetTiles(
    obj,        # Single-cell Seurat object
    'spatial',  # Name of dimesional reduction where x/y coordinates are stored

    embeddings = 'harmony',  # (Optional) Name of dimensional reduction where pre-computed single-cell embeddings are stored
    group.by = 'sample_id',  # (Optional) Name of meta.data column that provides sample IDs. If missing, treated as a single sample.

    # Additional Tessera algorithm parameters
    prune_thresh_quantile = 0.99, prune_min_cells = 1, ...
)
obj = res$obj                # Seurat object of single-cells (with cell-to-tile mapping)
tile_obj = res$tile_obj      # Seurat object of Tessera tiles
```

## Vignettes: 
(1) Quickstart (approx. runtime: <10 sec)

https://github.com/korsunskylab/tessera/blob/main/vignettes/vignette_basic.ipynb

(2) Walkthrough (approx. runtime: <10 sec)

https://github.com/korsunskylab/tessera/blob/main/vignettes/vignette_stepthrough.ipynb
