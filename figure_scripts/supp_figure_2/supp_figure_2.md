Supplemental figure 2
================

``` r
library(reticulate)
library(gtools)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(glue)
library(magrittr)
library(VennDiagram)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(ggrepel)

use_python("/projects/home/nealpsmith/.conda/envs/old_peg_github/bin/python")
```

``` python
import getpass
import pegasus as pg
```

    ## WARNING:param.Parameterized: Use method 'warning' via param namespace 
    ## WARNING:param.main: pandas could not register all extension types imports failed with the following error: cannot import name 'ABCIndexClass' from 'pandas.core.dtypes.generic' (/projects/home/nealpsmith/.conda/envs/old_peg_github/lib/python3.7/site-packages/pandas/core/dtypes/generic.py)

``` python
import scanpy as sc
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import anndata
import math
import seaborn as sns
import matplotlib.colors as clr
from pylab import cm
import matplotlib as mpl
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.sparse import csr_matrix
from collections import Counter
import wot
import pickle

from cellrank.external.kernels import WOTKernel
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.tl.estimators import GPCCA
from mpl_toolkits.axes_grid1 import make_axes_locatable, Size
import wot

mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
# mpl.rcParams['pdf.fonttype'] = 42

# Set a colormap
gene_colormap = clr.LinearSegmentedColormap.from_list('gene_cmap', ["#e0e0e1", '#4576b8', '#02024a'], N=200)

cmap = cm.get_cmap('YlOrRd', 110)    # PiYG
hex_list = []
for i in range(cmap.N):
    rgba = cmap(i)
    # rgb2hex accepts rgb or rgba
    hex_list.append(mpl.colors.rgb2hex(rgba))

colors = [c for n, c in enumerate(hex_list) if n%10 == 0]
colors = colors [1:11] # First one is too dim
days = ["0", "2", "5", "10", "15", "20", "25", "30", "45", "60"]

day_col_dict = dict(zip(days, colors))
day_cmap = clr.LinearSegmentedColormap.from_list('day_cmap', colors, N=len(colors))

# Set a switcher up so the script will run on any computer
def file_path(user = getpass.getuser()):
    switcher = {
            "nealp": "C:/Users/nealp/Documents/Dropbox (Partners HealthCare)/Chloe&Mazen/Collaborator_projects/Kupper_TRM/neal_analysis/all_data_analysis",
            "neal": "/home/neal/Documents/Dropbox (Partners HealthCare)/Chloe&Mazen/Collaborator_projects/Kupper_TRM/neal_analysis/all_data_analysis",
            "nealpsmith": "/projects/home/nealpsmith/projects/kupper/all_data_analysis"

    }
    if switcher.get(user):
        return(switcher.get(user))
    else :
        print("Add your local filepath to the switcher! run getpass.getuser() to get your ID")

filtered2_no_skin2 = pg.read_input(
    os.path.join(file_path(), "data", "integrated", "filtered2_no_skin2_harmonized_with_subclust.h5ad"))
```

    ## 2024-08-30 15:42:41,290 - pegasus - INFO - Time spent on 'read_input' = 5.08s.

``` python

# Add fle info to .obs for convenience
filtered2_no_skin2.obs["fle_x"] = filtered2_no_skin2.obsm["X_fle"][:,0]
filtered2_no_skin2.obs["fle_y"] = filtered2_no_skin2.obsm["X_fle"][:,1]

# Want to make a panel for each day
fig, ax = plt.subplots(nrows=2, ncols=5, figsize=(16, 5))
ax = ax.ravel()
tmpts = []
for num, tmpt in enumerate(sorted(set(filtered2_no_skin2.obs["day"]), key=int)):
    tmpts = tmpts + [tmpt]

    data = filtered2_no_skin2[filtered2_no_skin2.obs["day"].isin(tmpts)]

    for days in sorted(set(data.obs["day"]), key=int):
        col = [day_col_dict[tmp] for tmp in data.obs["day"][data.obs["day"] == days]]
        x = data.obsm["X_fle"][:, 0][data.obs["day"] == days]
        y = data.obsm["X_fle"][:, 1][data.obs["day"] == days]
        ax[num].scatter(x, y, s=0.1, c=col, label=days)
    ax[num].set_xticks([])
    ax[num].set_yticks([])
    ax[num].set_title("Day {tmpt}".format(tmpt=tmpt))

fig.suptitle("FLE clustering by day", size=15, fontweight="bold")
fig.text(0.5, 0.05, 'FLE1', va='center', size=12, fontweight="bold")
fig.text(0.08, 0.5, 'FLE2', va='center', rotation='vertical', size=12, fontweight="bold")
fig
```

<img src="supp_figure_2_files/figure-gfm/fig_S3A-1.png" width="1536" />

``` python

cl_8_9_10 = filtered2_no_skin2[filtered2_no_skin2.obs["new_clusters"].isin(["8", "9", "10"])]

genes = ["Ccl5", "Ifit1", "Isg20", "Btg1", "Hcst", "Sell"]

ncols = int(round(np.sqrt(len(genes))))
fig, ax = plt.subplots(ncols = ncols, nrows = ncols + 1, figsize = (7, 10))
ax = ax.ravel()

for num, gene in enumerate(genes) :
    plot_df = pd.DataFrame(cl_8_9_10[:,gene].X.toarray(), columns = [gene], index = cl_8_9_10.obs_names)
    plot_df["x"] = cl_8_9_10.obsm["X_fle"][:, 0]
    plot_df["y"] = cl_8_9_10.obsm["X_fle"][:, 1]
    ax[num].hexbin(plot_df["x"], plot_df["y"], C=plot_df[gene], cmap=gene_colormap, gridsize=250)
    ax[num].get_xaxis().set_ticks([])
    ax[num].get_yaxis().set_ticks([])
    ax[num].spines['top'].set_visible(False)
    ax[num].spines['right'].set_visible(False)
    ax[num].set_title(gene, size = 25)
for noplot in range(num + 1, len(ax)) :
    ax[noplot].axis("off")
fig.tight_layout()
fig
```

<img src="supp_figure_2_files/figure-gfm/fig_S3C_cl8-3.png" width="672" />

``` python

genes = ["Stmn1", "Mki67", "Pclaf", "Birc5", "Ezh2", "Cdk1"]

ncols = int(round(np.sqrt(len(genes))))
fig, ax = plt.subplots(ncols = ncols, nrows = ncols + 1, figsize = (7, 10))
ax = ax.ravel()

for num, gene in enumerate(genes) :
    plot_df = pd.DataFrame(cl_8_9_10[:,gene].X.toarray(), columns = [gene], index = cl_8_9_10.obs_names)
    plot_df["x"] = cl_8_9_10.obsm["X_fle"][:, 0]
    plot_df["y"] = cl_8_9_10.obsm["X_fle"][:, 1]
    ax[num].hexbin(plot_df["x"], plot_df["y"], C=plot_df[gene], cmap=gene_colormap, gridsize=250)
    ax[num].get_xaxis().set_ticks([])
    ax[num].get_yaxis().set_ticks([])
    ax[num].spines['top'].set_visible(False)
    ax[num].spines['right'].set_visible(False)
    ax[num].set_title(gene, size = 25)
for noplot in range(num + 1, len(ax)) :
    ax[noplot].axis("off")
fig.tight_layout()
fig
```

<img src="supp_figure_2_files/figure-gfm/fig_S3C_cl9-5.png" width="672" />

``` python

genes = ["Batf3", "Fabp5", "Fut7", "Il2ra", "Cxcr4", "Havcr2"]

ncols = int(round(np.sqrt(len(genes))))
fig, ax = plt.subplots(ncols = ncols, nrows = ncols + 1, figsize = (7, 10))
ax = ax.ravel()

for num, gene in enumerate(genes) :
    plot_df = pd.DataFrame(cl_8_9_10[:,gene].X.toarray(), columns = [gene], index = cl_8_9_10.obs_names)
    plot_df["x"] = cl_8_9_10.obsm["X_fle"][:, 0]
    plot_df["y"] = cl_8_9_10.obsm["X_fle"][:, 1]
    ax[num].hexbin(plot_df["x"], plot_df["y"], C=plot_df[gene], cmap=gene_colormap, gridsize=250)
    ax[num].get_xaxis().set_ticks([])
    ax[num].get_yaxis().set_ticks([])
    ax[num].spines['top'].set_visible(False)
    ax[num].spines['right'].set_visible(False)
    ax[num].set_title(gene, size = 25)
for noplot in range(num + 1, len(ax)) :
    ax[noplot].axis("off")
fig.tight_layout()
fig
```

<img src="supp_figure_2_files/figure-gfm/fig_S3C_cl10-7.png" width="672" />
