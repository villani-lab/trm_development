Supplemental figure 1
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

mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
# mpl.rcParams['pdf.fonttype'] = 42

# Set a colormap
gene_colormap = clr.LinearSegmentedColormap.from_list('gene_cmap', ["#e0e0e1", '#4576b8', '#02024a'], N=200)

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

filtered2_no_skin2 = pg.read_input(
    os.path.join(file_path(), "data", "integrated", "filtered2_no_skin2_harmonized_with_subclust.h5ad"))
```

<img src="supp_figure_1_files/figure-gfm/unnamed-chunk-2-1.png" width="672" />

``` python

genes =  ["Cd3d", "Cd8a", "Trbc2"]

ncols = int(round(np.sqrt(len(genes))))
fig, ax = plt.subplots(ncols = ncols, nrows = ncols, figsize = (7, 7))
ax = ax.ravel()

for num, gene in enumerate(genes) :
    plot_df = pd.DataFrame(filtered2_no_skin2[:,gene].X.toarray(), columns = [gene], index = filtered2_no_skin2.obs_names)
    plot_df["x"] = filtered2_no_skin2.obsm["X_umap"][:, 0]
    plot_df["y"] = filtered2_no_skin2.obsm["X_umap"][:, 1]
    hb = ax[num].hexbin(plot_df["x"], plot_df["y"], C=plot_df[gene], cmap=gene_colormap, gridsize=250, edgecolors = "none")
    _ = ax[num].get_xaxis().set_ticks([])
    _ = ax[num].get_yaxis().set_ticks([])
    _ = ax[num].spines['top'].set_visible(False)
    _ = ax[num].spines['right'].set_visible(False)
    _ = ax[num].set_title(gene, size = 25)
    cb = fig.colorbar(hb, ax=ax[num], shrink=.75, aspect=10)
    _ = cb.ax.set_title("logCPM")
for noplot in range(num + 1, len(ax)) :
    ax[noplot].axis("off")
_ = fig.text(0.5, 0.03, 'UMAP1', va='center', size = 15)
_ = fig.text(0.03, 0.5, 'UMAP2', va='center', rotation='vertical', size = 15)
# plt.subplots_adjust(left = 0.1, bottom = 0.1)
figure = plt.gcf()
# figure.tight_layout()
figure.set_size_inches(5, 5)
figure
```

<img src="supp_figure_1_files/figure-gfm/fig_S2A-3.png" width="480" />

``` python
annot_info = pd.read_csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/supplemental_tables/cluster_annotations.csv")
annot_info["cluster_number"] = [str(n) for n in annot_info["cluster_number"]]
annot_dict = dict(zip(annot_info["cluster_number"], annot_info["annotation"]))

filtered2_no_skin2.obs["annotation"] = [annot_dict[cl] for cl in filtered2_no_skin2.obs["new_clusters"]]
# Append a "C" to the clusters as well
filtered2_no_skin2.obs["new_clusters"] = [f"C{cl}" for cl in filtered2_no_skin2.obs["new_clusters"]]

# Cleaner cluster one
sorted_clusts = sorted(set(filtered2_no_skin2.obs["annotation"]), key = lambda s: int(s.split(":")[0].replace("C", "")))
col_dict = dict(zip(sorted_clusts,filtered2_no_skin2.uns["new_clusters_colors"]))

legend_elements = [Line2D([0], [0], marker='o', color=col_dict[cl], label=cl,
                          markerfacecolor=col_dict[cl], markersize=7, lw=0) for cl in
                   sorted_clusts]

fig, ax = plt.subplots(1)
plot = sc.pl.umap(filtered2_no_skin2, color="new_clusters",
                  legend_loc="on data", show=False, ax=ax,
                  title="", legend_fontoutline=5)
lgd = ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(0.9, 0.5), frameon=False, borderpad=5)
figure = plt.gcf()
figure.set_size_inches(8, 5)
figure.tight_layout()
figure
```

<img src="supp_figure_1_files/figure-gfm/fig_S2B-5.png" width="768" />

``` python
filtered2_no_skin2.obs["tissue"] = [n if n != "LN" else "dLN" for n in filtered2_no_skin2.obs["tissue"]]
col_dict = dict(zip(sorted(set(filtered2_no_skin2.obs["tissue"])),
                    filtered2_no_skin2.uns["tissue_colors"]))

legend_elements = [Line2D([0], [0], marker='o', color=col_dict[cl], label=cl,
                          markerfacecolor=col_dict[cl], markersize=7, lw=0) for cl in
                   sorted(set(filtered2_no_skin2.obs["tissue"]))]

fig, ax = plt.subplots(1)
plot = sc.pl.umap(filtered2_no_skin2, color="tissue",
                  legend_loc="on data", show=False, ax=ax,
                  title="", legend_fontoutline=5)
lgd = ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
fig = plt.gcf()
fig.set_size_inches(4.2, 3.7)
fig.tight_layout()
fig
```

<img src="supp_figure_1_files/figure-gfm/S2C-7.png" width="403" />

``` python

embed = "umap"

fig, ax = plt.subplots()
plot_df = filtered2_no_skin2.obs[["day"]]
int_dict = {k: v for k, v in zip(sorted(set(plot_df["day"]), key=int), range(len(set(plot_df["day"]))))}
plot_df["day_int"] = [int_dict[d] for d in plot_df["day"]]
plot_df["x"] = filtered2_no_skin2.obsm[f"X_{embed}"][:, 0]
plot_df["y"] = filtered2_no_skin2.obsm[f"X_{embed}"][:, 1]
ax.hexbin(plot_df["x"], plot_df["y"], C=plot_df["day_int"], cmap=day_cmap, gridsize=500, edgecolors = "none")
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
ax.set_xlabel("{embed}1".format(embed=embed.upper()))
ax.set_ylabel("{embed}2".format(embed=embed.upper()))
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.07)
ax2 = cax.imshow(np.array([[int(min(plot_df["day"])), int(max(plot_df["day"]))]]), cmap=day_cmap,
                 interpolation="nearest", aspect=4)
cax.set_title("Day")
cb = plt.colorbar(ax2, cax=cax)
tick_locator = mpl.ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()

fig
```

<img src="supp_figure_1_files/figure-gfm/S2D-9.png" width="672" />

``` r
cell_data = read.csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/filtered2_no_skin2_harmonized_with_subclust_obs.csv",
                     row.names = 1)
cell_data$tissue[cell_data$tissue == "LN"] <- "dLN"
cell_data$tissue[cell_data$tissue == "Naive"] <- "dLN"

cell_data$new_clusters <- paste("C", cell_data$new_clusters, sep = "")
prop_info <- cell_data %>%
  dplyr::select(tissue,new_clusters, day) %>%
  dplyr::filter(tissue != "Naive") %>%
  group_by(tissue, day, new_clusters) %>%
  summarise(n_cells = n()) %>%
  group_by(day, tissue) %>%
  mutate("n_cell_tissue_day" = sum(n_cells)) %>%
  mutate("perc" = n_cells / n_cell_tissue_day)
```

    ## `summarise()` has grouped output by 'tissue', 'day'. You can override using the
    ## `.groups` argument.

``` r
prop_info$day <- factor(prop_info$day, levels = unique(prop_info$day))
prop_info$new_clusters <- factor(prop_info$new_clusters, levels = mixedsort(unique(prop_info$new_clusters)))

ggplot(prop_info, aes(y = day, x = perc, fill = new_clusters)) +
  geom_bar(stat = "identity") +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b',
                               '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a',
                               '#ff9896')) +
  facet_wrap(~tissue) +
  theme_classic(base_size = 20) +
  ylab("Day") + xlab("cell fraction") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

![](supp_figure_1_files/figure-gfm/figure_S2E-11.png)<!-- -->

``` python

# Fix the percent mito
# This isn't working with markdown...don't know why
# raw_data = pg.read_input(os.path.join(file_path(), "data", "integrated", "all_data.h5ad"))
#
# pg.qc_metrics(raw_data, percent_mito = 20, mito_prefix = "mt-")
# mito_dict = dict(zip(raw_data.obs_names, raw_data.obs["percent_mito"]))
# filtered2_no_skin2.obs['percent_mito'] = [mito_dict[c] for c in filtered2_no_skin2.obs_names]

mito_info = pd.read_csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/mito_percents.csv", index_col=0)
mito_dict = dict(zip(mito_info.index.values, mito_info["percent"]))
filtered2_no_skin2.obs['percent_mito'] = [mito_dict[c] for c in filtered2_no_skin2.obs_names]

meta = ['n_genes', 'percent_mito']
violin_dat = filtered2_no_skin2.obs[["annotation", "percent_mito", "n_genes"]]

fig, ax = plt.subplots(nrows=1, ncols=2)
ax = ax.ravel()
for num, m in enumerate(meta) :
    _ = sns.violinplot(y = "annotation", x = m, color = "grey",
                   data = violin_dat, inner = None, scale = "width",
                   ax = ax[num], cut = 0, order = sorted_clusts)
    for violin in ax[num].collections:
        violin.set_alpha(0.8)

    if m == "n_genes" :
        ax[num].set_xlabel("# genes")
    else :
        ax[num].set_xlabel("% mitochondrial UMIs")
    if num > 0 :
        ax[num].set_yticks([])
    ax[num].set_ylabel("")
figure = plt.gcf()
figure.set_size_inches(8, 3)
figure.tight_layout()
figure
```

<img src="supp_figure_1_files/figure-gfm/figure_S2F-1.png" width="768" />
