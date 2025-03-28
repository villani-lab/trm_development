Supplemental figure 6
================

``` r
library(reticulate)
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
import seaborn as sns
import matplotlib.colors as clr
from pylab import cm
import matplotlib as mpl
from matplotlib.lines import Line2D
from collections import Counter

mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
# mpl.rcParams['pdf.fonttype'] = 42

# Set a colormap
gene_colormap = clr.LinearSegmentedColormap.from_list('gene_cmap', ["#e0e0e1", '#4576b8', '#02024a'], N=200)

cmap = cm.get_cmap('YlGnBu', 140)    # PiYG
hex_list = []
for i in range(cmap.N):
    rgba = cmap(i)
    # rgb2hex accepts rgb or rgba
    hex_list.append(mpl.colors.rgb2hex(rgba))

colors = [c for n, c in enumerate(hex_list) if n%10 == 0]
colors = colors [1:13] # First one is too dim
days = ["0", "3", "4", "5", "6", "7", "10", "14", "21", "32", "60", "90"]

day_col_dict = dict(zip(days, colors))
# Also make as a colormap
kurd_day_cmap = clr.LinearSegmentedColormap.from_list('day_cmap', colors, N=len(colors))

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
kupper_day_cmap = clr.LinearSegmentedColormap.from_list('day_cmap', colors, N=len(colors))

# Set a switcher up so the script will run on any computer
def file_path(user = getpass.getuser()):
    switcher = {
            "nealp": "C:/Users/nealp/Documents/Dropbox (Partners HealthCare)/Chloe&Mazen/Collaborator_projects/Kupper_TRM/neal_analysis/",
            "neal": "/home/neal/Documents/Dropbox (Partners HealthCare)/Chloe&Mazen/Collaborator_projects/Kupper_TRM/neal_analysis/",
            "nealpsmith": "/projects/home/nealpsmith/projects/kupper/"

    }
    if switcher.get(user):
        return(switcher.get(user))
    else :
        print("Add your local filepath to the switcher! run getpass.getuser() to get your ID")

skin1_data = pg.read_input(os.path.join(file_path(), "all_data_analysis", "data", "integrated", "skin1_subcluster.h5ad"))
```

    ## 2024-08-30 17:00:26,210 - pegasus - INFO - Time spent on 'read_input' = 1.37s.

``` python
gut1_data = pg.read_input(os.path.join(file_path(), "kurd_paper", "data", "gut1_data_subcluster.h5ad"))
```

    ## 2024-08-30 17:00:29,167 - pegasus - INFO - Time spent on 'read_input' = 2.95s.

``` python

genes = ["Fos", "Fosb", "Fosl2", "Gem", "Junb", "Zfp36l1", "Nr4a2", "Nr4a1",
         "Dennd4a", "Ifrd1", "Rel", "Nr4a3", "Egr1", "Dusp1"]

ncols = int(round(np.sqrt(len(genes))))
fig, ax = plt.subplots(ncols = ncols, nrows = ncols + 1, figsize = (10, 10))
ax = ax.ravel()

for num, gene in enumerate(genes) :
    plot_df = pd.DataFrame(skin1_data[:,gene].X.toarray(), columns = [gene], index = skin1_data.obs_names)
    plot_df["x"] = skin1_data.obsm["X_umap"][:, 0]
    plot_df["y"] = skin1_data.obsm["X_umap"][:, 1]
    hb = ax[num].hexbin(plot_df["x"], plot_df["y"], C=plot_df[gene], cmap=gene_colormap, gridsize=100, edgecolors = "none")
    ax[num].get_xaxis().set_ticks([])
    ax[num].get_yaxis().set_ticks([])
    ax[num].spines['top'].set_visible(False)
    ax[num].spines['right'].set_visible(False)
    ax[num].set_title(gene, size = 25)
    cb = fig.colorbar(hb, ax=ax[num], shrink=.75, aspect=10)

for noplot in range(num + 1, len(ax)) :
    ax[noplot].axis("off")
```

    ## (0.0, 1.0, 0.0, 1.0)
    ## (0.0, 1.0, 0.0, 1.0)
    ## (0.0, 1.0, 0.0, 1.0)
    ## (0.0, 1.0, 0.0, 1.0)
    ## (0.0, 1.0, 0.0, 1.0)
    ## (0.0, 1.0, 0.0, 1.0)

``` python
fig.tight_layout()
fig
```

<img src="supp_figure_6_files/figure-gfm/fig_S7A_skin-1.png" width="960" />

``` python

fig, ax = plt.subplots(ncols = ncols, nrows = ncols + 1, figsize = (10, 10))
ax = ax.ravel()

for num, gene in enumerate(genes) :
    plot_df = pd.DataFrame(gut1_data[:,gene].X.toarray(), columns = [gene], index = gut1_data.obs_names)
    plot_df["x"] = gut1_data.obsm["X_umap"][:, 0]
    plot_df["y"] = gut1_data.obsm["X_umap"][:, 1]
    hb = ax[num].hexbin(plot_df["x"], plot_df["y"], C=plot_df[gene], cmap=gene_colormap, gridsize=100, edgecolors = "none")
    ax[num].get_xaxis().set_ticks([])
    ax[num].get_yaxis().set_ticks([])
    ax[num].spines['top'].set_visible(False)
    ax[num].spines['right'].set_visible(False)
    ax[num].set_title(gene, size = 25)
    cb = fig.colorbar(hb, ax=ax[num], shrink=.75, aspect=10)

for noplot in range(num + 1, len(ax)) :
    ax[noplot].axis("off")
```

    ## (0.0, 1.0, 0.0, 1.0)
    ## (0.0, 1.0, 0.0, 1.0)
    ## (0.0, 1.0, 0.0, 1.0)
    ## (0.0, 1.0, 0.0, 1.0)
    ## (0.0, 1.0, 0.0, 1.0)
    ## (0.0, 1.0, 0.0, 1.0)

``` python
fig.tight_layout()
fig
```

<img src="supp_figure_6_files/figure-gfm/fig_S7A_gut-3.png" width="960" />

``` python

plot_adata = pg.read_input("/projects/home/nealpsmith/projects/kupper/fitzpatrick_paper/data/filtered_data_10x.h5ad")
```

    ## 2024-08-30 17:00:48,743 - pegasus - INFO - Time spent on 'read_input' = 0.12s.

``` python
sorted_clusts = sorted(set(plot_adata.obs["leiden_labels"]))
col_dict = dict(zip(sorted_clusts,plot_adata.uns["leiden_labels_colors"]))

legend_elements = [Line2D([0], [0], marker='o', color=col_dict[cl], label=cl,
                          markerfacecolor=col_dict[cl], markersize=7, lw=0) for cl in
                   sorted_clusts]

genes = ["CD3D", "CD8A", "CD4"]
x_loc = np.min(plot_adata.obsm["X_umap"][:,0]) - np.min(plot_adata.obsm["X_umap"][:,0]) * 0.01
y_loc = np.min(plot_adata.obsm["X_umap"][:,1]) - 0.4

fig, ax = plt.subplots(ncols = 2, nrows = 2, figsize = (7.5, 5))
ax = ax.ravel()

for num, gene in enumerate(genes) :
    # Calculate n cells and percent
    n_cells = plot_adata[:, gene].X.count_nonzero()
    perc_cells = n_cells / len(plot_adata) * 100
    perc_cells = round(perc_cells, 2)

    plot_df = pd.DataFrame(plot_adata[:,gene].X.toarray(), columns = [gene], index = plot_adata.obs_names)
    plot_df["x"] = plot_adata.obsm["X_umap"][:, 0]
    plot_df["y"] = plot_adata.obsm["X_umap"][:, 1]
    hb = ax[num].hexbin(plot_df["x"], plot_df["y"], C=plot_df[gene], cmap=gene_colormap, gridsize=50, edgecolors = "none")
    ax[num].get_xaxis().set_ticks([])
    ax[num].get_yaxis().set_ticks([])
    ax[num].spines['top'].set_visible(False)
    ax[num].spines['right'].set_visible(False)
    ax[num].set_title(gene, size = 25)
    ax[num].text(x_loc, y_loc, f"{n_cells} cells ({perc_cells}%)", style="italic")
    cb = fig.colorbar(hb, ax=ax[num], shrink=.75, aspect=10)
    cb.ax.set_title("Log(CPM)")
plot = sc.pl.umap(plot_adata, color="leiden_labels",
                  legend_loc="on data", show=False, ax=ax[num + 1],
                  title="", legend_fontoutline=5)
plot.set_xlabel("")
plot.set_ylabel("")
lgd = ax[num + 1].legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

ax[num + 1].set_rasterization_zorder(2)

fig.text(0.5, 0.03, 'UMAP1', va='center', size = 15)
fig.text(0.03, 0.5, 'UMAP2', va='center', rotation='vertical', size = 15)
fig.tight_layout()
plt.subplots_adjust(left = 0.1, bottom = 0.1)
fig
```

<img src="supp_figure_6_files/figure-gfm/fig_s6c_features-5.png" width="720" />

``` python

cd8_clusts = ["C3", "C5"]

cd8_data = plot_adata[plot_adata.obs["leiden_labels"].isin(cd8_clusts)]

# Plot for overall data with CD4/CD8
genes = ["ITGAE", "ITGB2", "FOSL2", "FOSB", "JUNB", "FOS"]
x_loc = np.min(cd8_data.obsm["X_umap"][:,0]) - np.min(cd8_data.obsm["X_umap"][:,0]) * 0.01
y_loc = np.max(cd8_data.obsm["X_umap"][:,1]) - 1

fig, ax = plt.subplots(ncols = 2, nrows = 3, figsize = (8, 8))
ax = ax.ravel()

for num, gene in enumerate(genes) :
    # Calculate n cells and percent
    n_cells = cd8_data[:, gene].X.count_nonzero()
    perc_cells = n_cells / len(cd8_data) * 100
    perc_cells = round(perc_cells, 2)

    plot_df = pd.DataFrame(cd8_data[:,gene].X.toarray(), columns = [gene], index = cd8_data.obs_names)
    plot_df["x"] = cd8_data.obsm["X_umap"][:, 0]
    plot_df["y"] = cd8_data.obsm["X_umap"][:, 1]
    hb = ax[num].hexbin(plot_df["x"], plot_df["y"], C=plot_df[gene], cmap=gene_colormap, gridsize=30, edgecolors = "none")
    ax[num].get_xaxis().set_ticks([])
    ax[num].get_yaxis().set_ticks([])
    ax[num].spines['top'].set_visible(False)
    ax[num].spines['right'].set_visible(False)
    ax[num].set_title(gene, size = 25)
    ax[num].text(x_loc, y_loc, f"{n_cells} cells ({perc_cells}%)", style="italic")
    cb = fig.colorbar(hb, ax=ax[num], shrink=.75, aspect=10)
    cb.ax.set_title("Log(CPM)")
fig.text(0.5, 0.03, 'UMAP1', va='center', size = 15)
fig.text(0.03, 0.5, 'UMAP2', va='center', rotation='vertical', size = 15)
fig.tight_layout()
plt.subplots_adjust(left = 0.1, bottom = 0.1)
fig
```

<img src="supp_figure_6_files/figure-gfm/fig_s6c_ap1_features-7.png" width="768" />

``` python

fig, ax = plt.subplots(1)
plot = sc.pl.dotplot(cd8_data,  ["ITGAE", "ITGB2", "FOSL2", "FOSB", "JUNB", "FOS"], groupby="leiden_labels",
                         show=False, return_fig=True,
                         cmap="Blues", standard_scale = None, swap_axes=True, ax = ax)

axes_dict = plot.get_axes()
axes_dict["mainplot_ax"].set_axisbelow(True)
axes_dict["mainplot_ax"].grid()
figure = plt.gcf()
figure.set_size_inches(3, 5)
plt.subplots_adjust(left = 0.2, bottom = 0.1)
fig
```

<img src="supp_figure_6_files/figure-gfm/fig_s6c_dots-9.png" width="288" />
