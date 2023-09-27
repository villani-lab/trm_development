Figure 1
================

``` r
library(reticulate)
library(gtools)
library(tidyverse)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(magrittr)
library(glue)
library(gtools)

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
# import harmonypy as hm
# import neals_python_functions as nealsucks
from mpl_toolkits.axes_grid1 import make_axes_locatable
# import rpy2.robjects as ro
# from rpy2.robjects import pandas2ri
# from rpy2.robjects.conversion import localconverter
from scipy.sparse import csr_matrix
from collections import Counter
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

    ## 2023-09-27 12:49:44,487 - pegasus - INFO - Time spent on 'read_input' = 5.37s.

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
plot = sc.pl.embedding(filtered2_no_skin2, basis = "fle", color="new_clusters",
                  legend_loc="on data", show=False, ax=ax,
                  title="", legend_fontoutline=5)
lgd = ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(0.9, 0.5), frameon=False, borderpad=5)
figure = plt.gcf()
figure.set_size_inches(8, 5)
figure.tight_layout()
figure
```

<img src="figure_1_files/figure-gfm/fig_1B-1.png" width="768" />

``` python
filtered2_no_skin2.obs["tissue"] = [n if n != "LN" else "dLN" for n in filtered2_no_skin2.obs["tissue"]]
col_dict = dict(zip(sorted(set(filtered2_no_skin2.obs["tissue"])),
                    filtered2_no_skin2.uns["tissue_colors"]))

legend_elements = [Line2D([0], [0], marker='o', color=col_dict[cl], label=cl,
                          markerfacecolor=col_dict[cl], markersize=7, lw=0) for cl in
                   sorted(set(filtered2_no_skin2.obs["tissue"]))]

fig, ax = plt.subplots(1)
plot = sc.pl.embedding(filtered2_no_skin2, basis = "fle", color="tissue",
                  legend_loc="on data", show=False, ax=ax,
                  title="", legend_fontoutline=5)
lgd = ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
fig = plt.gcf()
fig.set_size_inches(4.2, 3.7)
fig.tight_layout()
fig
```

<img src="figure_1_files/figure-gfm/1C_umap-3.png" width="403" />

``` python

embed = "fle"

fig, ax = plt.subplots()
plot_df = filtered2_no_skin2.obs[["day"]]
int_dict = {k: v for k, v in zip(sorted(set(plot_df["day"]), key=int), range(len(set(plot_df["day"]))))}
plot_df["day_int"] = [int_dict[d] for d in plot_df["day"]]
```

    ## /projects/home/nealpsmith/.conda/envs/old_peg_github/bin/python:1: SettingWithCopyWarning: 
    ## A value is trying to be set on a copy of a slice from a DataFrame.
    ## Try using .loc[row_indexer,col_indexer] = value instead
    ## 
    ## See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy

``` python
plot_df["x"] = filtered2_no_skin2.obsm[f"X_{embed}"][:, 0]
plot_df["y"] = filtered2_no_skin2.obsm[f"X_{embed}"][:, 1]
ax.hexbin(plot_df["x"], plot_df["y"], C=plot_df["day_int"], cmap=day_cmap, gridsize=500, edgecolors = "none")
ax.get_xaxis().set_ticks([])
```

    ## []

``` python
ax.get_yaxis().set_ticks([])
```

    ## []

``` python
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

<img src="figure_1_files/figure-gfm/1D_umap-5.png" width="672" />

``` python

clust_perc_by_day = pd.DataFrame(filtered2_no_skin2.obs.groupby(["new_clusters", "day"])["day"].count()) \
        .rename(columns={"day": "n_cells"}) \
        .reset_index()
clust_perc_by_day["perc"] = clust_perc_by_day.groupby(["new_clusters"])["n_cells"].transform(lambda x: x / x.sum())

clust_perc_by_tissue = pd.DataFrame(filtered2_no_skin2.obs.groupby(["new_clusters", "tissue"])["tissue"].count()) \
        .rename(columns={"tissue": "n_cells"}) \
        .reset_index()
clust_perc_by_tissue["perc"] = clust_perc_by_tissue.groupby(["new_clusters"])["n_cells"].transform(
    lambda x: x / x.sum())
```

``` r
prop_info <- reticulate::py$clust_perc_by_day
prop_info$day <- factor(prop_info$day, levels = unique(prop_info$day))
prop_info$new_clusters <- factor(prop_info$new_clusters, levels = mixedsort(unique(prop_info$new_clusters)))

ggplot(prop_info, aes(y = new_clusters, x = perc, fill = day)) +
  geom_bar(stat = "identity") +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = c('#fff2ac', '#ffe48c', '#fed16e', '#feb54f', '#fd9941', '#fd7435',
                               '#f94728', '#e6211e', '#cc0a22', '#a80026')) +
  theme_classic(base_size = 20) +
  ylab("") + xlab("cell fraction") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

![](figure_1_files/figure-gfm/1C_D_bars-7.png)<!-- -->

``` r
prop_info <- reticulate::py$clust_perc_by_tissue
colnames(prop_info)[colnames(prop_info) == "tissue"] <- "source"
prop_info$source[prop_info$source == "LN"] <- "dLN"
prop_info$source <- factor(prop_info$source, levels = c("Naive", "Skin", "dLN"))
prop_info$new_clusters <- factor(prop_info$new_clusters, levels = mixedsort(unique(prop_info$new_clusters)))

ggplot(prop_info, aes(y = new_clusters, x = perc, fill = source)) +
  geom_bar(stat = "identity") +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = c('#1f77b4', '#ff7f0e', '#2ca02c')) +
  theme_classic(base_size = 20) +
  ylab("") + xlab("cell fraction") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

![](figure_1_files/figure-gfm/1C_D_bars-8.png)<!-- -->

``` r
count_mtx <- read.csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/filtered2_no_skin2_harmonized_new_clusters_pseudobulk_on_cluster_only_counts.csv",
                      row.names = 1)
meta_data <- read.csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/filtered2_no_skin2_harmonized_new_clusters_pseudobulk_on_cluster_only_meta.csv",
                      row.names = 1)

annotation_info <- read.csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/heatmap_info/filtered2_no_skin2_annotation_genes.csv", row.names = 1)
colnames(annotation_info) <- sapply(str_replace(colnames(annotation_info), "_", " "), str_to_title)
norm_counts <- apply(count_mtx, 2, function(c){
  n_total <- sum(c)
  per_100k <- (c * 1000000) / n_total
  return(per_100k)
})
norm_counts <- log1p(norm_counts)

all_genes <- apply(annotation_info, 2, function(c){
  genes <- c[0:(nrow(annotation_info) - 3)]
  genes <- genes[genes != ""]
  return(genes)
}) %>%
  unlist(use.names = FALSE)
heatmap_data <- norm_counts[unique(all_genes),]
heatmap_data <- t(scale(t(heatmap_data)))
colnames(heatmap_data) <- sapply(colnames(heatmap_data), function(x) paste("Cluster", strsplit(x, "_")[[1]][2]))

clustering = hclust(dist(t(heatmap_data), method = "euclidean"), method = "ward.D2")
col_hc <- as.dendrogram(clustering)

order <- clustering$labels[clustering$order]

gene_order <- lapply(order, function(x){
  genes <- annotation_info[[x]][1:(nrow(annotation_info) - 3)]
  genes <- genes[genes != ""]
  return(genes)
}) %>%
  unlist(use.names = FALSE)

heatmap_data = heatmap_data[gene_order,]


label_genes <- lapply(colnames(annotation_info), function(c){
  genes <- as.character(unlist(annotation_info[1:4,c]))
  return(genes)
})
names(label_genes) <- colnames(annotation_info)

# Some custom labelings for interesting genes
label_genes$`Cluster 6` <- c("Nfkbid", "Cd69", "Klf6", "Tnf")
label_genes$`Cluster 3` <- c("Nabp1", "Foxo1", "Il21r", "Cwc25")
label_genes$`Cluster 2` <- c("Lpxn", "Ifngr1", "Stat4", "Tnfrsf18")

annotation_genes <- unlist(label_genes, use.names = FALSE)

# Now lets organize the color info that will be used for annotations
col_info = annotation_info %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(-mean_genes) %>%
  rownames_to_column(var = "cluster") %>%
  reshape2::melt(id.vars = c("cluster", "col")) %>%
  select(-variable)


# Get the gene colors
gene_cols = c()
for (gene in annotation_genes){
  color = as.character(filter(col_info, value == gene)["col"][[1]])
  gene_cols = c(gene_cols, color)
}


mean_genes <- annotation_info["mean_genes",] %>%
  mutate_each(funs(as.numeric(as.character(.)))) %>%
  select(colnames(heatmap_data)) # To order them like they will be ordered in the heatmap (same as how GEX data was read in)

gene_col_fun <- colorRamp2(c(min(mean_genes), max(mean_genes)), c("#1d111d", "#bbe7c8"))
gene_bar <-  HeatmapAnnotation("mean # genes" = as.numeric(mean_genes), col = list("mean # genes" = gene_col_fun), show_legend = FALSE)
gene_lgd <- Legend(col_fun = gene_col_fun, title = "# genes", legend_height = unit(4, "cm"), title_position = "topcenter")


heatmap_col_fun = colorRamp2(c(min(heatmap_data), 0, max(heatmap_data)), c("purple", "black", "yellow"))
heatmap_lgd = Legend(col_fun = heatmap_col_fun, title = "z-score", legend_height = unit(4, "cm"), title_position = "topcenter")

lgd_list <- packLegend(heatmap_lgd, gene_lgd, column_gap = unit(1,"cm"), direction = "horizontal")


split <- c()
for (gene in rownames(heatmap_data)){
  for (cl in order){
    if(gene %in% annotation_info[[cl]]){
      split <- c(split, sub("Cluster", "", cl))
      break
    }
  }
}
split <- factor(split, levels = as.character(unique(split)))


# Get the cluster colors
col_label_colors <- c()
for (clust in colnames(heatmap_data)){
  color <- col_info %>%
    select(cluster, col) %>%
    distinct() %>%
    filter(cluster == clust)
  col_label_colors <- c(col_label_colors, as.character(color$col))
}

# Make block annotation
clust_cols <- c()
for (clust in order){
  color <- col_info %>%
    select(cluster, col) %>%
    distinct() %>%
    filter(cluster == clust)
  clust_cols <- c(clust_cols, as.character(color$col))
}

left_annotation =   HeatmapAnnotation(blk = anno_block(gp = gpar(fill = clust_cols, col = clust_cols)),
                                      which = "row", width = unit(1.5, "mm"))
heatmap_list = Heatmap(heatmap_data, name = "z-score", col = heatmap_col_fun, cluster_rows = FALSE, cluster_columns = col_hc,
                       clustering_method_columns = "ward.D2", clustering_distance_columns = "euclidean",
                       column_dend_reorder = FALSE, top_annotation = gene_bar, show_heatmap_legend = FALSE,
                       column_names_gp = gpar(col = col_label_colors, fontface = "bold"),
                       split = split, left_annotation = left_annotation, show_column_names = TRUE) +
  rowAnnotation(link = anno_mark(at = match(annotation_genes, rownames(heatmap_data)),labels = annotation_genes,
                                 labels_gp = gpar(col = gene_cols, fontsize = 12, fontface = "bold")))
draw(heatmap_list, heatmap_legend_list =lgd_list, padding = unit(c(0.5, 0.5, 2, 2), "cm"))
```

![](figure_1_files/figure-gfm/fig_1E-1.png)<!-- -->

``` python

markers = ["Sell", "Tcf7", "Tnfrsf18", "Gzmb", "Ctla2a", "Itgae", "Junb", "Ezh2", "Eif3b", "Tnf",
           "Cd69", "Ccl5", "Il7r", "Ifit3", "Isg15", "Mki67", "Stmn1", "Fut7", "Fabp5", "H2-Ab1",
           "H2-Eb1", "Ccr7"]

plot = sc.pl.dotplot(filtered2_no_skin2, markers, groupby="new_clusters",
                     show=False, return_fig=True, title="cannonical markers",
                     cmap = "Blues", standard_scale = "var")
axes_dict = plot.get_axes()
axes_dict["mainplot_ax"].set_axisbelow(True)
axes_dict["mainplot_ax"].grid()
fig = plt.gcf()
fig.set_size_inches(6, 5)
plt.subplots_adjust(left = 0.1, bottom = 0.2)
fig.tight_layout()
```

    ## /projects/home/nealpsmith/.conda/envs/old_peg_github/bin/python:1: UserWarning: This figure includes Axes that are not compatible with tight_layout, so results might be incorrect.

``` python
fig
```

<img src="figure_1_files/figure-gfm/1F-1.png" width="576" />
