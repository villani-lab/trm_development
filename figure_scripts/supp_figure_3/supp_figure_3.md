Supplemental figure 3
================

``` r
library(reticulate)
library(gtools)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(glue)
library(magrittr)
library(fgsea)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)

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

    ## 2024-08-30 15:52:06,555 - pegasus - INFO - Time spent on 'read_input' = 5.01s.

``` r
# Load in the counts and metadata
count_data_kupper <- read.csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/skin1_subcluster_pseudobulk_on_time_counts.csv",
                       row.names = 1)
meta_data_kupper <- read.csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/skin1_subcluster_pseudobulk_on_time_meta.csv",
                    row.names = 1)

# Remove the duplacte day 30 sample
count_data_kupper <- count_data_kupper %>%
  dplyr::select(-samp_19_D30_Skin_30)
colnames(count_data_kupper) <- sapply(colnames(count_data_kupper), function(x) strsplit(x, "_")[[1]][3])

# Calculate variance across genes
norm_res_kupper <- apply(count_data_kupper, 2, function(c){
  n_total <- sum(c)
  per_100k <- (c * 1000000) / n_total
  return(per_100k)
})
norm_res_kupper <- log1p(norm_res_kupper)

# calculate the variance
vars <- apply(norm_res_kupper, 1, var)
set_v <- sort(vars, decreasing = TRUE)[1:500]
hist(vars, breaks = 1000)
abline(v = min(set_v), col = "red", lwd = 2, lty = 2)
```

![](supp_figure_3_files/figure-gfm/fig_s3a_skin_hmap-1.png)<!-- -->

``` r
top_var_norm <- norm_res_kupper[names(set_v),]
top_var_norm <- t(scale(t(top_var_norm)))

clustering = hclust(dist(t(top_var_norm), method = "euclidean"), method = "ward.D2")
col_hc <- as.dendrogram(clustering)
col_hc[[1]] <- rev(col_hc[[1]])
col_hc[[2]][[1]][[2]] <- rev(col_hc[[2]][[1]][[2]])
col_hc[[2]][[1]][[2]][[1]] <- rev(col_hc[[2]][[1]][[2]][[1]])
col_hc[[2]][[2]] <- rev(col_hc[[2]][[2]])

row_clustering <- hclust(dist(top_var_norm, method = "euclidean"))
row_hc <- as.dendrogram(row_clustering)

# col_hc[[1]][[1]] <- rev(col_hc[[1]])
# temp <- col_hc[[2]][[1]]
# col_hc[[2]][[1]] <- col_hc[[2]][[2]]
# Day bar
days <- colnames(top_var_norm)

# Now make the day annotation
day_colors = list("Day" = c('D0'= '#fff2ac',
                            'D2'= '#ffe48c',
                            'D5'= '#fed16e',
                            'D10'= '#feb54f',
                            'D15'= '#fd9941',
                            'D20'= '#fd7435',
                            'D25'= '#f94728',
                            'D30'= '#e6211e',
                            'D45'= '#cc0a22',
                            'D60'= '#a80026'))

days <- factor(days, levels = c("D0", "D2", "D5", "D10", "D15", "D20", "D25", "D30", "D45", "D60"))
days <- droplevels(days)
day_bar = HeatmapAnnotation("Day" = days, col = day_colors,
                            show_legend = FALSE)
heatmap_col_fun = colorRamp2(c(min(top_var_norm), 0, max(top_var_norm)), c("purple", "black", "yellow"))

day_fill <- sapply(levels(days), function(x) day_colors$Day[as.character(x)])
day_legend <- Legend(labels = levels(days), legend_gp = gpar(fill = day_fill), title = "timepoint",
                     title_position = "topcenter",
                     nr = 3)
heatmap_lgd = Legend(col_fun = heatmap_col_fun, title = "z-score", legend_height = unit(4, "cm"),
                     title_position = "topcenter", direction = "horizontal")
lgd_list <- packLegend(heatmap_lgd, day_legend, column_gap = unit(1,"cm"), direction = "horizontal",
                       max_height = unit(14, "cm"))

hmap = Heatmap(top_var_norm, name = "z-score", col = heatmap_col_fun,
               top_annotation = day_bar, show_column_names = FALSE,
               show_row_names = FALSE, cluster_rows = row_hc,
               cluster_columns = col_hc, show_heatmap_legend = FALSE, row_split = 3,
               heatmap_legend_param = list(legend_direction = "horizontal"))
draw(hmap, heatmap_legend_list = lgd_list, heatmap_legend_side = "bottom")
```

![](supp_figure_3_files/figure-gfm/fig_s3a_skin_hmap-2.png)<!-- -->

``` r
# What if we cut the tree to get all the patterns
gene_groups <- cutree(row_clustering, 3) %>%
  as.data.frame() %>%
  `colnames<-`(c("group"))

plot_list <- list()
for (g in unique(gene_groups$group)){
  genes <- gene_groups %>%
    dplyr::filter(group == g) %>%
    rownames(.)

  plot_data <- reshape2::melt(top_var_norm) %>%
    `colnames<-`(c("gene", "day", "z_score")) %>%
    dplyr::filter(gene %in% genes) %>%
    mutate(day_int =  as.numeric(sub("D", "", .$day)))
  plot_data$day_int <- factor(plot_data$day_int)
  n_genes <- length(genes)
  p <- ggplot(plot_data, aes(x = day_int, y = z_score)) +
    geom_line(color = "grey", aes(group = gene)) +
    geom_point(pch = 21, fill = "grey", alpha = 0.5) +
    # geom_smooth(color = "red",method = "loess") +
    ggtitle(glue("group {g}; # genes  = {n_genes}")) +
    xlab("day") +
    ylab("z-score") +
    theme_classic(base_size = 20) +
    theme(plot.title = element_text(size = 15))

  plot_list <- c(plot_list, list(p))
}
ggarrange(plotlist = plot_list)
```

![](supp_figure_3_files/figure-gfm/fig_s3a_skin_lines-1.png)<!-- -->

``` r
# Wonder if we see the same thing in the gut dataset
count_data_kurd <- read.csv("/projects/home/nealpsmith/projects/kupper/kurd_paper/data/gut1_data_subcluster_pseudobulk_on_time_counts.csv",
                         row.names = 1)
meta_data_kurd <- read.csv("/projects/home/nealpsmith/projects/kupper/kurd_paper/data/gut1_data_subcluster_pseudobulk_on_time_meta.csv",
                    row.names = 1)
norm_res_kurd <- apply(count_data_kurd, 2, function(c){
  n_total <- sum(c)
  per_100k <- (c * 1000000) / n_total
  return(per_100k)
})
norm_res_kurd <- log1p(norm_res_kurd)

colname_to_tmpt <- lapply(colnames(norm_res_kurd), function(x) paste("D", tail(strsplit(x, "_")[[1]], n = 1), sep = ""))
tmpts <- unique(sapply(colnames(norm_res_kurd), function(x) tail(strsplit(x, "_")[[1]], n = 1)))

norm_res_kurd <- lapply(tmpts, function(tpt){
  cols = as.data.frame(norm_res_kurd[,grepl(glue("*_{tpt}$"), colnames(norm_res_kurd))])
  if(ncol(cols) > 1){
   mean = rowMeans(cols) %>%
    as.data.frame() %>%
    `colnames<-`(glue("D{tpt}")) %>%
     rownames_to_column("gene")
  return(mean)
  } else {
    cols <- cols %>%
      `colnames<-`(glue("D{tpt}")) %>%
     rownames_to_column("gene")
    return(cols)
  }
}) %>%
  purrr::reduce(left_join, by = "gene") %>%
  column_to_rownames("gene")

vars <- apply(norm_res_kurd, 1, var)
set_v <- sort(vars, decreasing = TRUE)[1:500]
hist(vars, breaks = 1000)
abline(v = min(set_v), col = "red", lwd = 2, lty = 2)
```

![](supp_figure_3_files/figure-gfm/fig_s3a_gut_hmap-1.png)<!-- -->

``` r
top_var_norm <- norm_res_kurd[names(set_v),]
top_var_norm <- t(scale(t(top_var_norm)))

clustering = hclust(dist(t(top_var_norm), method = "euclidean"), method = "ward.D2")
col_hc <- as.dendrogram(clustering)
col_hc[[2]] <- rev(col_hc[[2]])
col_hc[[2]][[2]][[1]] <- rev(col_hc[[2]][[2]][[1]])
col_hc[[2]][[1]] <- rev(col_hc[[2]][[1]])

row_clustering <- hclust(dist(top_var_norm, method = "euclidean"))
row_hc <- as.dendrogram(row_clustering)

# Day bar
days <- colnames(top_var_norm)

# Now make the day annotation
day_colors = list("Day" = c('D0'= '#f5fbc2',
                            'D3'= '#e7f6b1',
                            'D4'= '#d1edb3',
                            'D5'= '#b1e1b6',
                            'D6'= '#88d0ba',
                            'D7'= '#63c3bf',
                            'D10'= '#40b5c4',
                            'D14'= '#2ba0c2',
                            'D21'= '#1e88bc',
                            'D32'= '#216aae',
                            'D60'= '#2350a1',
                            'D90'= '#253896'))

days <- factor(days, levels = c("D0", "D3", "D4", "D5", "D6", "D7", "D10", "D14", "D21", "D32", "D60", "D90"))
days <- droplevels(days)
day_bar = HeatmapAnnotation("Day" = days, col = day_colors,
                            show_legend = FALSE)
heatmap_col_fun = colorRamp2(c(min(top_var_norm), 0, max(top_var_norm)), c("purple", "black", "yellow"))

day_fill <- sapply(levels(days), function(x) day_colors$Day[as.character(x)])
day_legend <- Legend(labels = levels(days), legend_gp = gpar(fill = day_fill), title = "timepoint",
                     title_position = "topcenter",
                     nr = 3)
heatmap_lgd = Legend(col_fun = heatmap_col_fun, title = "z-score", legend_height = unit(4, "cm"),
                     title_position = "topcenter", direction = "horizontal")
lgd_list <- packLegend(heatmap_lgd, day_legend, column_gap = unit(1,"cm"), direction = "horizontal",
                       max_height = unit(14, "cm"))

hmap = Heatmap(top_var_norm, name = "z-score", col = heatmap_col_fun,
               top_annotation = day_bar, show_column_names = FALSE,
               show_row_names = FALSE, cluster_rows = row_hc,
               cluster_columns = col_hc, show_heatmap_legend = FALSE, row_split = 3)

draw(hmap, heatmap_legend_list = lgd_list, heatmap_legend_side = "bottom")
```

![](supp_figure_3_files/figure-gfm/fig_s3a_gut_hmap-2.png)<!-- -->

``` r
# What if we cut the tree to get all the patterns
gene_groups <- cutree(row_clustering, 3) %>%
  as.data.frame() %>%
  `colnames<-`(c("group"))

plot_list <- list()
for (g in unique(gene_groups$group)){
  genes <- gene_groups %>%
    dplyr::filter(group == g) %>%
    rownames(.)

  plot_data <- reshape2::melt(top_var_norm) %>%
    `colnames<-`(c("gene", "day", "z_score")) %>%
    dplyr::filter(gene %in% genes) %>%
    mutate(day_int =  as.numeric(sub("D", "", .$day)))
  plot_data$day_int <- factor(plot_data$day_int)
  n_genes <- length(genes)
  p <- ggplot(plot_data, aes(x = day_int, y = z_score)) +
    geom_line(color = "grey", aes(group = gene)) +
    geom_point(pch = 21, fill = "grey", alpha = 0.5) +
    # geom_smooth(color = "red",method = "loess") +
    ggtitle(glue("group {g}; # genes  = {n_genes}")) +
    xlab("day") +
    ylab("z-score") +
    theme_classic(base_size = 20) +
    theme(plot.title = element_text(size = 15))

  plot_list <- c(plot_list, list(p))
}
ggarrange(plotlist = plot_list)
```

![](supp_figure_3_files/figure-gfm/fig_s3a_gut_lines-1.png)<!-- -->

``` r
# Read in the linear modeling results, scrip for this is in fig 3 script
res_skin <- read.csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/pseudobulk/lin_modelling_subclusters/skin1_by_day_all_genes.csv", row.names = 1)


mouse_to_human <- read.csv("/projects/home/nealpsmith/data/useful/human_to_mouse_genes.csv")
# Make the gene name column the same
colnames(res_skin)[colnames(res_skin) == "gene"] <- "MGI.symbol"

res2_skin <- res_skin %>%
  dplyr::select(MGI.symbol, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(MGI.symbol) %>%
  summarize(stat=mean(stat)) %>%
  left_join(mouse_to_human, by = "MGI.symbol") %>%
  dplyr::select(HGNC.symbol, stat) %>%
  na.omit()

skin_ranks <- deframe(res2_skin)

gene_sets <- gmtPathways("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/msigdb_symbols.gmt")

set <- "KEGG|HALLMARK"

gene_set_group = gene_sets[grep(set, names(gene_sets))]

skin_gsea <- fgsea(pathways = gene_set_group, stats=skin_ranks, nperm=10000) %>%
   dplyr::select(pathway, pval, padj, NES) %>%
  `colnames<-`(c("pathway", "skin_pval", "skin_padj", "skin_NES"))
```

    ## Warning in fgsea(pathways = gene_set_group, stats = skin_ranks, nperm = 10000):
    ## You are trying to run fgseaSimple. It is recommended to use fgseaMultilevel. To
    ## run fgseaMultilevel, you need to remove the nperm argument in the fgsea function
    ## call.

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (3.2% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize,
    ## gseaParam, : There are duplicate gene names, fgsea may produce unexpected
    ## results.

``` r
## Now the gut ###
res_gut <- read.csv("/projects/home/nealpsmith/projects/kupper/kurd_paper/data/pseudobulk/lin_modelling_subclusters/gut1_by_day_all_genes.csv", row.names = 1)
colnames(res_gut)[colnames(res_gut) == "gene"] <- "MGI.symbol"

res2_gut <- res_gut %>%
  dplyr::select(MGI.symbol, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(MGI.symbol) %>%
  summarize(stat=mean(stat)) %>%
  left_join(mouse_to_human, by = "MGI.symbol") %>%
  dplyr::select(HGNC.symbol, stat) %>%
  na.omit()

gut_ranks <- deframe(res2_gut)

set <- "KEGG|HALLMARK"

gene_set_group = gene_sets[grep(set, names(gene_sets))]

gut_gsea <- fgsea(pathways = gene_set_group, stats=gut_ranks, nperm=10000) %>%
  dplyr::select(pathway, pval, padj, NES) %>%
  `colnames<-`(c("pathway", "gut_pval", "gut_padj", "gut_NES"))
```

    ## Warning in fgsea(pathways = gene_set_group, stats = gut_ranks, nperm = 10000):
    ## You are trying to run fgseaSimple. It is recommended to use fgseaMultilevel. To
    ## run fgseaMultilevel, you need to remove the nperm argument in the fgsea function
    ## call.

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (3.91% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize,
    ## gseaParam, : There are duplicate gene names, fgsea may produce unexpected
    ## results.

``` r
all_gsea_info <- skin_gsea %>%
  dplyr::left_join(gut_gsea, by = "pathway")

shared_up <- all_gsea_info %>%
  dplyr::filter(skin_padj < 0.05, gut_padj < 0.05, skin_NES > 1, gut_NES > 1) %>%
  .$pathway
shared_down <- all_gsea_info %>%
  dplyr::filter(skin_padj < 0.05, gut_padj < 0.05, skin_NES < -1, gut_NES < -1) %>%
  .$pathway
shared_all <- c(shared_up, shared_down)


skin_res <- skin_gsea %>%
  as_tibble() %>%
  dplyr::filter(!pathway %in% shared_all) %>%
  arrange(desc(skin_NES))
n_plot <- nrow(skin_res[skin_res$skin_padj < 0.05 & abs(skin_res$skin_NES) > 1,])
ggplot(skin_res[skin_res$skin_padj < 0.05 & abs(skin_res$skin_NES) > 1,],
       aes(reorder(pathway, skin_NES), skin_NES)) +
  coord_flip() +
  geom_col() +
  labs(x="", y="Normalized Enrichment Score",
       title= "Skin")+
  theme_classic(base_size = 15)
```

![](supp_figure_3_files/figure-gfm/fig_S3B_skin-1.png)<!-- -->

``` r
gut_res <- gut_gsea %>%
  as_tibble() %>%
  dplyr::filter(!pathway %in% shared_all) %>%
  arrange(desc(gut_NES))
n_plot <- nrow(gut_res[gut_res$gut_padj < 0.05 & abs(gut_res$gut_NES) > 1,])

ggplot(gut_res[gut_res$gut_padj < 0.05 & abs(gut_res$gut_NES) > 1,],
       aes(reorder(pathway, gut_NES), gut_NES)) +
  coord_flip() +
  geom_col() +
  labs(x="", y="Normalized Enrichment Score",
       title= "Gut")+
  theme_classic(base_size = 15)
```

![](supp_figure_3_files/figure-gfm/fig_S3B_gut-1.png)<!-- -->

``` r
test_pathway <- "HALLMARK_APOPTOSIS"
pathway_to_test <- gene_sets[test_pathway]
rank_list <- list("skin" = skin_ranks, "gut" = gut_ranks)
plot_list <- list()
gsea_res <- data.frame()
for (i in 1:length(rank_list)){

  ranks <- rank_list[[i]]

  fgsea_res <- fgsea(pathways = pathway_to_test, stats = ranks, nperm = 10000)
  fgsea_res$tissue <- names(rank_list[i])
  gsea_res <- rbind(gsea_res, fgsea_res)
  annot <- names(rank_list[i])

  nes <- round(fgsea_res$NES[fgsea_res$pathway == test_pathway], 3)
  pval <- round(fgsea_res$pval[fgsea_res$pathway == test_pathway], 3)
  n_genes <- fgsea_res$size[fgsea_res$pathway == test_pathway]

  rnk <- rank(-ranks)
  ord <- order(rnk)

  statsAdj <- ranks[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ 1)
  statsAdj <- statsAdj / max(abs(statsAdj))

  pathway <- unname(as.vector(na.omit(match(pathway_to_test[[test_pathway]], names(statsAdj)))))
  pathway <- sort(pathway)

  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                            returnAllExtremes = TRUE)

  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops

  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))

  diff <- (max(tops) - min(bottoms)) / 8

  x=y=NULL

  p <- ggplot(toPlot, aes(x = x, y = y)) +
    # geom_point(color="blue", size=0.1) +
    geom_line(color="blue") +
    geom_hline(yintercept=0, colour="black") +
    geom_segment(data=data.frame(x=pathway),
                     mapping=aes(x=x, y=-0.15,
                                 xend=x, yend=-0.25),
                     size=0.4) +
    scale_y_continuous(expand = c(0.05,0.05)) +
    scale_x_continuous(breaks = seq(0, 9000, 3000)) +
    xlab("Rank") + ylab("Enrichment score") +
    geom_text(aes(label = "")) +
    annotate("text", label = glue("NES : {nes}"), x = length(ranks) - 2000, y  =0.9) +
    annotate("text", label = glue("p-value : {pval}"), x = length(ranks) - 2000, y = 0.8) +
    annotate("text", label = glue("# genes : {n_genes}"), x = length(ranks) - 2000, y = 0.7) +
    ggtitle(glue("{test_pathway} : {annot}")) +
    theme_classic(base_size = 15) +
    theme(plot.title = element_text(face = "bold", size = 10))
  plot_list <- c(plot_list, list(p))
}
```

    ## Warning in fgsea(pathways = pathway_to_test, stats = ranks, nperm = 10000): You
    ## are trying to run fgseaSimple. It is recommended to use fgseaMultilevel. To run
    ## fgseaMultilevel, you need to remove the nperm argument in the fgsea function
    ## call.

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (3.2% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize,
    ## gseaParam, : There are duplicate gene names, fgsea may produce unexpected
    ## results.

    ## Warning in fgsea(pathways = pathway_to_test, stats = ranks, nperm = 10000): You
    ## are trying to run fgseaSimple. It is recommended to use fgseaMultilevel. To run
    ## fgseaMultilevel, you need to remove the nperm argument in the fgsea function
    ## call.

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (3.91% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize,
    ## gseaParam, : There are duplicate gene names, fgsea may produce unexpected
    ## results.

``` r
ggarrange(plotlist = plot_list, ncol = 1)
```

![](supp_figure_3_files/figure-gfm/fig_S3C-1.png)<!-- -->

``` r
plot_df <- all_gsea_info %>%
  dplyr::filter(grepl("HALLMARK", .$pathway))
plot_df$pathway <- sapply(plot_df$pathway, function(x) sub("HALLMARK_", "", x))
ggplot(plot_df, aes(x = skin_NES, y = gut_NES)) +
  geom_point(data = plot_df[plot_df$gut_NES < 0 | plot_df$skin_NES < 0 |
                                    plot_df$skin_padj < 0.1 | plot_df$gut_NES < 0,], color = "grey") +
  geom_point(data = plot_df[plot_df$gut_padj < 0.1 & plot_df$skin_padj < 0.1 &
                                    plot_df$gut_NES > 0 & plot_df$skin_NES > 0,], color = "red") +
  geom_text_repel(data = plot_df[plot_df$gut_padj < 0.1 & plot_df$skin_padj < 0.1 &
                                    plot_df$gut_NES > 0 & plot_df$skin_NES > 0,], aes(label = pathway), size = 5,
                  min.segment.length = 0.1) +
  xlab("Skin Normalized enrichment score") + ylab("siIEL normalized enrichment score") +
  ggtitle("Hallmark Pathways") +
  theme_classic(base_size = 20)
```

![](supp_figure_3_files/figure-gfm/fig_S3D-1.png)<!-- -->

``` r
plot_df <- all_gsea_info %>%
  dplyr::filter(grepl("KEGG", .$pathway))
plot_df$pathway <- sapply(plot_df$pathway, function(x) sub("KEGG_", "", x))
ggplot(plot_df, aes(x = skin_NES, y = gut_NES)) +
  geom_point(data = plot_df[plot_df$gut_NES < 0 | plot_df$skin_NES < 0 |
                                    plot_df$skin_padj < 0.1 | plot_df$gut_NES < 0,], color = "grey") +
  geom_point(data = plot_df[plot_df$gut_padj < 0.1 & plot_df$skin_padj < 0.1 &
                                    plot_df$gut_NES > 0 & plot_df$skin_NES > 0,], color = "red") +
  geom_text_repel(data = plot_df[plot_df$gut_padj < 0.1 & plot_df$skin_padj < 0.1 &
                                    plot_df$gut_NES > 0 & plot_df$skin_NES > 0,], aes(label = pathway), size = 5,
                  min.segment.length = 0.1) +
  xlab("Skin Normalized enrichment score") + ylab("siIEL normalized enrichment score") +
  ggtitle("KEGG Pathways") +
  theme_classic(base_size = 20)
```

![](supp_figure_3_files/figure-gfm/fig_S3D-2.png)<!-- -->
