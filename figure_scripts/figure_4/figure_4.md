Figure 4
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
library(UpSetR)
library(affy)
library(mogene10sttranscriptcluster.db)
library(fgsea)

godsnot_102 = c(
  # "#FFFF00", Yellow isn't great
  #"#1CE6FF", # This one is kinda ugly too
  "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
  "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
  "#5A0007", "#809693", "#6A3A4C", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
  "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
  "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
  "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
  "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
  "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
  "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
  "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
  "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
  "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
  "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72"
)

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

mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
# mpl.rcParams['pdf.fonttype'] = 42

# Set a colormap
gene_colormap = clr.LinearSegmentedColormap.from_list('gene_cmap', ["#e0e0e1", '#4576b8', '#02024a'], N=200)

# Set a switcher up so the script will run on any computer
def file_path(user = getpass.getuser()):
    switcher = {
            "nealp": "C:/Users/nealp/Documents/Dropbox (Partners HealthCare)/Chloe&Mazen/Collaborator_projects/Kupper_TRM/neal_analysis",
            "neal": "/home/neal/Documents/Dropbox (Partners HealthCare)/Chloe&Mazen/Collaborator_projects/Kupper_TRM/neal_analysis",
            "nealpsmith": "/projects/home/nealpsmith/projects/kupper"

    }
    if switcher.get(user):
        return(switcher.get(user))
    else :
        print("Add your local filepath to the switcher! run getpass.getuser() to get your ID")
```

``` r
padj_cutoff = 0.1
perc_cells_cutoff = 5
slope_cutoff = 0.15

# Get the Kupper gene set
skin_de <- read.csv(glue("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/updated_linear_modeling/skin1_lin_model_by_day_padj_{padj_cutoff}_perc_cells_{perc_cells_cutoff}_slope_{slope_cutoff}.csv"))
gut_de <- read.csv(glue("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/updated_linear_modeling/gut_lin_model_by_day_padj_{padj_cutoff}_perc_cells_{perc_cells_cutoff}_slope_{slope_cutoff}.csv"))

skin_mem <- skin_de[skin_de$log2FoldChange > 0,]
gut_mem <- gut_de[gut_de$log2FoldChange > 0,]

skin_genes <- skin_mem$gene
gut_genes <- gut_mem$gene
kupper_trm <- intersect(skin_genes, gut_genes)

ln_de <- read.csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/pseudobulk/lin_modelling_subclusters/lymph_by_day_all_genes.csv",
                    row.names = 1)

spleen_de <- read.csv("/projects/home/nealpsmith/projects/kupper/kurd_paper/data/pseudobulk/lin_modelling_subclusters/spleen1_by_day_all_genes.csv",
                   row.names = 1)
ln_mem <- ln_de[ln_de$log2FoldChange > 0.15 & ln_de$padj < 0.1 & ln_de$percent_cells > 5,]
spleen_mem <- spleen_de[spleen_de$log2FoldChange > 0.15 & spleen_de$padj < 0.1 & spleen_de$percent_cells > 5,]
ln_genes <- ln_mem$gene
spleen_genes <- spleen_mem$gene
kupper_tcm <- intersect(ln_mem$gene, spleen_mem$gene)

milner_sigs <- read.csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/trm_tcm_signatures_milner_2017.csv")
milner_trm <- milner_sigs$core_trm_signature

mackay_trm <- c("Cd244a", "Cdh1", "Chn2", "Ctla4", "Hpgds", "Hspa1a", "Icos", "Inpp4b", "Itga1",
                "Itgae", "Litaf", "Zfp683", "Pcdh10", "Nr4a1", "Nr4a2", "Qpct", "Rgs1", "Rgs2",
                "Sik1", "Skil", "Tmem123", "Vps37b", "Xcl1")

kupper_trm_filt <- kupper_trm[!kupper_trm %in% kupper_tcm]

upset_list <- list("Trm-consensus" = kupper_trm_filt, "Milner (2017)" = milner_trm, "Mackay (2013)" = mackay_trm)

upset(fromList(upset_list), sets = names(upset_list), order.by = "freq",
      point.size = 7, line.size = 1.5, keep.order = TRUE,
      sets.x.label = "# genes", mainbar.y.label = "# overlapping genes",
      text.scale = c(3.5, 3.5, 2.5, 1.8, 3.5, 0))
grid.text("Trm gene signatures",x = 0.75, y=0.98, gp=gpar(fontsize=25))
```

![](figure_4_files/figure-gfm/fig_4A-1.png)<!-- -->

``` r
data_files <- list.files("/projects/home/nealpsmith/projects/kupper/mackay_2013/data/comp_data", full.names = TRUE)

all_data <- ReadAffy(filenames = data_files)
eset <- rma(all_data)
```

    ## Background correcting
    ## Normalizing
    ## Calculating Expression

``` r
array_df <- data.frame(exprs(eset))
array_df$probe <- rownames(array_df)

Annot <- data.frame(SYMBOL=sapply(contents(mogene10sttranscriptclusterSYMBOL), paste, collapse=","),
                    DESC=sapply(contents(mogene10sttranscriptclusterGENENAME), paste, collapse=","),
                    ENTREZID=sapply(contents(mogene10sttranscriptclusterENTREZID), paste, collapse=","),
                    ENSEMBLID=sapply(contents(mogene10sttranscriptclusterENSEMBL), paste, collapse=","))
Annot$probe <- rownames(Annot)

array_df %<>%
  dplyr::left_join(Annot %>% dplyr::select(SYMBOL, probe), by = "probe")
kupper_sig_data <- array_df[array_df$SYMBOL %in% kupper_trm_filt,]

heatmap_df <- kupper_sig_data %>%
  distinct()
heatmap_df <- heatmap_df[!duplicated(heatmap_df$SYMBOL),]
rownames(heatmap_df) <- heatmap_df$SYMBOL
heatmap_df$SYMBOL <- NULL
heatmap_df$probe <- NULL

colnames(heatmap_df) <- c("Gut_trm1", "Gut_trm2", "Gut_trm3", "Lung_trm1", "Lung_trm2", "Lung_trm3",
                          "naive_1", "naive_2", "naive_3", "Tcm1", "Tcm2", "Tcm3", "Tem1", "Tem2", "Tem3",
                          "skin_trm1", "skin_trm2", "skin_trm3")
# heatmap_df <- heatmap_df[,!grepl("naive", colnames(heatmap_df))]

heatmap_df <- t(scale(t(heatmap_df)))

 # Now the clustering
clustering = hclust(dist(t(heatmap_df), method = "euclidean"), method = "ward.D2")
col_hc <- as.dendrogram(clustering)

# Lets get the top bar
celltype_cols <- list("population" = c("Naive" = "blue",
                                       "Tcm" = "green",
                                       "Tem" = "yellow",
                                       "Trm" = "red"),
                      "location" = c("Skin" = "#a80026",
                                     "Gut" = "#036b38",
                                     "Lung" = "orange",
                                     "Blood" = "#ADD8E6"))

celltypes <- sapply(colnames(heatmap_df), function(n){
  if(grepl("trm", tolower(n))){
    return("Trm")
  } else if (grepl("naive", tolower(n))){
    return("Naive")
  } else if (grepl("tcm", tolower(n))){
    return("Tcm")
  } else if (grepl("tem", tolower(n))){
    return("Tem")
  }
})

# Now annitomical location bar
location <- sapply(colnames(heatmap_df), function(n){
  if(grepl("skin", tolower(n))){
    return("Skin")
  } else if (grepl("gut", tolower(n))){
    return("Gut")
  } else if (grepl("lung", tolower(n))){
    return("Lung")
  } else {
    return("Blood")
  }
})

annotation_bars = HeatmapAnnotation("population" = celltypes,
                                    "location" = location,
                                    col = celltype_cols,
                            show_legend = FALSE, show_annotation_name = FALSE)


# Heatmap colors
heatmap_col_fun = colorRamp2(c(min(heatmap_df), 0, max(heatmap_df)), c("purple", "black", "yellow"))

kupper_unique <- kupper_trm[!kupper_trm %in% c(milner_trm, mackay_trm)]

label_genes <- rownames(heatmap_df)[rownames(heatmap_df) %in% kupper_unique]
core_genes <- c("Xcl1", "Sik1", "Rgs1", "Inpp4b")

label_genes <- c(label_genes, core_genes)

cols <- c()
for (g in label_genes){
  if (g %in% core_genes){
    cols <- c(cols, "red")
  } else {
    cols <- c(cols, "black")
  }
}

# gene_cols <- sapply(rownames(heatmap_df), function(x) ifelse(x %in% kupper_unique, "red", "black"))

# Legends
heatmap_lgd = Legend(col_fun = heatmap_col_fun, title = "z-score", legend_height = unit(4, "cm"),
                     title_position = "topcenter")

pop_lgd <- Legend(labels = names(celltype_cols$population),
                       legend_gp = gpar(fill = unlist(celltype_cols$population)),
                          title = "Population",
                          legend_height = unit(4, "cm"),
                          title_position = "topcenter")
location_lgd <- Legend(labels = names(celltype_cols$location),
                       legend_gp = gpar(fill = unlist(celltype_cols$location)),
                          title = "Location",
                          legend_height = unit(4, "cm"),
                          title_position = "topcenter")

lgd_list <- packLegend(heatmap_lgd, pop_lgd, location_lgd, column_gap = unit(1,"cm"), direction = "vertical")


hmap = Heatmap(heatmap_df, name = "z-score", col = heatmap_col_fun,
               show_column_names = FALSE,
               top_annotation = annotation_bars,
               # row_names_gp = gpar(col = gene_cols, fontsize = 6),
               show_heatmap_legend = FALSE,
               cluster_columns = col_hc) +
  rowAnnotation(link = anno_mark(at = match(label_genes, rownames(heatmap_df)),labels = label_genes,
                                 labels_gp = gpar(col = cols, fontsize = 8, fontface = "bold")))

draw(hmap,heatmap_legend_list =lgd_list)
```

![](figure_4_files/figure-gfm/fig_3G-1.png)<!-- -->

``` python

adata = pg.read_input("/projects/home/nealpsmith/projects/kupper/dealmeida_2022/data/skin_allo_cells.h5ad")
# adata = all_data.to_anndata()
```

    ## 2024-08-30 13:10:03,360 - pegasus - INFO - Time spent on 'read_input' = 0.20s.

``` python
sc.pl.umap(adata, color = ["genotype_anno", "annot"], use_raw = False,
           cmap = gene_colormap, ncols = 2, show = False)

for f in ["genotype_anno", "annot"] :
    print(f)
    col_dict = dict(zip(sorted(set(adata.obs[f])),
                        adata.uns[f"{f}_colors"]))

    legend_elements = [Line2D([0], [0], marker='o', color=col_dict[cl], label=cl,
                              markerfacecolor=col_dict[cl], markersize=7, lw=0) for cl in
                       sorted(set(adata.obs[f]))]
    if f == "genotype_anno" :
        ll = None
    else :
        ll = "on data"
    fig, ax = plt.subplots(1)
    plot = sc.pl.embedding(adata, basis = "umap", color=f,
                      legend_loc=ll, show=False, ax=ax,
                      title="", legend_fontoutline=5)
    lgd = ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    fig = plt.gcf()
    fig.set_size_inches(5, 3)
    fig.tight_layout()
    fig
    plt.close()
```

<img src="figure_4_files/figure-gfm/fig_4C-1.png" width="1514" />

``` python

dot_genes = ["MS4A1", "CD3D", "CD4", "IL7R", "CD8A", "TRDC", "TRGC2", "NKG7", "FOXP3", "TIGIT", "CD52"]
plot = sc.pl.dotplot(adata, dot_genes, groupby="annot",
                     show=False, return_fig=True, title="",
                     categories_order =  ['B cells', 'CD4 T', 'CD8 T', 'GD T', 'NK', 'Treg', 'CD52hi T'],
                     cmap = "Blues", standard_scale = "var")
axes_dict = plot.get_axes()
axes_dict["mainplot_ax"].set_axisbelow(True)
axes_dict["mainplot_ax"].grid()
fig = plt.gcf()
fig.set_size_inches(6, 5)
plt.subplots_adjust(left = 0.2, bottom = 0.2)
fig.tight_layout()
```

    ## /projects/home/nealpsmith/.conda/envs/old_peg_github/bin/python:1: UserWarning: This figure includes Axes that are not compatible with tight_layout, so results might be incorrect.

``` python
fig
```

<img src="figure_4_files/figure-gfm/fig_4c_dots-3.png" width="576" />

``` r
multi_gsea_plot <- function(ranks, gene_sets, title = "", color_list = NULL, annot_offset = 5000){
  # Perform the GSEA
  fgsea_res <- fgsea(gene_sets, stats = ranks, nperm = 10000)

  # Plots
  plot_df <- data.frame()
  pathway_df <- data.frame()
  annot <- c()
  for (gset in names(gene_sets)){
    nes <- round(fgsea_res$NES[fgsea_res$pathway == gset], 3)
    pval <- round(fgsea_res$pval[fgsea_res$pathway == gset], 3)
    n_genes <- fgsea_res$size[fgsea_res$pathway == gset]
    annot <- c(annot, glue("{gset}\n NES : {nes}\n pval : {pval}"))

    rnk <- rank(-ranks)
    ord <- order(rnk)

    statsAdj <- ranks[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ 1)
    statsAdj <- statsAdj / max(abs(statsAdj))

    pathway <- unname(as.vector(na.omit(match(gset_list[[gset]], names(statsAdj)))))
    pathway <- sort(pathway)
    pathway_df <- rbind(pathway_df, data.frame(x = pathway, gs = gset))

    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                              returnAllExtremes = TRUE)

    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops

    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
    toPlot$gset <- gset
    plot_df <- rbind(plot_df, toPlot)

  }
  # col_list <- c("red", "blue")
  # names(col_list) <- names(gset_list)
  p <-ggplot(plot_df, aes(x = x, y = y, group = gset, color = gset)) +
      # geom_point(color="blue", size=0.1) +
      geom_line() +
      geom_hline(yintercept=0, colour="black") +
      scale_y_continuous(expand = c(0.05,0.05)) +
      xlab("Rank") + ylab("Enrichment score") +
      # ggtitle(cl) +
      theme_classic(base_size = 12)
  if (is.null(color_list)){
    color_list <- scales::hue_pal()(length(gene_sets))
    names(color_list) <- names(gene_sets)
  }
  p <- p + scale_color_manual(values= color_list)

  y = max(plot_df$y) + 0.15
  for(an in annot){
    p <- p + annotate("text", label = an, x = length(ranks) - annot_offset, y  =y, hjust = 0)
    y <- y - 0.25
  }
  line_plots = list()
  for (gset in unique(pathway_df$gs)){
    p1 <- ggplot(pathway_df %>% dplyr::filter(gs == gset)) +
      geom_segment(mapping=aes(x=x, y=-0.15,
                               xend=x, yend=-0.25, color = gset),
                   size=0.4) +
      scale_color_manual(values= color_list[[gset]]) +
      theme_classic(base_size = 20) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
            legend.text = element_blank(), legend.title = element_blank()) +
      guides(color = guide_legend(override.aes = list(linetype = 0))) # Funky way to totally hide legend

    line_plots <- c(line_plots, list(p1))
  }
  plot_list <- c(list(p), line_plots)
  all_plots = ggarrange(plotlist = plot_list, ncol = 1, align = "hv", heights = c(1, 0.4, 0.4))
  return(all_plots)
}


fc_data <- read.csv("/projects/home/nealpsmith/projects/kupper/dealmeida_2022/data/cd8_t_host_vs_donor_fc.csv")
trm_genes <- read.csv("/projects/home/nealpsmith/data/useful/kupper_trm_human_genes.csv", row.names = 1)

# Get the Milner gene set
# Lets load in the milner et al gene set
mouse_to_human <- read.csv("/projects/home/nealpsmith/data/useful/human_to_mouse_genes.csv")
milner_sigs <- read.csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/trm_tcm_signatures_milner_2017.csv")

# colnames(mouse_to_human) <- c("HGNC.symbol", "core_trm_signature")
milner_sigs %<>%
  dplyr::left_join(mouse_to_human, by = c("core_trm_signature" = "MGI.symbol"))
colnames(milner_sigs)[colnames(milner_sigs) == "HGNC.symbol"] <- "core_trm_human_genes"

milner_trm <- milner_sigs$core_trm_human_genes
milner_trm <- milner_trm[!is.na(milner_trm)]
gset_list <- list("Temporal Trm" = trm_genes$HGNC.symbol,
                  "milner Trm" = milner_trm)

ranks <- fc_data %>%
  dplyr::filter(perc_donor > 2 | perc_host > 2) %>%
  dplyr::select(featurekey, log_fc_host) %>%
  na.omit() %>%
  arrange(desc(log_fc_host))  %>%
  deframe(.)

fgsea_res <- fgsea(gset_list, stats = ranks, nperm = 10000)
```

    ## Warning in fgsea(gset_list, stats = ranks, nperm = 10000): You are trying to run
    ## fgseaSimple. It is recommended to use fgseaMultilevel. To run fgseaMultilevel,
    ## you need to remove the nperm argument in the fgsea function call.

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (0.02% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

``` r
# Make the multi plot
multi_gsea_plot(ranks = ranks, gene_sets = gset_list,
                color_list = c("Temporal Trm" =  "#7570B3", "milner Trm" = "#E7298A"),
                annot_offset = 2000)
```

    ## Warning in fgsea(gene_sets, stats = ranks, nperm = 10000): You are trying to run fgseaSimple. It is recommended to use fgseaMultilevel. To run fgseaMultilevel, you need to remove the nperm argument in the fgsea function call.
    
    ## Warning in fgsea(gene_sets, stats = ranks, nperm = 10000): There are ties in the preranked stats (0.02% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

![](figure_4_files/figure-gfm/fig_4D-5.png)<!-- -->

``` r
# Which genes are distinct from our leading edge
kupper_le <- fgsea_res[fgsea_res$pathway == "Temporal Trm",]$leadingEdge[[1]]
milner_le <- fgsea_res[fgsea_res$pathway == "milner Trm",]$leadingEdge[[1]]
common <- intersect(kupper_le, milner_le)
kupper_unique <- kupper_le[!kupper_le %in% milner_le]
kupper_unique[!kupper_unique %in% gset_list$`milner Trm`]
```

    ##  [1] "RGCC"    "DDX3X"   "CDKN1A"  "RGS2"    "NFE2L2"  "TSC22D3" "CLK1"   
    ##  [8] "SLC38A2" "AFF4"    "PLCXD2"  "ID3"     "METRNL"  "GADD45A" "PDE7A"  
    ## [15] "BTG1"    "TBX21"

``` r
milner_unique <- milner_le[!milner_le %in% gset_list$`consensus Trm`]
milner_unique[!milner_unique %in% gset_list$`consensus Trm`]
```

    ##  [1] "RGS1"     "DUSP4"    "GZMB"     "PPP1R15A" "DNAJB1"   "CREM"    
    ##  [7] "CCL4"     "DNAJA1"   "DUSP1"    "GZMH"     "JUNB"     "NR4A2"   
    ## [13] "GEM"      "TNFAIP3"  "BTG3"     "CD69"     "TNF"      "DNAJB6"  
    ## [19] "UBE2S"    "ISG20"    "PRDX6"    "ATF3"     "ICOS"     "SPTY2D1" 
    ## [25] "ELL2"     "IFNG"

``` r
draw.pairwise.venn(length(kupper_le), length(milner_le), length(common),
                   category = c("Temporal Trm genes", "Milner genes"), fill = c("#7570B3", "#E7298A"),
                   lwd = rep(0, 2), cex = rep(unit(1, "cm"), 3),
                   cat.pos = c(180, 180), cat.cex = rep(unit(1, "cm"), 2))
```

![](figure_4_files/figure-gfm/fig_4E-1.png)<!-- -->

    ## (polygon[GRID.polygon.514], polygon[GRID.polygon.515], polygon[GRID.polygon.516], polygon[GRID.polygon.517], text[GRID.text.518], text[GRID.text.519], text[GRID.text.520], text[GRID.text.521], text[GRID.text.522])
