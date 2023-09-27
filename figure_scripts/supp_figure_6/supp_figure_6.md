Supplemental figure 6
================

``` r
library(reticulate)
library(gtools)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(glue)
library(magrittr)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(ggrepel)
library(ggpubr)

use_python("/projects/home/nealpsmith/.conda/envs/old_peg_github/bin/python")
```

``` r
# Lets load in the milner et al gene set
milner_sigs <- read.csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/trm_tcm_signatures_milner_2017.csv")

milner_trm <- milner_sigs$core_trm_signature

# Lets look at the overlap with the "core signature" we have from skin + gut
skin_de <- read.csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/pseudobulk/lin_modelling_subclusters/skin1_by_day_padj_0.1_fc_0.2_perc_5.csv",
                    row.names = 1)

gut_de <- read.csv("/projects/home/nealpsmith/projects/kupper/kurd_paper/data/pseudobulk/lin_modelling_subclusters/gut1_by_day_padj_0.1_fc_0.2_perc_5.csv",
                   row.names = 1)

# memory first
skin_mem <- skin_de[skin_de$log2FoldChange > 0,]
gut_mem <- gut_de[gut_de$log2FoldChange > 0,]

skin_genes <- skin_mem$gene
gut_genes <- gut_mem$gene
kupper_trm <- intersect(skin_genes, gut_genes)

milner_kupper_trm_overlap <- intersect(kupper_trm, milner_trm)


# Load in the counts and metadata
count_data <- read.csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/skin1_subcluster_pseudobulk_on_time_counts.csv",
                       row.names = 1)
meta_data <- read.csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/skin1_subcluster_pseudobulk_on_time_meta.csv",
                    row.names = 1)

meta_temp <- meta_data[meta_data$n_cells > 100,]

# Remove the extra D30 sample, seems unfair to only have 2 samples for one timepoint
meta_temp <- meta_temp[rownames(meta_temp) != "samp_19_D30_Skin_30",]
# Need to adjust the int to the lowest point
meta_temp$day_int <- meta_temp$day_int - 2

count_temp <- count_data[,rownames(meta_temp)]

# Make sure the genes are detected in enough samples
n_samp <- rowSums(count_temp != 0)
count_temp <- count_temp[n_samp > round(nrow(meta_temp) / 2),]
count_temp <- count_temp[!grepl("Gm[0-9]|Rpl|Rps|mt-|-ps", rownames(count_temp)),]

# Okay now we can run DESeq
dds <- DESeqDataSetFromMatrix(countData = count_temp,
                              colData = meta_temp,
                              design = ~day_int)
dds <- DESeq(dds)

vsd <- vst(dds, blind=TRUE)
norm_res <- assay(vsd)

# Now lets get the gut data
count_data <- read.csv("/projects/home/nealpsmith/projects/kupper/kurd_paper/data/gut1_data_subcluster_pseudobulk_on_time_counts.csv",
                       row.names = 1)
meta_data <- read.csv("/projects/home/nealpsmith/projects/kupper/kurd_paper/data/gut1_data_subcluster_pseudobulk_on_time_meta.csv",
                    row.names = 1)

meta_temp <- meta_data[meta_data$n_cells > 100,]
count_temp <- count_data[,rownames(meta_temp)]

# Make sure the genes are detected in enough samples
count_temp <- count_temp[!grepl("Gm[0-9]|Rpl|Rps|mt-|-ps", rownames(count_temp)),]

# Okay now we can run DESeq
dds <- DESeqDataSetFromMatrix(countData = count_temp,
                              colData = meta_temp,
                              design = ~day_int)
vsd <- vst(dds, blind=TRUE)
norm_res_kurd <- assay(vsd)

# Okay lets get the overlapping genes
milner_exclusive <- milner_trm[!milner_trm %in% kupper_trm]
found_genes <- intersect(milner_exclusive, intersect(rownames(norm_res), rownames(norm_res_kurd)))

heatmap_mtx <- norm_res[found_genes,]
heatmap_mtx <- t(scale(t(heatmap_mtx)))

colnames(heatmap_mtx) <- sapply(colnames(heatmap_mtx), function(x) strsplit(x, "_")[[1]][3])
heatmap_mtx <- heatmap_mtx[,mixedsort(colnames(heatmap_mtx))]
# Day bar
days <- colnames(heatmap_mtx)

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
                            show_legend = FALSE, show_annotation_name = FALSE)

heatmap_col_fun = colorRamp2(c(min(heatmap_mtx), 0, max(heatmap_mtx)), c("purple", "black", "yellow"))
# Lets make the legends
day_fill <- sapply(levels(days), function(x) day_colors$Day[as.character(x)])
day_skin_legend <- Legend(labels = levels(days), legend_gp = gpar(fill = day_fill), title = "timepoint (Skin)")
heatmap_lgd = Legend(col_fun = heatmap_col_fun, title = "z-score", legend_height = unit(4, "cm"), title_position = "topcenter")


kupper_hmap = Heatmap(heatmap_mtx, name = "z-score", col = heatmap_col_fun,
               top_annotation = day_bar, show_column_names = FALSE,
               show_row_names = FALSE, row_names_gp = gpar(cex = 0.5),
               cluster_columns = FALSE, show_heatmap_legend = FALSE,
                      column_title = "Skin")


### Okay lets add the Gut to the heatmap ###
heatmap_mtx_kurd <- norm_res_kurd[found_genes,]
colname_to_tmpt <- lapply(colnames(heatmap_mtx_kurd), function(x) paste("D", tail(strsplit(x, "_")[[1]], n = 1), sep = ""))
tmpts <- unique(sapply(colnames(heatmap_mtx_kurd), function(x) tail(strsplit(x, "_")[[1]], n = 1)))

heatmap_mtx_kurd <- lapply(tmpts, function(tpt){
  cols = as.data.frame(heatmap_mtx_kurd[,grepl(glue("*_{tpt}$"), colnames(heatmap_mtx_kurd))])
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

heatmap_mtx_kurd <- t(scale(t(heatmap_mtx_kurd)))
heatmap_mtx_kurd <- heatmap_mtx_kurd[,mixedsort(colnames(heatmap_mtx_kurd))]

# Day bar
days <- colnames(heatmap_mtx_kurd)

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

heatmap_col_fun = colorRamp2(c(min(heatmap_mtx_kurd), 0, max(heatmap_mtx_kurd)), c("purple", "black", "yellow"))
# Lets make the legends
day_fill <- sapply(levels(days), function(x) day_colors$Day[as.character(x)])
day_gut_legend <- Legend(labels = levels(days), legend_gp = gpar(fill = day_fill), title = "timepoint (siIEL)")
heatmap_lgd = Legend(col_fun = heatmap_col_fun, title = "z-score", legend_height = unit(4, "cm"), title_position = "topcenter")

lgd_list <- packLegend(heatmap_lgd, day_skin_legend, day_gut_legend, column_gap = unit(1,"cm"), direction = "vertical",
                     max_height = unit(18, "cm"))


kurd_hmap = Heatmap(heatmap_mtx_kurd, name = "z-score", col = heatmap_col_fun,
               top_annotation = day_bar, show_column_names = FALSE,
               show_row_names = TRUE, row_names_gp = gpar(cex = 0.5),
               cluster_columns = FALSE, show_heatmap_legend = FALSE,
                    column_title = "siIEL")

hmap_list <- kupper_hmap + kurd_hmap
draw(hmap_list, heatmap_legend_list = lgd_list)
```

![](supp_figure_6_files/figure-gfm/fig_S6A-1.png)<!-- -->

``` r
### First the siIEL ###
gut1_aucell <- read.csv('/projects/home/nealpsmith/projects/kupper/pyscenic/res/gut1/gut1.auc_with_timepoint.csv', row.names = 1)
colnames(gut1_aucell) <- sapply(colnames(gut1_aucell), function(x) sub("...", "", x, fixed = TRUE))
gut1_aucell$day <- paste("d", gut1_aucell$day, sep = "")
gut_day_int <- list("d4" = 0, "d7" = 1, "d10" = 2, "d14" = 3, "d21" = 4, "d32" = 5, "d60" = 6, "d90" = 7)

# Filter to just regulons that meet criteria
gut1_regulon_genes <-read.csv("/projects/home/nealpsmith/projects/kupper/kurd_paper/data/pyscenic/gut1/gut1_regulon_genes.csv",
                            row.names = 1)
gut_keep <- intersect(colnames(gut1_regulon_genes), colnames(gut1_aucell))
gut1_aucell <- gut1_aucell[,c(gut_keep, "day")]

gut1_aucell$day_int <- sapply(gut1_aucell$day, function(x) gut_day_int[[x]])

mean_auc_gut <- gut1_aucell %>%
  reshape2::melt(id.vars = c("day", "day_int")) %>%
  group_by(day, day_int, variable) %>%
  summarise(mean_aucell = mean(value)) %>%
  mutate(tissue = "gut")

# Use lienar model to calculate p-values
lm_gut <- lapply(unique(mean_auc_gut$variable), function(reg){
  reg_df <- mean_auc_gut %>%
    dplyr::filter(variable == reg)
  lin_model <- summary(lm(mean_aucell ~ day_int, data = reg_df))
  stat_info <- data.frame(regulon = reg,
                          gut_r_squared = lin_model$r.squared,
                          gut_p_val = as.numeric(lin_model$coefficients[,4][2]),
                          gut_slope = as.numeric(lin_model$coefficients[,1][2]))


  return(stat_info)
}) %>%
  do.call(rbind, .)

lm_gut$gut_adj_pval <- p.adjust(lm_gut$gut_p_val, method = "fdr")

### Now the Skin ###
skin1_aucell <- read.csv('/projects/home/nealpsmith/projects/kupper/pyscenic/res/skin1/skin1.auc_with_timepoint.csv', row.names = 1)
colnames(skin1_aucell) <- sapply(colnames(skin1_aucell), function(x) sub("...", "", x, fixed = TRUE))
skin1_aucell$day <- paste("d", skin1_aucell$day, sep = "")
skin_day_int <- list("d5" = 0, "d10" = 1, "d15" = 2, "d20" = 3, "d25" = 4, "d30" = 5, "d45" = 6, "d60" = 7)

skin1_regulon_genes <- read.csv("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/integrated/pyscenic/skin1/skin1_regulon_genes.csv",
                            row.names = 1)
skin_keep <- intersect(colnames(skin1_regulon_genes), colnames(skin1_aucell))
skin1_aucell <- skin1_aucell[,c(skin_keep, "day")]

skin1_aucell$day_int <- sapply(skin1_aucell$day, function(x) skin_day_int[[x]])

mean_auc_skin <- skin1_aucell %>%
  reshape2::melt(id.vars = c("day", "day_int")) %>%
  group_by(day, day_int, variable) %>%
  summarise(mean_aucell = mean(value)) %>%
  mutate(tissue = "skin")

# Use lienar model to calculate p-values
lm_skin <- lapply(unique(mean_auc_skin$variable), function(reg){
  reg_df <- mean_auc_skin %>%
    dplyr::filter(variable == reg)
  lin_model <- summary(lm(mean_aucell ~ day_int, data = reg_df))
  stat_info <- data.frame(regulon = reg,
                          skin_r_squared = lin_model$r.squared,
                          skin_p_val = as.numeric(lin_model$coefficients[,4][2]),
                          skin_slope = as.numeric(lin_model$coefficients[,1][2]))
  return(stat_info)
}) %>%
  do.call(rbind, .)

lm_skin$skin_adj_pval <- p.adjust(lm_skin$skin_p_val, method = "fdr")

# Get rid of the ones that are sig in both for this plot
skin_sig <- lm_skin[lm_skin$skin_adj_pval < 0.1 & lm_skin$skin_slope > 0,]
gut_sig <- lm_gut[lm_gut$gut_adj_pval < 0.1 & lm_gut$gut_slope > 0,]
overlap <- intersect(as.character(skin_sig$regulon), as.character(gut_sig$regulon))
skin_sig_unique <- skin_sig$regulon[!skin_sig$regulon %in% overlap]
gut_sig_unique <- gut_sig$regulon[!gut_sig$regulon %in% overlap]

plot_df <- mean_auc_skin %>%
  dplyr::filter(!variable %in% overlap)

plot_df$sig <- ifelse(plot_df$variable %in% skin_sig_unique, "sig", "not_sig")
plot_df$timepoint <- sapply(plot_df$day, function(x) as.numeric(sub("d", "", x)))
plot_df$timepoint <- factor(plot_df$timepoint)
plot_df %<>%
  group_by(variable) %>%
  mutate(dif_from_start = mean_aucell / mean_aucell[day == "d5"]) %>%
  mutate(label = ifelse(day == "d60" & sig == "sig", as.character(variable), NA))

plot_df$sig <- factor(plot_df$sig, levels = c("not_sig", "sig"))
p_skin <- ggplot(plot_df, aes(x = timepoint, y = dif_from_start, group = variable, color = sig)) +
  geom_point(pch = 21, aes(fill = sig)) +
  scale_fill_manual(values = c("grey", "red")) +
  scale_color_manual(values = c("grey", "red")) +
  geom_line() +
  scale_y_log10() +
  annotation_logticks(side = "l", outside = TRUE) +
  coord_cartesian(clip = "off") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("Day") +
  # ylab("mean AUCell fold-change from baseline") +
  ggtitle("Skin") +
  geom_label_repel(data = plot_df[plot_df$day == "d60",], aes(label = label), color = "black", direction = "x", size = 6) +
  theme_classic(base_size = 20) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

plot_df <- mean_auc_gut %>%
  dplyr::filter(!variable %in% overlap)

plot_df$sig <- ifelse(plot_df$variable %in% gut_sig_unique, "sig", "not_sig")
plot_df$timepoint <- sapply(plot_df$day, function(x) as.numeric(sub("d", "", x)))
plot_df$timepoint <- factor(plot_df$timepoint)
plot_df %<>%
  group_by(variable) %>%
  mutate(dif_from_start = mean_aucell / mean_aucell[day == "d4"])

label_regs <- c("Rb1", "Foxo3", "Myc", "Elf2")
plot_df <- plot_df %>%
  mutate(label = ifelse(day == "d90" & variable %in% label_regs, as.character(variable), NA))

plot_df$sig <- factor(plot_df$sig, levels = c("not_sig", "sig"))
p_gut <- ggplot(plot_df, aes(x = timepoint, y = dif_from_start, group = variable, color = sig)) +
  geom_point(pch = 21, aes(fill = sig)) +
  scale_fill_manual(values = c("grey", "red")) +
  scale_color_manual(values = c("grey", "red")) +
  geom_line() +
  scale_y_log10() +
  annotation_logticks(side = "l", outside = TRUE) +
  coord_cartesian(clip = "off") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("Day") +
  # ylab("mean AUCell fold-change from baseline") +
  ggtitle("siIEL") +
  geom_label_repel(data = plot_df[plot_df$day == "d90",], aes(label = label), color = "black", size = 6) +
  theme_classic(base_size = 20) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

figure <- ggarrange(plotlist = list(p_skin, p_gut), common.legend = TRUE)
annotate_figure(figure, left = textGrob("mean AUCell fold-change from baseline", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                    bottom = textGrob("Day", gp = gpar(cex = 1.3)))
```

![](supp_figure_6_files/figure-gfm/fig_3E-1.png)<!-- -->

``` r
plot_df <- lm_skin %>%
  dplyr::left_join(lm_gut)
plot_df$regulon <- gsub("reg_", "", plot_df$regulon)

ggplot(plot_df, aes(x = skin_slope, y = gut_slope)) +
  geom_point(data = plot_df[plot_df$gut_adj_pval > 0.1 | plot_df$skin_adj_pval > 0.1,], color = "grey", size = 2) +
  geom_point(data = plot_df[plot_df$gut_adj_pval < 0.1 &
                              plot_df$skin_adj_pval < 0.1 & plot_df$skin_slope > 0 &
                              plot_df$gut_slope > 0,], pch = 21, fill = "red", size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab("Skin regression slope") +
  ylab("siIEL regression slope") +
  geom_text_repel(data = plot_df[plot_df$gut_adj_pval < 0.1 & plot_df$skin_adj_pval < 0.1,], aes(label = regulon), size = 8) +
  theme_classic(base_size = 20)
```

![](supp_figure_6_files/figure-gfm/fig_S6C-1.png)<!-- -->
