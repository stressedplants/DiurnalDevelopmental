# Plot all 'core clock' genes as defined in Wang et al 2024 paper  ====

clock_highlights_df <- as.data.frame(
  read_excel("data/gene_lists/wang_et_al_2024_genes.xlsx")
)

clock_genes_plot_list <- lapply(1:nrow(clock_highlights_df), function(row_i) {
  combined_days_plot(clock_highlights_df[row_i, 2], 
                     clock_highlights_df[row_i, 1]) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank())
})

all_clock_grid <- cowplot::plot_grid(plotlist = clock_genes_plot_list,
                                     align = "hv",
                                     ncol = 3)

svg("outputs/clock/core_genes.svg",
    width = 3 * standard_pwidth, height = 7 * standard_pheight)
print(all_clock_grid)
dev.off()

# Plot the three that are important for the story to go into the main figure ====

ggsave(
  "outputs/clock/CCA1_plot.svg",
  plot = combined_days_plot("AT2G46830", "CCA1"),
  width = standard_pwidth,
  height = standard_pheight
)

ggsave(
  "outputs/clock/ELF3_plot.svg",
  plot = combined_days_plot("AT2G25930", "ELF3"),
  width = standard_pwidth,
  height = standard_pheight
)

# Create a main function to run the cluster on different sets of genes and params ====

# use to print out the percentages of diurnal genes on plots

combined_JTK_BH.Q <- data.frame(
  "A_BH.Q" = day_A_JTK$BH.Q,
  "B_BH.Q" = day_B_JTK$BH.Q,
  "C_BH.Q" = day_C_JTK$BH.Q
)
row.names(combined_JTK_BH.Q) <- day_A_JTK$CycID

# needs to be a data.frame so we can pull out single rows and still use `apply`
combined_JTK_signif <- data.frame((combined_JTK_BH.Q < 5e-2))

run_clustering_core_genes <- function(
  target_genes,  # vector - targets of the gene from ChIP-seq data
  output_folder,  # string - prefix of all file outputs - should end with /
  n_clusters,  # integer - number of clusters to produce
  color_palette = "Spectral",  # string - set colours to use for brewer.pal,
  line_alpha = 0.25,  # applies to timeseries plots
  base_size = 11,  # applies to timeseries plots
  timeseries_linewidth = 0.25,
  mean_linewidth = 0.25,
  custom_cluster_order = 1:n_clusters,  # change timeseries plotting order
  ...  # params passed to ggplot2::theme for the individual timeseries plots
) {
  
  # 0) Prep data
  filtered_targets <- intersect(target_genes, row.names(final_gene_abun_ord))
  expr_df <- t(scale_over_days(t(log_median_df[filtered_targets,])))
  dist_mat <- dist(expr_df)
  
  # 1) Create a hierarchical clustering tree and cut
  hc_output <- hclust(dist_mat, method = "ward.D2")
  cutree_output <- cutree(hc_output, k = n_clusters)
  
  # 2) Print a heatmap
  cluster_colors <- brewer.pal(n = n_clusters, name = color_palette)
  names(cluster_colors) <- as.character(1:n_clusters)
  
  heatmap_output <- scaled_medians_heatmap(
    filtered_targets,
    cluster_rows = gclus::reorder.hclust(hc_output, dist_mat),
    cutree_rows = n_clusters,
    annotation_row = data.frame("groups" = as.character(cutree_output)),
    annotation_col = data.frame("sampling_day" = sapply(
      names(log_median_df), function(str_in) {str_sub(str_in, 1, 1)}
    )),
    annotation_colors = list("groups" = cluster_colors,
                             "sampling_day" = day_colours)
  )
  
  # save both an svg and png
  svg(filename = paste0(output_folder, "heatmap.svg"),
      width = 6, height = 6, pointsize = 12)
  draw(heatmap_output)
  dev.off()
  
  png(filename = paste0(output_folder, "heatmap.png"),
      width = 6, height = 6, units = "in", res = 600, pointsize = 12)
  draw(heatmap_output)
  dev.off()
  
  # 3) Print timeseries plots for each cluster
  total_plot_list <- list()
  
  # helpful for removing axes in plots
  blank_x <- theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  
  blank_y <- theme(axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  
  for (plot_i in 1:n_clusters) {
    
    group_i <- custom_cluster_order[plot_i]
    
    # find vec of genes from this group
    group_genes <- names(cutree_output)[cutree_output == group_i]
    
    # calculate the percentage of genes which are diurnal in each group
    group_JTK <- combined_JTK_signif[group_genes, ]
    
    proportion_diurnal <- apply(group_JTK, 2, 
                                function(col_in) {sum(col_in) / length(col_in)})
    
    # generate list of plots per day
    ts_plot_list <- scaled_timeseries_groups( 
      group_genes,
      line_alpha = line_alpha,
      base_size = base_size,
      timeseries_linewidth = timeseries_linewidth,
      mean_linewidth = mean_linewidth,
      ...
    )
    
    # if not the final plot, remove all the x axes
    if (plot_i != n_clusters) {
      ts_plot_list <- lapply(ts_plot_list, function(plot_in) {
        plot_in + blank_x
      })
    }
    
    # add new plot within a list - otherwise the ggplot object is expanded
    total_plot_list <- c(
      total_plot_list, 
      list(ts_plot_list[[1]] + 
        annotate(geom = 'text', 
                 label = paste0(
                   sprintf(proportion_diurnal[["A_BH.Q"]]*100, fmt = "%#.1f"), 
                   "%"
                 ), 
                 x = Inf, y = -Inf, hjust = 1.2, vjust = -0.5)
    ))
    
    total_plot_list <- c(
      total_plot_list, 
      list(ts_plot_list[[2]] + blank_y +
        annotate(geom = 'text', 
                 label = paste0(
                   sprintf(proportion_diurnal[["B_BH.Q"]]*100, fmt = "%#.1f"), 
                   "%"
                 ), 
                 x = Inf, y = -Inf, hjust = 1.2, vjust = -0.5)
    ))
    
    total_plot_list <- c(
      total_plot_list, 
      list(ts_plot_list[[3]] + blank_y +
        annotate(geom = 'text', 
                 label = paste0(
                   sprintf(proportion_diurnal[["C_BH.Q"]]*100, fmt = "%#.1f"), 
                   "%"
                 ), 
                 x = Inf, y = -Inf, hjust = 1.2, vjust = -0.5)
    ))

  }
  
  timeseries_joint_plot <- cowplot::plot_grid(
    plotlist = total_plot_list,
    align = "hv",
    ncol = 3
  )

  ggsave(
    paste0(output_folder, "ts_plot.svg"),
    timeseries_joint_plot,
    width = 5,
    height = 5.5
  )
  
  # 4) Run GO term analysis on each cluster
  gprofiler_sub <- lapply(1:n_clusters, function(clust_i) {
    names(cutree_output)[cutree_output == clust_i]
  })
  names(gprofiler_sub) <- sapply(1:n_clusters, function(clust_i) {
    paste0("cluster_", clust_i)
  })
  
  GO_output <- gprofiler2::gost(
    gprofiler_sub, 
    organism = "athaliana",
    multi_query = TRUE
  )
  
  write.xlsx2(
    separate_columns(GO_output$result,
                     names(gprofiler_sub)),
    paste0(output_folder, "GO_submission.xlsx"),
    row.names = FALSE
  )
  
  # 5) Save xlsx with the cluster allocation for each gene
  write.xlsx2(
    as.data.frame(sort(cutree_output)), 
    file = paste0(output_folder, "cluster_allocations.xlsx")
  )
}

# Run for CCA1 and ELF3 - they produce interesting outputs ====

# 5 clusters is a balance between separating highly rhythmic and non-rhythmic 
# genes while keeping the clusters large enough for gProfiler

ELF3_df <- read_excel("data/chip_seq/ELF3_ELF4_LUX_Ezer_2017.xlsx",
                      sheet = "ELF3_22")

if (length(list.files("outputs/clock/ELF3/")) == 0) {
  run_clustering_core_genes(
    ELF3_df$`Locus Identifier`,
    "outputs/clock/ELF3/",
    5,
    custom_cluster_order = c(5, 2, 1, 4, 3),
    plot.margin = margin(l = -0.3, b = -0.3, t = -0.3, unit = "in")
  )
}

CCA1_col_types <- rep("guess", 32)
CCA1_col_types[2] <- "text"

CCA1_LD_df <- read_excel("data/chip_seq/CCA1_Nagel_2015_LD.xlsx",
                         skip = 2,
                         col_names = TRUE,
                         col_types = CCA1_col_types)

if (length(list.files("outputs/clock/CCA1/")) == 0) {
  run_clustering_core_genes(
    tools::file_path_sans_ext(CCA1_LD_df$`Nearest PromoterID`),
    "outputs/clock/CCA1/",
    5,
    line_alpha = 0.09,
    custom_cluster_order = c(2, 5, 1, 4, 3),
    plot.margin = margin(l = -0.3, b = -0.3, t = -0.3, unit = "in")
  )
}
