# Manually create a list of cluster names ====

cluster_labels <- c("epidermal",
                    "mesophyll group 1",
                    "stressed mesophyll", 
                    "mesophyll group 2",
                    "mesophyll group 3",
                    "bundle sheath",
                    "phloem",
                    "vascular",
                    "stressed epidermal",
                    "companion",
                    "unknown group 1",
                    "unknown group 2",
                    "guard",
                    "sieve",
                    "hydathode",
                    "myrosinase",
                    "unknown group 3")

names(cluster_labels) <- paste("cl", c(0:16), sep = "")

# Load results from Daphne ====

cibersort_res <- read.delim("data/CIBERSORTx_Job145_Results.txt")
row.names(cibersort_res) <- cibersort_res$Mixture

nice_named_cibersort <- cibersort_res[, names(cluster_labels)]
names(nice_named_cibersort) <- cluster_labels

cibersort_by_row <- t(nice_named_cibersort)

# Run MetaCycle on this data (unscaled) ====

if (length(list.files("pre_compute/metacycle_cibersort")) == 0) {
  
  message("Running JTK_cycle on Cibersort data")
  
  # save data without ZT1 for application of metaCycle
  ZT1_names <- c("A1R1", "A1R2", "A1R3",
                 "B1R1", "B1R2", "B1R3",
                 "C1R1", "C1R2", "C1R3")
  
  
  cibersort_wout_ZT1 <- cibersort_by_row[
    , !(colnames(cibersort_by_row) %in% ZT1_names)
  ]
  
  cibersort_A <- cibersort_wout_ZT1[, 1:18]
  write.csv(cibersort_A, "pre_compute/metacycle_cibersort/A.csv", quote = FALSE)
  
  cibersort_B <- cibersort_wout_ZT1[, 19:36]
  write.csv(cibersort_B, "pre_compute/metacycle_cibersort/B.csv", quote = FALSE)
  
  cibersort_C <- cibersort_wout_ZT1[, 37:54]
  write.csv(cibersort_C, "pre_compute/metacycle_cibersort/C.csv", quote = FALSE)
  
  # run metacycle
  mc_timepoints <- rep(seq(0, 20, by = 4), each = 3)
  
  MetaCycle::meta2d(infile = "pre_compute/metacycle_cibersort/A.csv", 
                    filestyle = "csv",
                    outdir = "pre_compute/metacycle_cibersort/", 
                    timepoints = mc_timepoints,
                    cycMethod = c("JTK"), 
                    outIntegration = "noIntegration")
  
  MetaCycle::meta2d(infile = "pre_compute/metacycle_cibersort/B.csv", 
                    filestyle = "csv",
                    outdir = "pre_compute/metacycle_cibersort/", 
                    timepoints = mc_timepoints,
                    cycMethod = c("JTK"), 
                    outIntegration = "noIntegration")
  
  MetaCycle::meta2d(infile = "pre_compute/metacycle_cibersort/C.csv", 
                    filestyle = "csv",
                    outdir = "pre_compute/metacycle_cibersort/", 
                    timepoints = mc_timepoints,
                    cycMethod = c("JTK"), 
                    outIntegration = "noIntegration")
  
  rm(ZT1_names,
     cibersort_A,
     cibersort_B,
     cibersort_C,
     mc_timepoints,
     cibersort_wout_ZT1)

}

# Create data frames for plotting ====

# A lot of this code is based on the plots for gene expression in 0
# but rewrite here so that doesn't have to be modified

cibersort_median_df <- data.frame(matrix(nrow = nrow(cibersort_by_row), ncol = 21))

# colnames are correct in all_means_mat - removing any reference to reps
colnames(cibersort_median_df) <- colnames(all_means_mat)
row.names(cibersort_median_df) <- row.names(cibersort_by_row)

for (row_i in 1:nrow(cibersort_median_df)) {
  for (col_j in 1:ncol(cibersort_median_df)) {
    col_indicies <- c(1,2,3) + 3*(col_j  - 1)
    
    cibersort_median_df[row_i, col_j] <- median(
      cibersort_by_row[row_i, col_indicies]
    )
  }
}

# scale and remove NAs before plotting
# since a few cell types have 0 medians for every timepoint

cibersort_scaled <- na.omit(t(scale(t(cibersort_median_df))))

cibersort_truncated <- apply(cibersort_scaled, c(1,2),
                             function(num_in) {
                               if (num_in >= 3) {
                                 return(3)
                               } else if (num_in <= -3) {
                                 return(-3)
                               } else {
                                 return(num_in)
                               }
                             })

# Plotting heatmaps ====

# truncated data
svg(filename = "outputs/cibersort/heatmap.svg", pointsize = 12,
    height = 7, width = 7.5)
draw(pheatmap(cibersort_truncated, 
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         gaps_col = c(7, 14),
         annotation_col = data.frame("sampling_day" = sapply(
           colnames(cibersort_truncated), function(str_in) {str_sub(str_in, 1, 1)}
         )),
         color = inferno(10),
        annotation_colors = list("sampling_day" = day_colours)
))
dev.off()

# Line plot functions ====

cibersort_plot <- function(cell_type,
                           title_text = cell_type,
                           text_size = 10,
                           point_scaling_factor = 2,
                           line_scaling_factor = 1.25,
                           x_extend = 0.5 # use to expand the x-axis beyond (0,24)
) {
  
  plotting_df <- structure_data_from_labels(cell_type, data_df = cibersort_by_row)
  
  # make median df
  # include `.groups = "keep"` to suppress warning message
  median_by_time_df <- plotting_df %>%
    group_by(time, series) %>%
    summarize(med_value = median(value, na.rm = TRUE), .groups = "keep")
  
  # make basic plot
  out_plot <- ggplot() +
    geom_rect(mapping = aes(xmin = 16, xmax = Inf, ymin = -Inf, ymax = Inf),
              fill = "#cccacd")
  
  out_plot <- out_plot + 
    geom_line(
      aes(x = time, y = med_value, color = series, linetype = series),
      data = median_by_time_df,
      linewidth = line_scaling_factor
    ) +
    geom_point(
      aes(x = time, y = value, color = series, shape = series),
      data = plotting_df,
      size = point_scaling_factor
    ) +
    scale_color_manual(values = day_colours) +
    ggtitle(title_text)
  
  # fix x and y limits, and text size
  out_plot <- out_plot + 
    scale_x_continuous(limits = c((-1) * x_extend, 24 + x_extend), 
                       breaks = seq(0, 24, by = 4),
                       expand = expansion(mult = c(.01, 0))) +
    theme_classic(base_size = text_size) 
  
  # fix text size and content
  out_plot <- out_plot +
    theme(legend.position = "none",
          text = element_text(size = text_size),
          plot.title = element_text(hjust = 0.5)) +
    ylab("Cibersort proportion") +
    xlab("Time (h)") +
    guides(color = guide_legend("Day"),
           linetype = guide_legend("Day"),
           shape = guide_legend("Day"))
  
  return(out_plot)
}

cibersort_all_plots <- lapply(row.names(cibersort_by_row), function(cell_type_in) {
  cibersort_plot(
    cell_type_in,
    text_size = 14, point_scaling_factor = 4
  ) +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.title.x = element_blank())
})

# Plot all cell types ====

cibersort_plot_grid <- cowplot::plot_grid(plotlist = cibersort_all_plots,
                                          align = "hv",
                                          ncol = 3)

svg("outputs/cibersort/all_cell_types.svg", width = 12, height = 16)
print(cibersort_plot_grid)
dev.off()

# Plot only the mesophylls ====

cibersort_mesophylls_plot <- cowplot::plot_grid(
  cibersort_plot("mesophyll group 1"),
  cibersort_plot("mesophyll group 2") +
    theme(axis.title.y = element_blank()),
  cibersort_plot("mesophyll group 3") +
    theme(axis.title.y = element_blank()),
  align = "hv",
  ncol = 3
)

svg(filename = "outputs/cibersort/mesophylls.svg", pointsize = 12,
    height = standard_pheight, width = 3 * standard_pwidth)
print(cibersort_mesophylls_plot)
dev.off()
