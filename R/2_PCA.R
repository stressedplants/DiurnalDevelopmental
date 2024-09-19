# Run PCA for further analysis ====

pca_output_final <- prcomp(t(log_processed_ge))
plant_pca_summary <- summary(pca_output_final)$importance

# Save a scree plot ====

scree_df <- data.frame(
  PC_num = 1:10,
  PC_var_percent = plant_pca_summary[2, 1:10]*100
)

scree_plot <- ggplot(scree_df, aes(x = PC_num, y = PC_var_percent)) +
  geom_point(size = 3) + 
  geom_line() +
  theme_classic(base_size = large_base_size) +
  scale_x_continuous(name = "Principal component",
                     breaks = 1:10) +
  ylab("Percentage variance explained")

ggsave("outputs/pca/scree_plot.svg", scree_plot,
       width = large_pwidth, height = large_pheight)

rm(scree_df, scree_plot)

# Prepare PCA data for plotting =====

pca_plant_for_plot <- data.frame(pca_output_final$x[, 1:3])

# create colours for different categories
sample_categories <- sapply(row.names(pca_plant_for_plot), function(i) {
  substring(i, first = 1, last = 1)
})
names(sample_categories) <- sample_categories

pca_plant_for_plot$cats <- as.factor(sample_categories)

timepoints_ordered <- str_match(
  row.names(pca_plant_for_plot), "[ABC](.*?)R"
)[,2]

pca_plant_for_plot$timepoint <- as.factor(timepoints_ordered)

rm(sample_categories, timepoints_ordered)

# PC1 / PC2 samples ====

PC1_PC2_labelled <- ggplot(pca_plant_for_plot, aes(x = PC1, y = PC2)) + 
  geom_vline(aes(xintercept = 0)) + 
  geom_hline(aes(yintercept = 0)) +
  theme_minimal(base_size = standard_text_size) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_point(aes(color = cats, shape = cats)) +
  scale_x_continuous(limits = c(-110, 120)) +
  scale_color_manual(values = day_colours) +
  xlab(paste0("PC1 (", 
              sprintf(plant_pca_summary[2, 1]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC2 (",
              sprintf(plant_pca_summary[2, 2]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  geom_text_repel(aes(label = timepoint),
                  color = dark_labels_col,
                  max.overlaps = 70,
                  box.padding = 0,
                  force = 1.5,
                  data = pca_plant_for_plot,
                  seed = 123,
                  size = 3.75)

ggsave("./outputs/pca/PC1_PC2_labelled.svg", PC1_PC2_labelled,
       width = standard_pwidth, 
       height = standard_pheight, 
       bg = "white")

rm(PC1_PC2_labelled)

# PC2 / PC3 samples ====

PC2_PC3_labelled <- ggplot(pca_plant_for_plot, aes(x = PC2, y = PC3)) + 
  geom_vline(aes(xintercept = 0)) + 
  geom_hline(aes(yintercept = 0)) +
  theme_minimal(base_size = standard_text_size) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_point(aes(color = cats, shape = cats)) +
  scale_color_manual(values = day_colours) +
  xlab(paste0("PC2 (", 
              sprintf(plant_pca_summary[2, 2]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC3 (",
              sprintf(plant_pca_summary[2, 3]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  geom_text_repel(aes(label = timepoint),
                  color = dark_labels_col,
                  max.overlaps = 70,
                  box.padding = 0,
                  force = 1.5,
                  data = pca_plant_for_plot,
                  seed = 123,
                  size = 3.75)

ggsave("./outputs/pca/PC2_PC3_labelled.svg", PC2_PC3_labelled,
       width = standard_pwidth, 
       height = standard_pheight, 
       bg = "white")

rm(PC2_PC3_labelled)

# Prepare data for PCA loadings plots ====

pca_loadings_for_plot <- as.data.frame(
  pca_output_final$rotation[, 1:3]
)

PC1_quantiles <- quantile(pca_loadings_for_plot[, 1], 
                               probs = c(0.1, 0.9))

# Function for plotting about loadings ====

highlight_genes_on_loadings <- function(
    ids_to_highlight,
    first_pc = 1,
    second_pc = 2,
    text_size = standard_text_size,
    gene_names = NULL,  # custom since the automatic gene names can be bad,
    repel_force = 3,
    label_seed = 123,
    label_size = 4,
    axes_width = 0.3,
    point_size = 0.6,
    add_PC1_percentiles = FALSE # create vertical lines for PC1 quantiles
) {
  
  # prepare separate data frames for background and highlighted genes
  pca_loadings_colours <- rep("black", nrow(pca_output_final$rotation))
  names(pca_loadings_colours) <- row.names(pca_output_final$rotation)
  pca_loadings_colours[ids_to_highlight] <- "red"
  pca_loadings_for_plot$col <- pca_loadings_colours
  
  # filtered data for labels
  labels_loadings <- pca_loadings_for_plot[ids_to_highlight, ]
  
  if (!is.null(gene_names)) labels_loadings$gnames <- gene_names
  
  background_gene_df <- pca_loadings_for_plot %>% filter(col == "black")
  highlight_gene_df <- pca_loadings_for_plot %>% filter(col == "red")
  
  # make plot (with two separate geom_point funcs) and remove grid lines
  pca_loadings_plot <- ggplot(pca_loadings_for_plot,
                              aes(x = .data[[paste0("PC", first_pc)]], 
                                  y = .data[[paste0("PC", second_pc)]])) +
    theme_minimal(base_size = text_size) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    geom_vline(aes(xintercept = 0), linewidth = axes_width) + 
    geom_hline(aes(yintercept = 0), linewidth = axes_width)
  
  pca_loadings_plot <- pca_loadings_plot +
    geom_vline(aes(xintercept = 0)) + 
    geom_hline(aes(yintercept = 0)) +
    geom_point(data = background_gene_df,
               aes(x = .data[[paste0("PC", first_pc)]], 
                   y = .data[[paste0("PC", second_pc)]]),
               shape = 19,
               size = point_size,
               colour = "black",
               alpha = 0.2) +
    geom_point(data = highlight_gene_df,
               aes(x = .data[[paste0("PC", first_pc)]], 
                   y = .data[[paste0("PC", second_pc)]]),
               shape = 19,
               size = point_size,
               colour = "red") +
    xlab(paste0("PC",
                first_pc, 
                " (", 
                sprintf(plant_pca_summary[2, first_pc]*100, fmt = "%#.1f"),
                "% of variance)")) +
    ylab(paste0("PC",
                second_pc,
                " (",
                sprintf(plant_pca_summary[2, second_pc]*100, fmt = "%#.1f"),
                "% of variance)"))
  
  if (add_PC1_percentiles) {
    if (!(first_pc == 1)) stop("Must use PC1 to add percentiles")
    pca_loadings_plot <- pca_loadings_plot +
      geom_vline(xintercept = PC1_quantiles[1], 
                 linetype = "dashed", colour = "red") +
      geom_vline(xintercept = PC1_quantiles[2], 
                 linetype = "dashed", colour = "red")
  }

  if (!is.null(gene_names)) {
    pca_loadings_plot <- pca_loadings_plot +
      geom_text_repel(aes(label = gnames),
                      data = labels_loadings, 
                      force = repel_force,
                      max.overlaps = length(gene_names),
                      size = label_size,
                      seed = label_seed,
                      colour = "red",
                      segment.color = "red")
  }
  
  return(pca_loadings_plot)
}

# PC1 and PC2 loadings plot ====

highlights_df <- read_excel("data/gene_lists/PC1_genes_to_highlight.xlsx")

PC1_PC2_highlights_plot <- highlight_genes_on_loadings(
  highlights_df$`TAIR ID`,
  gene_names = highlights_df$`Gene short name`,
  add_PC1_percentiles = TRUE
)

ggsave("outputs/pca/PC1_PC2_highlighted_genes.svg", PC1_PC2_highlights_plot,
       width = standard_pwidth, height = standard_pheight)

rm(PC1_PC2_highlights_plot, highlights_df)

# PC2 and PC3 loadings plot ====

clock_highlights_df <- read_excel("data/gene_lists/wang_et_al_2024_genes.xlsx")

PC2_PC3_highlights_plot <- highlight_genes_on_loadings(
  clock_highlights_df$`TAIR ID`,
  first_pc = 2,
  second_pc = 3,
  gene_names = clock_highlights_df$`Gene short name`,
  label_size = 3,
  repel_force = 3
)

ggsave("outputs/pca/PC2_PC3_highlighted_genes.svg", PC2_PC3_highlights_plot,
       width = standard_pwidth, height = standard_pheight)

rm(clock_highlights_df,
   PC2_PC3_highlights_plot)

# The top / bottom genes from PC1 / PC2 loadings ====

CCA1_plot <- combined_days_plot("AT2G46830", "CCA1")

ELF4_plot <- combined_days_plot("AT2G40080", "ELF4")

PR1_plot <- combined_days_plot("AT2G14610", "PR1")

AT2G10940_plot <- combined_days_plot("AT2G10940", "AT2G10940")

eg_plots <- cowplot::plot_grid(
  PR1_plot,
  ELF4_plot,
  AT2G10940_plot,
  CCA1_plot,
  ncol = 2,
  align = 'hv'
)

ggsave(
  "outputs/pca/Fig1_eg.svg",
  eg_plots,
  width = standard_pwidth * 2, height = standard_pheight * 2
)

rm(eg_plots,
   CCA1_plot,
   ELF4_plot,
   PR1_plot,
   AT2G10940_plot)

# GO term overrep of the top / bottom 10% of genes by PC1 loading ====

# this depends on external database 
# so don't rerun if it already exists in the outputs folder

GO_query_file <- "outputs/pca/GO_query.xlsx"

if (!file.exists(GO_query_file)) {
  
  message("Running gProfiler to get PC1 high / low loadings GO analysis")
  
  gene_ids_pca_loadings <- row.names(pca_loadings_for_plot)
  
  PC1_low <- gene_ids_pca_loadings[
    which(pca_loadings_for_plot[, 1] <= PC1_quantiles[1]) 
  ]
  
  PC1_high <- gene_ids_pca_loadings[
    which(pca_loadings_for_plot[, 1] >= PC1_quantiles[2]) 
  ]
  
  pca_GO_term_submission <- list(
    "PC1_low" = PC1_low,
    "PC1_high" = PC1_high
  )
  
  pca_GO_output_multiquery <- gprofiler2::gost(
    pca_GO_term_submission, 
    organism = "athaliana", 
    multi_query = TRUE
  )
  
  write.xlsx2(
    separate_columns(pca_GO_output_multiquery$result,
                     names(pca_GO_term_submission)),
    GO_query_file,
    row.names = FALSE
  )
  
  pca_GO_output_shortlink <- gprofiler2::gost(
    pca_GO_term_submission, 
    organism = "athaliana", 
    multi_query = TRUE,
    as_short_link = TRUE
  )
  
  write_lines(pca_GO_output_shortlink, "outputs/pca/GO_query_link.txt")
  
  rm(PC1_low, PC1_high, pca_GO_output_shortlink,
     pca_GO_output_multiquery, pca_GO_term_submission)
  
}

rm(GO_query_file)

# Function to add histograms to top / bottom of PCA loading plot ====

histograms_on_highlights <- function(
  genes_to_highlight,
  first_pc = 1,
  second_pc = 2,
  num_top_breaks = 4,
  num_side_breaks = 3,
  density_alpha = .4,
  ... # extra params passed to highlight_genes_on_loadings
) {
  scatter_plot <- highlight_genes_on_loadings(
    genes_to_highlight,
    first_pc = first_pc,
    second_pc = second_pc,
    ...
  )
  
  x_range <- ggplot_build(scatter_plot)$layout$panel_params[[1]]$x.range
  y_range <- ggplot_build(scatter_plot)$layout$panel_params[[1]]$y.range
  
  # get loadings specifically for highlighted genes - for histograms
  highlighted_loadings <- pca_loadings_for_plot[genes_to_highlight, ]
  
  hist_PC1 <- ggplot(highlighted_loadings, aes(x = .data[[paste0("PC", first_pc)]])) + 
    geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "white") +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_density(alpha = density_alpha, fill = "red") +
    scale_x_continuous(limits = x_range, expand = c(0, 0)) +
    scale_y_continuous(n.breaks = num_top_breaks) +
    theme_minimal() +
    theme(axis.title.x = element_blank())
  
  hist_PC2 <- ggplot(highlighted_loadings, aes(x = .data[[paste0("PC", second_pc)]])) +
    geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "white") +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_density(alpha = density_alpha, fill = "red") +
    scale_x_continuous(limits = y_range, expand = c(0, 0)) +
    scale_y_continuous(n.breaks = num_side_breaks) +
    theme_minimal() +
    theme(axis.title.y = element_blank()) +
    coord_flip()
  
  combined_plot <- cowplot::plot_grid(
    hist_PC1, NULL, scatter_plot, hist_PC2,
    align = 'hv',
    rel_widths = c(3.5, 1),
    rel_heights = c(1, 3.5)
  )
  
  return(combined_plot)
}

# Plot senescence / kinase related genes on PCA ====

# find genes annotated with a few GO terms...

annotated_query_file <- "outputs/loading_distributions/annotated_genes.Rdata"

if (!file.exists(annotated_query_file)) {
  
  message("Running gProfiler to retrieve genes associated to GO terms")
  
  # protein serine / threonine kinase activity
  ser_thr_kinase_df <- gprofiler2::gconvert("GO:0004674", 
                                            organism = "athaliana")
  
  # leaf senescence
  leaf_sen_df <- gprofiler2::gconvert("GO:0010150", organism = "athaliana")
  
  # response to light stimulus
  light_stim_df <- gprofiler2::gconvert("GO:0009416", "athaliana")
  
  # photosynthesis
  photosynthesis_df <- gprofiler2::gconvert("GO:0015979", "athaliana")
  
  save(ser_thr_kinase_df, 
       leaf_sen_df,
       light_stim_df,
       photosynthesis_df,
       file = annotated_query_file)
} else {
  load(annotated_query_file)
}

# ... that also overlap with filtered data
ser_thr_intersect <- intersect(ser_thr_kinase_df$target, 
                               row.names(pca_loadings_for_plot))

# load the 4 high PC1 kinases to highlight
kinase_highlight_df <- as.data.frame(
  read_excel("data/gene_lists/kinases_to_highlight.xlsx")
)
row.names(kinase_highlight_df) <- kinase_highlight_df$ID

# set NA names if they are not in the kinase_highlight_df - so the labels 
# don't show up on the plot
kinase_names_including_NA <- as.vector(
  kinase_highlight_df[ser_thr_intersect, "Name"]
)

ggsave(
  "outputs/loading_distributions/ser_thr_plot_PC1_PC2.svg",
  histograms_on_highlights(ser_thr_intersect, 
                           gene_names = kinase_names_including_NA,
                           first_pc = 1,
                           second_pc = 2),
  width = large_pwidth, 
  height = large_pheight
)


ggsave(
  "outputs/loading_distributions/ser_thr_plot_PC2_PC3.svg",
  histograms_on_highlights(ser_thr_intersect, 
                           gene_names = kinase_names_including_NA,
                           first_pc = 2,
                           second_pc = 3),
  width = large_pwidth, 
  height = large_pheight
)

# Run t-tests to check for difference to 0 ====

t.test(pca_loadings_for_plot[ser_thr_intersect, "PC1"])
t.test(pca_loadings_for_plot[ser_thr_intersect, "PC2"])
t.test(pca_loadings_for_plot[ser_thr_intersect, "PC3"])

# Other kinase examples ====

WAK_CRK_examples <- cowplot::plot_grid(
  combined_days_plot(kinase_highlight_df[1, "ID"],
                     kinase_highlight_df[1, "Name"]) +
    theme(axis.title.x = element_blank()),
  combined_days_plot(kinase_highlight_df[2, "ID"], 
                     kinase_highlight_df[2, "Name"]) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank()),
  combined_days_plot(kinase_highlight_df[3, "ID"],
                     kinase_highlight_df[3, "Name"]),
  combined_days_plot(kinase_highlight_df[4, "ID"], 
                     kinase_highlight_df[4, "Name"]) +
    theme(axis.title.y = element_blank()),
  align = "hv"
)

ggsave("outputs/loading_distributions/WAK_CRK_examples.svg",
       WAK_CRK_examples,
       width = standard_pwidth * 2, height = standard_pheight * 2)

rm(ser_thr_intersect, ser_thr_kinase_df,
   WAK_CRK_examples,
   kinase_highlight_df, kinase_names_including_NA)

# Remove PCA related things ====

rm(pca_loadings_for_plot,
   plant_pca_summary,
   pca_output_final,
   pca_plant_for_plot,
   highlight_genes_on_loadings,
   histograms_on_highlights,
   annotated_query_file,
   PC1_quantiles)

   