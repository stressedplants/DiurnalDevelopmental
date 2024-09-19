# Prepare CV vs. mean expression (log2) ====

# include ALL replicates!
means_all_reps <- rowMeans(log_processed_ge)
names(means_all_reps) <- row.names(log_processed_ge)

sds_all_reps <- rowSds(log_processed_ge)
names(sds_all_reps) <- row.names(log_processed_ge)

relative_sds_all_reps <- sds_all_reps / means_all_reps
names(relative_sds_all_reps) <- row.names(log_processed_ge)

housekeeping_df <- data.frame(
  gene_id = row.names(log_processed_ge),
  mean_expr = means_all_reps,
  cv = relative_sds_all_reps
)

# Load existing housekeeping genes and add categories to DF ====

development_data_2005 <- read_excel(
  path = "data/other_datasets/Czechowski_et_al_2005_Supplement.xls",
  sheet = "Tabelle1",
  range = "A5:H105"
)
development_data_2005$`AGI code` <- stringr::str_to_upper(
  development_data_2005$`AGI code`
)
development_data_2005 <- development_data_2005 %>% filter(
  `AGI code` %in% row.names(log_processed_ge)
)

diurnal_data_2005 <- read_excel(
  path = "data/other_datasets/Czechowski_et_al_2005_Supplement.xls",
  sheet = "Tabelle1",
  range = "BD5:BK105"
)
diurnal_data_2005$`AGI code` <- stringr::str_to_upper(
  diurnal_data_2005$`AGI code`
)
diurnal_data_2005 <- diurnal_data_2005 %>% filter(
  `AGI code` %in% row.names(log_processed_ge)
)

dev_only <- setdiff(development_data_2005$`AGI code`,
                    diurnal_data_2005$`AGI code`)
diu_only <- setdiff(diurnal_data_2005$`AGI code`,
                    development_data_2005$`AGI code`)
both_dev_diu <- intersect(development_data_2005$`AGI code`,
                          diurnal_data_2005$`AGI code`)

labelled_groups <- rep("none", length(housekeeping_df$gene_id))
names(labelled_groups) <- housekeeping_df$gene_id

labelled_groups[dev_only] <- "development"
labelled_groups[diu_only] <- "diurnal"
labelled_groups[both_dev_diu] <- "both"

housekeeping_df$gene_group <- labelled_groups

# Which are the recommended "new" housekeeping genes? ====

# AT1G67200 is a pseudogene so not included in the analysis
manual_gnames <- c("AT3G03070", "ATRPB13.6", "B14/NDUFA6")
names(manual_gnames) <- c("AT3G03070", "AT3G52090", "AT3G12260")

suggested_new_genes_df <- housekeeping_df[names(manual_gnames), ]
suggested_new_genes_df[names(manual_gnames), "gname"] <- manual_gnames

# Which are the two "old" examples? ====

manual_old_gnames <- c("PP2A subunit 3", "YLS8")
names(manual_old_gnames) <- c("AT1G13320", "AT5G08290")

old_housekeeping_df <- housekeeping_df[names(manual_old_gnames), ]
old_housekeeping_df[names(manual_old_gnames), "gname"] <- manual_old_gnames

# Combined plot ====

high_expr_threshold <- log2(11)

colours_for_gene_groups <- brewer.pal(5, name = "Dark2")
hk_point_size <- 1

combined_housekeeping_plot <- ggplot(housekeeping_df, 
                                     aes(x = mean_expr, y = cv,
                                         group = gene_group, colour = gene_group,
                                         shape = gene_group, alpha = gene_group)) +
  # plot this first so it goes in the background
  geom_point(data = housekeeping_df %>% 
               filter(mean_expr >= high_expr_threshold) %>%
               filter(gene_group == "none"),
             size = hk_point_size) +
  # and this in the foreground
  geom_point(data = housekeeping_df %>% 
               filter(mean_expr >= high_expr_threshold) %>%
               filter(gene_group != "none"),
             size = 2*hk_point_size) +
  scale_alpha_manual(values = c("none" = 0.2, "development" = 1, 
                       "diurnal" = 1, "both" = 1)) +
  scale_shape_manual(values = c("none" = 16, "development" = 15, 
                       "diurnal" = 17, "both" = 19)) +
  scale_colour_manual(values = c("none" = "black", 
                                 "development" = colours_for_gene_groups[4], 
                                 "diurnal" = colours_for_gene_groups[2], 
                                 "both" = colours_for_gene_groups[3])) +
  # add labels for two "standard" housekeeping genes (PP2A subunit 3, LHCA2)
  geom_text_repel(data = old_housekeeping_df, aes(label = gname),
                   seed = 123, min.segment.length = 0, nudge_y = -0.3) +
  # add suggested new housekeeping genes
  geom_text_repel(data = suggested_new_genes_df, aes(label = gname),
                   seed = 123, min.segment.length = 0, nudge_y = -0.1) +
  theme_classic() +
  theme(text = element_text(size = standard_text_size)) +
  xlab("Mean expression across all timepoints (log-transformed TPM)") +
  ylab("Coefficient of variation (log scale)") +
  scale_y_log10()

ggsave("outputs/housekeeping/combined_plot.svg", 
       width = large_pwidth, height = large_pheight)

# Save all examples with the same y scales ====

all_example_genes <- rbind(old_housekeeping_df, suggested_new_genes_df)

ggsave(
  "outputs/housekeeping/examples_fixed_y_scale.svg",
  cowplot::plot_grid(
    plotlist = lapply(1:nrow(all_example_genes), function(i) {
      gene_id <- all_example_genes[i, "gene_id"]
      gname <- all_example_genes[i, "gname"]
      combined_days_plot(gene_id, title_text = gname)  +
        theme(axis.title.y = element_blank(),
              axis.title.x = element_blank()) + 
        scale_y_continuous(limits = c(4, 9))
    }),
    ncol = 2
  ),
  width = 2 * standard_pwidth,
  height = 3 * standard_pheight
)

ggsave(
  "outputs/housekeeping/examples_variable_y_scale.svg",
  cowplot::plot_grid(
    plotlist = lapply(1:nrow(all_example_genes), function(i) {
      gene_id <- all_example_genes[i, "gene_id"]
      gname <- all_example_genes[i, "gname"]
      combined_days_plot(gene_id, title_text = gname)  +
        theme(axis.title.y = element_blank(),
              axis.title.x = element_blank())
    }),
    ncol = 2
  ),
  width = 2 * standard_pwidth,
  height = 3 * standard_pheight
)
