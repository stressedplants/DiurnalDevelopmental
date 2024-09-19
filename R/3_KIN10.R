# Components of SnRK1 complex (including KIN10) ====

KINs_lims <- c(3.25, 5.5)

KIN_alphas_plot <- cowplot::plot_grid(
  combined_days_plot("AT3G01090", "KIN10") + 
    scale_y_continuous(limits = KINs_lims),
  combined_days_plot("AT3G29160", "KIN11") + 
    scale_y_continuous(limits = KINs_lims) +
    theme(axis.title.y = element_blank()),
  align = "v"
)

ggsave("outputs/KIN10_related/KIN10_KIN11_plot.svg",
       KIN_alphas_plot,
       width = 2 * standard_pwidth, height = standard_pheight)

rm(KINs_lims, KIN_alphas_plot)

KIN_beta_lims <- c(1.5, 6)

KIN_betas_plot <- cowplot::plot_grid(
  combined_days_plot("AT5G21170", "KINβ1") + 
    scale_y_continuous(limits = KIN_beta_lims),
  combined_days_plot("AT4G16360", "KINβ2") + 
    scale_y_continuous(limits = KIN_beta_lims) +
    theme(axis.title.y = element_blank()),
  combined_days_plot("AT2G28060", "KINβ3") + 
    scale_y_continuous(limits = KIN_beta_lims) +
    theme(axis.title.y = element_blank()),
  align = "v",
  ncol = 3
)

ggsave(
  "outputs/KIN10_related/KIN_betas.svg",
  KIN_betas_plot,
  width = 3*standard_pwidth, height = standard_pheight
)

rm(KIN_beta_lims, KIN_betas_plot)

KIN_gamma_lims <- c(2.75, 6.25)

KIN_gamma_plot <- cowplot::plot_grid(
  combined_days_plot("AT1G09020", "SNF4 / KINβγ") + 
    scale_y_continuous(limits = KIN_gamma_lims),
  combined_days_plot("AT3G48530", "KINγ") + 
    scale_y_continuous(limits = KIN_gamma_lims) +
    theme(axis.title.y = element_blank()),
  align = "v"
)

ggsave("outputs/KIN10_related/KIN_gamma.svg",
       KIN_gamma_plot,
       width = 2 * standard_pwidth, height = standard_pheight)

rm(KIN_gamma_plot, KIN_gamma_lims)

# Set values for KIN10 heatmaps ====

heatmap_KIN10_width <- 6
heatmap_KIN10_height <- 7

# Function to read in KIN10 downstream targets ====

process_kinase_targs <- function(sheet_name) {
  df_out <- read_excel(
    "data/other_datasets/KIN S3.xls",
    sheet = sheet_name
  ) %>% remove_missing()
  colnames(df_out) <- df_out[1,]
  df_out <- df_out[-1,]
  df_out$`AGI number` <- toupper(df_out$`AGI number`)
  
  return(df_out)
}

# Process increasing KIN10 targets ====

kinase_targs_inc <- process_kinase_targs("INCREASED")

# filter out any genes that did not pass initial filtering
inc_gene_expr <- log_median_df[
  intersect(kinase_targs_inc$`AGI number`, row.names(final_gene_abun_ord)), 
] %>% remove_missing()

inc_hc <- hclust(dist(t(scale_over_days(t(inc_gene_expr)))),
                 method = "ward.D2")

inc_cutree <- cutree(inc_hc, k = 2)

KIN10_increasing_heatmap <- scaled_medians_heatmap(
  row.names(inc_gene_expr),
  cluster_rows = inc_hc,
  cutree_rows = 2,
  annotation_row = data.frame("groups" = inc_cutree),
  annotation_col = data.frame("sampling_day" = sapply(
    names(inc_gene_expr), function(str_in) {str_sub(str_in, 1, 1)}
  )),
  annotation_colors = list("groups" = c("1" = "black", "2" = "red"),
                           "sampling_day" = day_colours),
  
)

svg(filename = "outputs/KIN10_related/inc_heatmap.svg",
    pointsize = standard_text_size,
    height = heatmap_KIN10_height,
    width = heatmap_KIN10_width)
draw(KIN10_increasing_heatmap)
dev.off()

KIN10_inc_file <- "outputs/KIN10_related/inc_terms.xlsx"

if (!file.exists(KIN10_inc_file)) {
  
  KIN10_inc_sub <- list(
    "Cluster_1" = names(inc_cutree)[inc_cutree == 1],
    "Cluster_2" = names(inc_cutree)[inc_cutree == 2]
  )
  
  KIN10_inc_GO <- gprofiler2::gost(
    KIN10_inc_sub, 
    organism = "athaliana",
    multi_query = TRUE
  )
  
  write.xlsx2(
    separate_columns(KIN10_inc_GO$result,
                     names(KIN10_inc_sub)),
    KIN10_inc_file,
    row.names = FALSE
  )

  rm(KIN10_inc_GO, KIN10_inc_sub)
}

rm(KIN10_inc_file, kinase_targs_inc,
   inc_gene_expr, inc_hc, inc_cutree,
   KIN10_increasing_heatmap)

# Process decreasing KIN10 targets ====

kinase_targs_dec <- process_kinase_targs("DECREASED")

# filter out any genes that did not pass initial filtering
dec_gene_expr <- log_median_df[
  intersect(kinase_targs_dec$`AGI number`, row.names(final_gene_abun_ord)), 
] %>% remove_missing()

dec_hc <- hclust(dist(t(scale_over_days(t(dec_gene_expr)))),
                 method = "ward.D2")
dec_cutree <- cutree(dec_hc, k = 2)

KIN10_decreasing_heatmap <- scaled_medians_heatmap(
  row.names(dec_gene_expr),
  cluster_rows = dec_hc,
  cutree_rows = 2,
  annotation_row = data.frame("groups" = dec_cutree),
  annotation_col = data.frame("sampling_day" = sapply(
    names(dec_gene_expr), function(str_in) {str_sub(str_in, 1, 1)}
  )),
  annotation_colors = list("groups" = c("1" = "black", "2" = "red"),
                           "sampling_day" = day_colours)
)

svg(filename = "outputs/KIN10_related/dec_heatmap.svg",
    pointsize = standard_text_size,
    height = heatmap_KIN10_height, 
    width = heatmap_KIN10_width)
draw(KIN10_decreasing_heatmap)
dev.off()

KIN10_dec_file <- "outputs/KIN10_related/dec_terms.xlsx"

if (!file.exists(KIN10_dec_file)) {
  
  KIN10_dec_sub <- list(
    "Cluster_1" = names(dec_cutree)[dec_cutree == 1],
    "Cluster_2" = names(dec_cutree)[dec_cutree == 2]
  )
  
  KIN10_dec_GO <- gprofiler2::gost(
    KIN10_dec_sub, 
    organism = "athaliana",
    multi_query = TRUE
  )
  
  write.xlsx2(
    separate_columns(KIN10_dec_GO$result,
                     names(KIN10_dec_sub)),
    KIN10_dec_file,
    row.names = FALSE
  )
  
  rm(KIN10_dec_GO, KIN10_dec_sub)
}

rm(KIN10_dec_file, kinase_targs_dec,
   dec_gene_expr, dec_hc, dec_cutree,
   KIN10_decreasing_heatmap)