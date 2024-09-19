# basic clustering and filtering of samples ========

correlation_and_filtering_file <- "pre_compute/correlation_and_filtering.Rdata"

if (!file.exists(correlation_and_filtering_file)) {
  
  # do initial clustering based on correlation
  cor_gene_abun_ord <- apply(log_gene_abun_ord, 2, function(i) {
    apply(log_gene_abun_ord, 2, function(j){
      cor(i, j)   
    })
  })
  
  diag(cor_gene_abun_ord) <- NA
  svg("./outputs/filtering/initial_correlation_heatmap.svg",
      width = 16,
      height = 16)
  heatmap_initial <- pheatmap(
    cor_gene_abun_ord,
    color = inferno(10),
    cellwidth = 13, cellheight = 13,
    fontsize = 15
  )
  draw(heatmap_initial)
  dev.off()

  # low filtering with non-zero samples
  
  low_threshold_lims <- seq(from = 0, to = 2, by = 0.1)
  num_removed_low <- sapply(low_threshold_lims, function(i){
    filtered_expr_df <- as.data.frame(median_df) %>%
      filter_all(all_vars(. <= i))
    return(nrow(filtered_expr_df))
  })
  
  low_df <- data.frame(limits = low_threshold_lims, nums = num_removed_low)
  
  low_plot <- ggplot(low_df, aes(x = limits, y = nums)) +
    geom_point() +
    geom_vline(xintercept = 0.5, color = 'red') +
    theme_minimal(base_size = 12) +
    theme(axis.line = element_line(linewidth = 1, colour = "black")) +
    xlab("Lower threshold limit") + 
    ylab("Number of genes removed from further analysis")
  ggsave("outputs/filtering/low_threshold.svg", plot = low_plot, height = 5)
  
  # finally remove all genes with all TPM values less than 0.5 ...
  final_low_removed_genes <- row.names(
    as.data.frame(median_df) %>% filter_all(all_vars(. <= 0.5))
  )
  
  message(paste0("Removed ", length(final_low_removed_genes), 
                 " genes before further analysis, due to low expression "))
  
  final_gene_abun_ord <- gene_abun_ord[
    !(row.names(gene_abun_ord) %in% final_low_removed_genes),
  ]
  
  # do final clustering heatmap ====
  
  log_processed_ge <- log2(final_gene_abun_ord + 1)
  
  cor_post_final <- apply(log_processed_ge, 2, function(i) {
    apply(log_processed_ge, 2, function(j){
      cor(i, j)   
    })
  })
  diag(cor_post_final) <- NA
  
  # set up data frames to annotate the rows and columns by day
  anno_days_df <- data.frame(
    day_sampled = sapply(row.names(cor_post_final), function(i) {
      substring(i, first = 1, last = 1)
    })
  )
  row.names(anno_days_df) <- row.names(cor_post_final)
  
  svg("./outputs/filtering/after_filtering_correlation_heatmap.svg",
      width = 16,
      height = 16)
  heatmap_post_final <- pheatmap(
    cor_post_final,
    color = inferno(10),
    cellwidth = 13, cellheight = 13,
    fontsize = 15,
    annotation_row = anno_days_df,
    annotation_col = anno_days_df,
    annotation_colors = list(day_sampled = day_colours)
  )
  draw(heatmap_post_final)
  dev.off()
  
 save(final_gene_abun_ord,
      log_processed_ge,
      file = correlation_and_filtering_file)
  rm(cor_gene_abun_ord,
     cor_post_final,
     log_gene_abun_ord,
     low_df,
     low_plot,
     final_low_removed_genes,
     low_threshold_lims,
     anno_days_df,
     num_removed_low,
     heatmap_initial,
     heatmap_post_final)
} else {
  message("Loading correlation and filtering data from file")
  load(correlation_and_filtering_file)
}

# download relevant biomaRt data and store this for loading ======

# Load the file if it already exists - to be consistent across runs
biomart_data_file <- "./pre_compute/biomart_data.Rdata"

if (!file.exists(biomart_data_file)) {
  message("Downloading gene names and GO terms from TAIR")
  
  # BiomaRt set up - only need to run once and then save
  ensembl_arabidopsis <- useEnsemblGenomes(biomart = "plants_mart",
                                           dataset = "athaliana_eg_gene")
  
  attributes_to_retrieve <- c('external_gene_name', 'external_gene_source',
                              'external_transcript_name', 
                              'external_transcript_source_name',
                              'external_synonym',
                              'tair_locus', 
                              'go_id', 'name_1006', 'definition_1006',
                              'go_linkage_type', 'namespace_1003')
  
  genes_to_retrieve <- row.names(gene_abun_ord)
  
  BM_download <- getBM(attributes = attributes_to_retrieve,
                       filters = 'tair_locus',
                       values = genes_to_retrieve,
                       mart = ensembl_arabidopsis)
  
  message("Creating GO term dataframe")
  # some GO terms do not have names attached - just remove these for now
  # but could manually add these later
  GO_names_df <- unique(
    BM_download[, c("go_id", "name_1006", "namespace_1003")]
  ) %>% filter(go_id != "") %>% filter(name_1006 != "")
  GO_names_df$namespace_1003 <- as.factor(GO_names_df$namespace_1003)
  
  row.names(GO_names_df) <- GO_names_df$go_id
  
  message("Producing gene name vec")
  BM_simplified <- distinct(
    BM_download[, c('external_gene_name', 'tair_locus')]
  )
  row.names(BM_simplified) <- BM_simplified$tair_locus
  
  gene_name_id <- sapply(row.names(gene_abun_ord), function(gene_id){
    attempt_gene_name <- BM_simplified[gene_id, 'external_gene_name']
    if (is.na(attempt_gene_name)) {
      return(gene_id)
    } else {
      return(attempt_gene_name) 
    }
  })
  
  message("Saving BiomaRt related files")

  save(BM_download,
       GO_names_df,
       gene_name_id,
       file = biomart_data_file)
  rm(ensembl_arabidopsis,
     attributes_to_retrieve,
     genes_to_retrieve,
     BM_simplified)
} else {
  message("Loading biomart data from file")
  load(biomart_data_file)
}

# Run MetaCycle ====

if (length(list.files("pre_compute/metacycle")) == 0) {
  
  message("Running JTK_cycle")
  
  # save data without ZT1 for application of metaCycle
  ZT1_names <- c("A1R1", "A1R2", "A1R3",
                 "B1R1", "B1R2", "B1R3",
                 "C1R1", "C1R2", "C1R3")
  
  log_ge_wout1 <- log_processed_ge[
    , !(colnames(log_processed_ge) %in% ZT1_names)
  ]
  
  # save separate csv files for metacycle data
  gene_abun_A <- log_ge_wout1[, 1:18]
  write.csv(gene_abun_A, "pre_compute/metacycle/gene_abun_A.csv", quote = FALSE)
  
  gene_abun_B <- log_ge_wout1[, 19:36]
  write.csv(gene_abun_B, "pre_compute/metacycle/gene_abun_B.csv", quote = FALSE)
  
  gene_abun_C <- log_ge_wout1[, 37:54]
  write.csv(gene_abun_C, "pre_compute/metacycle/gene_abun_C.csv", quote = FALSE)
  
  # run metacycle
  mc_timepoints <- rep(seq(0, 20, by = 4), each = 3)
  
  MetaCycle::meta2d(infile = "pre_compute/metacycle/gene_abun_A.csv", 
                    filestyle = "csv",
                    outdir = "pre_compute/metacycle/", 
                    timepoints = mc_timepoints,
                    cycMethod = c("JTK"), 
                    outIntegration = "noIntegration")
  
  MetaCycle::meta2d(infile = "pre_compute/metacycle/gene_abun_B.csv", 
                    filestyle = "csv",
                    outdir = "pre_compute/metacycle/", 
                    timepoints = mc_timepoints,
                    cycMethod = c("JTK"), 
                    outIntegration = "noIntegration")
  
  MetaCycle::meta2d(infile = "pre_compute/metacycle/gene_abun_C.csv", 
                    filestyle = "csv",
                    outdir = "pre_compute/metacycle/", 
                    timepoints = mc_timepoints,
                    cycMethod = c("JTK"), 
                    outIntegration = "noIntegration")
  
  rm(gene_abun_A,
     gene_abun_B,
     gene_abun_C,
     log_ge_wout1,
     mc_timepoints,
     ZT1_names)
}

# load results
message("Loading JTK_cycle results")
day_A_JTK <- read.csv("pre_compute/metacycle/JTKresult_gene_abun_A.csv")
day_B_JTK <- read.csv("pre_compute/metacycle/JTKresult_gene_abun_B.csv")
day_C_JTK <- read.csv("pre_compute/metacycle/JTKresult_gene_abun_C.csv")

# save TPM tables ====

TPM_output_dir <- "outputs/TPM_tables/"

if (length(list.files(TPM_output_dir)) == 0) {
  message("Saving gene expression tables")
  write.csv(gene_abun_ord, paste0(TPM_output_dir, "TPM_master.csv"))
  write.csv(median_df, paste0(TPM_output_dir, "TPM_medians.csv"))
  write.csv(final_gene_abun_ord, paste0(TPM_output_dir, "TPM_filtered.csv"))
}

# Garbage collect ====

gc()
