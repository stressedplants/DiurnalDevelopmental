# Input data and provide basic plotting functions and analysis

message('Loading packages')

library(matrixStats)
library(tidyverse)

# plotting
library(gridExtra)
library(patchwork)
library(grid)
library(ggrepel)
library(plotly)
library(svglite)

# heatmaps
library(ComplexHeatmap)
library(RColorBrewer)
library(viridis)
library(clValid)

# data IO
library(tximport)
library(biomaRt)
library(readxl)
library(xlsx)

# Colours for plotting consistently =====

day_colours <- c("#648FFF", "#DC267F", "#FFB000")
names(day_colours) <- c("A", "B", "C")

# Standard size for line plots ====

# these are the standard sizes for the line plots (produced by ggplot)
# excluding the figure legend, which is added manually in SVGs
standard_pwidth <- 3.5
standard_pheight <- 2

# other values for standard line plots
standard_text_size <- 10
standard_point_scaling_factor <- 4
standard_line_scaling_factor <- 1.25

# these are for large line plots
# only use when the plot is not being tiled horizontally
large_pwidth <- 7
large_pheight <- 5
large_base_size <- 15

dark_labels_col <- "#3b3b3b"

# Functions for use throughout =========

lp_norm <- function(a, b, p = 2) {
  (sum(abs((a - b))^p))^(1/p)
} 

# Useful for creating Excel outputs from gprofiler multi queries
# This function splits a data frame with *vector values* into separate columns.
separate_columns <- function(
    df_combined,
    col_suffixes,
    cols_to_split = c("p_values", "significant", 
                      "query_sizes", "intersection_sizes")  # for gprofiler
    
){
  
  df_output <- df_combined[, !(colnames(df_combined) %in% cols_to_split)]
  
  for (col_name in cols_to_split) {
    for (col_suf_i in 1:length(col_suffixes)) {
      # get a vector of all 'col_suf_i'th elements from 'col_name' column
      new_col <- unlist(lapply(df_combined[[col_name]], function(vec_in){
        return(vec_in[col_suf_i])
      }))
      
      new_name <- paste0(col_name, "_", col_suffixes[col_suf_i])
      # add this to the output data frame
      df_output[[new_name]] <- new_col
    }
  }
  
  return(df_output)
}

# Load data from Salmon files ======

if (file.exists("pre_compute/salmon_data.Rdata")) {
  load("pre_compute/salmon_data.Rdata")
} else {
  folder_list <- list.dirs(path="./data/quants", recursive = FALSE)
  #initalise data table with P002 - just use this for rownames
  first_fold <- "./data/quants/A0R1_quant/"
  count_raw_init <- read.table(paste(first_fold, "/quant.sf", sep=""), header=TRUE)
  
  # bit of a hack - treats .d as file extension
  gene_only <- tools::file_path_sans_ext(count_raw_init$Name)
  
  # use this data frame for input to tximport, with tx2gene
  trans_to_gene <- data.frame(count_raw_init$Name, gene_only)
  colnames(trans_to_gene) <- c("TXNAME", "GENEID")
  
  # create list of files for tximport
  # this is dictionary ordered - i.e. 12 is before 4
  file_list <- sapply(folder_list, function(i){paste0(i, "/quant.sf")})
  
  # remove the prefix and suffix of the sample name
  sample_name_list <- sapply(folder_list, function(i) {
    temp_str <- str_remove(i, "./data/quants/")
    temp_str <- str_remove(temp_str, "_quant")
    return(temp_str)
  })
  
  # this helps tximport process everything correctly
  names(file_list) <- sample_name_list
  
  # save two copies of the input data - either genes only or with all transcripts
  # note that these ARE NOT ordered by time - just a dictionary ordering
  txi_data_genes <- tximport(file_list, type="salmon", tx2gene = trans_to_gene)
  
  gene_abun <- txi_data_genes$abundance
  
  # reorder in terms of TIME not dictionary order
  ordered_list_names <- c('A0R1', 'A0R2', 'A0R3',
                          'A1R1', 'A1R2', 'A1R3',
                          'A4R1', 'A4R2', 'A4R3',
                          'A8R1', 'A8R2', 'A8R3',
                          'A12R1', 'A12R2', 'A12R3',
                          'A16R1', 'A16R2', 'A16R3',
                          'A20R1', 'A20R2', 'A20R3',
                          'B0R1', 'B0R2', 'B0R3',
                          'B1R1', 'B1R2', 'B1R3',
                          'B4R1', 'B4R2', 'B4R3',
                          'B8R1', 'B8R2', 'B8R3',
                          'B12R1', 'B12R2', 'B12R3',
                          'B16R1', 'B16R2', 'B16R3',
                          'B20R1', 'B20R2', 'B20R3',
                          'C0R1', 'C0R2', 'C0R3',
                          'C1R1', 'C1R2', 'C1R3',
                          'C4R1', 'C4R2', 'C4R3',
                          'C8R1', 'C8R2', 'C8R3',
                          'C12R1', 'C12R2', 'C12R3',
                          'C16R1', 'C16R2', 'C16R3',
                          'C20R1', 'C20R2', 'C20R3')
  gene_abun_ord <- gene_abun[,ordered_list_names]
  
  save(gene_abun_ord,
       txi_data_genes,
       ordered_list_names,
       file = "pre_compute/salmon_data.Rdata")
  
  rm(count_raw_init, file_list, folder_list,
     first_fold, gene_only, trans_to_gene, gene_abun)
}

# Create log values for gene expression =======

log_gene_abun_ord <- log2(gene_abun_ord + 1)

# Plot timeseries on top of each other =========

# takes a gene name then creates a new data frame with gene expression values
# and metadata from each sample
structure_data_from_labels <- function(gene_id, use_log = FALSE, data_df = NULL) {
  
  #  if not provided with a data frame, assume we should use gene abundance
  if (is.null(data_df)) {
    if (use_log) {
      data_df <- log_gene_abun_ord
    } else {
      data_df <- gene_abun_ord
    }
  }
  
  gene_vec <- data_df[gene_id, ]
  
  plotting_df <- data.frame(matrix(nrow = length(gene_vec),
                                   ncol = 4))
  colnames(plotting_df) <- c("series", "time", "rep", "value")
  
  # used to get the time and rep out of sample names
  split_string_by_non_digits <- function(string) {
    pattern <- "[^0-9]"
    split_string <- strsplit(string, pattern)
    return(unlist(split_string)[-1])
  }
  
  for (i in 1:nrow(plotting_df)) {
    
    sample_name <- colnames(data_df)[i]
    
    # series is always the first character of string
    plotting_df[i, "series"] <- substring(sample_name, 1, 1)
    
    # get other values from name
    useful_values <- split_string_by_non_digits(sample_name)
    plotting_df[i, "time"] <- useful_values[1]
    plotting_df[i, "rep"] <- useful_values[2]
    
    # add value
    plotting_df[i, "value"] <- gene_vec[i]
  }
  
  # fix data types
  plotting_df$series <- as.factor(plotting_df$series)
  plotting_df$time <- as.numeric(plotting_df$time)
  plotting_df$rep <- as.numeric(plotting_df$series)
  
  return(plotting_df)
}

combined_days_plot <- function(
    gene_id,
    title_text = gene_id,
    use_log = TRUE,  # standard to use log ge
    scale = FALSE,  # aka z-score
    remove_legend = TRUE,  # since mostly we add the legend in later
    text_size = standard_text_size,
    point_scaling_factor = standard_point_scaling_factor,
    line_scaling_factor = standard_line_scaling_factor,
    x_extend = 0.5 # use to expand the x-axis beyond (0,24)
) {
  
  plotting_df <- structure_data_from_labels(gene_id, use_log = use_log)
  
  # make median df
  # include `.groups = "keep"` to suppress warning message
  median_by_time_df <- plotting_df %>%
    group_by(time, series) %>%
    summarize(med_value = median(value, na.rm = TRUE), .groups = "keep")
  
  # take mean and sd from the *median* and apply this to all replicates
  if (scale) {
    medians_mean <- mean(median_by_time_df$med_value)
    medians_sd <- sd(median_by_time_df$med_value)
    median_by_time_df$med_value <- (median_by_time_df$med_value - 
                                      medians_mean) / medians_sd
    plotting_df$value <- (plotting_df$value - 
                            medians_mean) / medians_sd
  }
  
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
    theme(text = element_text(size = text_size),
          plot.title = element_text(hjust = 0.5)) +
    xlab("Time (h)") +
    guides(color = guide_legend("Day"),
           linetype = guide_legend("Day"),
           shape = guide_legend("Day"))
  
  if (use_log) {
    if (scale) {
      out_plot <- out_plot + ylab("Expression (TPM, log-transformed and scaled)") 
    } else {
      out_plot <- out_plot + ylab("Expression (TPM, log-transformed)")
    }
  } else {
    if (scale) {
      out_plot <- out_plot + ylab("Expression (TPM, scaled)") 
    } else {
      out_plot <- out_plot + ylab("Expression (TPM)")
    }
  }
  
  if (remove_legend) out_plot <- out_plot + theme(legend.position = "none")
  
  return(out_plot)
}

# Create mean and sd dataframes ====

mean_sd_file <- "pre_compute/mean_sd.Rdata"

if (!file.exists(mean_sd_file)) {

  time_labels <- c('0','1','4','8','12','16','20')
  
  a_tpms_orig <- gene_abun_ord[,1:21]
  a_tpms_avg <- as.data.frame(t(apply(a_tpms_orig, 1, function(row){
    new_mat <- matrix(row, nrow = 3)
    new_means <- c(apply(new_mat, 2, mean))
    return(new_means)
  }
  )))
  names(a_tpms_avg) <- sapply(time_labels, function(i){paste0('A',i)})
  
  b_tpms_orig <- gene_abun_ord[,22:42]
  b_tpms_avg <- as.data.frame(t(apply(b_tpms_orig, 1, function(row){
    new_mat <- matrix(row, nrow = 3)
    new_means <- c(apply(new_mat, 2, mean))
    return(new_means)
  }
  )))
  names(b_tpms_avg) <- sapply(time_labels, function(i){paste0('B',i)})
  
  c_tpms_orig <- gene_abun_ord[,43:63]
  c_tpms_avg <- as.data.frame(t(apply(c_tpms_orig, 1, function(row){
    new_mat <- matrix(row, nrow = 3)
    new_means <- c(apply(new_mat, 2, mean))
    return(new_means)
  }
  )))
  names(c_tpms_avg) <- sapply(time_labels, function(i){paste0('C',i)})
  
  all_means_mat <- cbind(a_tpms_avg, b_tpms_avg, c_tpms_avg)
  all_means_mat <- as.matrix(all_means_mat)
  
  # now sd
  a_tpms_sd <- as.data.frame(t(apply(a_tpms_orig, 1, function(row){
    new_mat <- matrix(row, nrow = 3)
    new_sds <- c(apply(new_mat, 2, sd))
    return(new_sds)
  }
  )))
  names(a_tpms_sd) <- sapply(time_labels, function(i){paste0('A',i)})
  
  b_tpms_sd <- as.data.frame(t(apply(b_tpms_orig, 1, function(row){
    new_mat <- matrix(row, nrow = 3)
    new_sds <- c(apply(new_mat, 2, sd))
    return(new_sds)
  }
  )))
  names(b_tpms_sd) <- sapply(time_labels, function(i){paste0('B',i)})
  
  c_tpms_sd <- as.data.frame(t(apply(c_tpms_orig, 1, function(row){
    new_mat <- matrix(row, nrow = 3)
    new_sds <- c(apply(new_mat, 2, sd))
    return(new_sds)
  }
  )))
  names(c_tpms_sd) <- sapply(time_labels, function(i){paste0('C',i)})
  
  all_sd_mat <- cbind(a_tpms_sd, b_tpms_sd, c_tpms_sd)
  all_sd_mat <- as.matrix(all_sd_mat)
  
  rm(a_tpms_avg, a_tpms_orig, a_tpms_sd,
     b_tpms_avg, b_tpms_orig, b_tpms_sd,
     c_tpms_avg, c_tpms_orig, c_tpms_sd)
  
  save(all_sd_mat, all_means_mat, time_labels,
       file = mean_sd_file)
} else {
  load(mean_sd_file)
}

# Create median dataframe ========

median_file <- "pre_compute/median.Rdata"

if (!file.exists(median_file)) {
  message("Creating medians dataframe")
  
  median_df <- data.frame(matrix(nrow = nrow(gene_abun_ord), ncol = 21))
  
  # colnames are correct in all_means_mat - removing any reference to reps
  colnames(median_df) <- colnames(all_means_mat)
  row.names(median_df) <- row.names(gene_abun_ord)
  
  for (row_i in 1:nrow(median_df)) {
    for (col_j in 1:ncol(median_df)) {
      col_indicies <- c(1,2,3) + 3*(col_j  - 1)
      
      median_df[row_i, col_j] <- median(
        gene_abun_ord[row_i, col_indicies]
      )
    }
    
    # progress counter
    if (row_i %% 1000 == 0) message(
      paste0("Medians calculation on ", row_i, " of ", nrow(median_df))
    )
  }
  
  rm(row_i, col_j, col_indicies)
  save(median_df, file = median_file)
} else {
  message("Loading median file from file")
  load(median_file)
}

# Create function for heatmaps / many curves at once ====

log_median_df <- log2(median_df + 1)

scale_over_days <- function(gene_expr_in, # columns are genes 
                            scale_over_all = TRUE # FALSE is scale per day
                            ) {
  
  if (scale_over_all) {
    scaled_gene_expression <- scale(gene_expr_in)
  } else {
    A_scaled <- scale(gene_expr_in[
      c("A0", "A1", "A4", "A8", "A12", "A16", "A20"), 
    ])
    
    B_scaled <- scale(gene_expr_in[
      c("B0", "B1", "B4", "B8", "B12", "B16", "B20"), 
    ])
    
    C_scaled <- scale(gene_expr_in[
      c("C0", "C1", "C4", "C8", "C12", "C16", "C20"), 
    ])
    
    scaled_gene_expression <- rbind(A_scaled, B_scaled, C_scaled)
  }
  
  return(scaled_gene_expression)
}

# can also be used to get groups for GO terms 
# take tree_col from the pheatmap object!
scaled_medians_heatmap <- function(gene_vec, 
                                   scale_over_all = TRUE,
                                   horizontal = TRUE,
                                   ...) {
  
  scaled_gene_expression <- scale_over_days(
    t(log_median_df[gene_vec, ]),
    scale_over_all = scale_over_all
  )
  
  if (horizontal) {
    ComplexHeatmap::pheatmap(
      t(scaled_gene_expression),
      cluster_cols = FALSE,
      gaps_col = c(7, 14),
      show_rownames = FALSE,
      color = inferno(10),
      ...
    )
  } else {
    ComplexHeatmap::pheatmap(
      scaled_gene_expression,
      cluster_rows = FALSE,
      gaps_row = c(7, 14),
      show_colnames = FALSE,
      color = inferno(10),
      ...
    )
  }

}

scaled_timeseries_groups <- function(
    gene_vec,
    line_alpha = 0.2,
    base_size = 14,
    timeseries_linewidth = 0.5,
    mean_linewidth = 0.5,
    ...  # passed to ggplot2::theme - to customise the individual plots
) {
  
  # need times as columns for ggplot - so transpose
  expr_df <- as_tibble(scale_over_days(
    t(log_median_df[gene_vec, ]),
    scale_over_all = TRUE
  ))
  
  ZTs <- c(0, 1, 4, 8, 12, 16, 20)
  expr_df$time <- rep(ZTs, 3)
  expr_df$day <- rep(c("A", "B", "C"), each = 7)
  
  plot_one_day <- function(day_in, ...) {
    plot_data <- expr_df %>% 
      filter(day == day_in) %>%
      dplyr::select(!day) %>%  # doesn't contain useful info
      gather(key = gene, value = value, -time)
    
    # calculate the mean over all genes and all timepoints to add to y axis
    mean_over_all <- mean(plot_data$value)
    
    out_plot <- ggplot(data = plot_data, aes(x = time, y = value, 
                                             group = gene)) +
        geom_line(alpha = line_alpha, colour = day_colours[day_in], 
                  linewidth = timeseries_linewidth) + 
        geom_hline(yintercept = mean_over_all, linetype = 2, 
                   linewidth = mean_linewidth) +
        theme_classic(base_size = base_size) +
        theme(...) +
        scale_x_continuous(breaks = ZTs)
    
    return(out_plot)
  }
  
  # make sure that they are on the same scale!
  A_plot <- plot_one_day("A", ...)
  B_plot <- plot_one_day("B", ...)
  C_plot <- plot_one_day("C", ...)
  
  lims <- c(min(layer_scales(A_plot)$y$range$range, 
                layer_scales(B_plot)$y$range$range, 
                layer_scales(C_plot)$y$range$range),
            max(layer_scales(A_plot)$y$range$range, 
                layer_scales(B_plot)$y$range$range, 
                layer_scales(C_plot)$y$range$range))
  
  plot_list <- list(
    A_plot + ylim(lims), 
    B_plot + ylim(lims),
    C_plot + ylim(lims)
  )
  
  return(plot_list)

}
