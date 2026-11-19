
#-----------------------------------------------------#
#-------                Logging              ---------
#-----------------------------------------------------#

# Get log path from Snakemake, fallback if missing
log_file <- tryCatch(snakemake@log[[1]], error = function(e) "logs/plot_heatmaps/default.log")

# Ensure log directory exists
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

# Redirect stdout & stderr to the log file
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")   # redirect stdout
sink(log_con, type = "message")  # redirect messages / stderr

# Print time & message each time
log_message <- function(...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(timestamp, " - ", ...)
}

#-----------------------------------------------------#
#------              Load Libraries            -------
#-----------------------------------------------------#

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

# Params
irds <- snakemake@input[["RDS"]]
oplt <- snakemake@output[["plt"]]
#ords <- snakemake@output[["rds"]]

#-----------------------------------------------------#
#------               Import Data              -------
#-----------------------------------------------------#

log_message("Importing RDS results...")

res <- readRDS(irds)

# model coefficients
res_tidy <- res %>% 
  select(trait_name, tidy) %>%
  unnest(tidy)

#----------#

# convert result to wide format
results_wide <- res_tidy %>%
  dplyr::filter(!str_detect(term, "(Intercept)")) %>%
  filter(term != "Ref.") %>%
  dplyr::select(trait_name, term, estimate) %>%
  pivot_wider(
    names_from = trait_name,
    values_from = estimate
    ) %>%
  # convert to matrix
  column_to_rownames("term") %>%
  as.matrix.data.frame()


show_dendo <- nrow(results_wide) > 1

#-----------------------------------------------------#
#------                  Heatmap               -------
#-----------------------------------------------------#

log_message("Creating heatmap...")

png(oplt, units = "in", res = 300, width = 19, height = 6)

pheatmap(
  results_wide,
  #color = hcl.colors(50, "Blue-Red 2"),
  cluster_cols = show_dendo,
  cluster_rows = show_dendo,
  clustering_method = "ward.D2",
  na_col = "white",
  border_color = NA,
  angle_col = "270"
  )

dev.off()

log_message("Stored heatmap here: ", oplt)
