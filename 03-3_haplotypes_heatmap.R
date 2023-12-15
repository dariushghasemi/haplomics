#!/usr/bin/Rscript


# The script created on 15/12/2023
# to depic the consequences of the
# variants at each replicated loci.

#----------#
# print time and date
Sys.time()

# date
today.date <- format(Sys.Date(), "%d-%b-%y")

#----------#
# taking variants file as input
args <- commandArgs(trailingOnly = TRUE)

rds_file <- args[1]

# taking the locus name
locus_name  <- gsub("\\d{2}-\\w{3}-\\d{2}_|_association_results_.+.RDS", "", basename(rds_file))
locus_name

#------------#
# directories
base.dir <- "/home/dghasemisemeskandeh/projects/haploAnalysis/output/result_associations"
results.full <- paste0(base.dir, "/", today.date, "_", locus_name, "_association_results_full1.RDS")
out.plot <- paste0(base.dir, "/output/plot_annotations/", today.date, "_", locus_name, "_plot_annotations.png")

#------------#
# function to install uninstalled required packages
is.installed <- function(package_name){
    is.element(package_name, installed.packages()[,1])
}

# check if package "Haplo.stats" is installed
if (!is.installed("pheatmap")){
   install.packages("pheatmap");
}

#------------#
library("tidyverse")
library("pheatmap")

#------------#
# preparing results for drawing heatmap
res_to_heat <- function(df){
  
  df %>%
    select(- haplotype) %>%
    unnest(tidy) %>%
    ungroup() %>%
    mutate(associated = ifelse(p.value <= 0.05, "Yes", "No"),
           term = str_replace(term,
                              "(?<=\\.)\\d{1,2}(?!\\d)",
                              sprintf("%03d", as.numeric(str_extract(term, "(?<=\\.)\\d{1,2}(?!\\d)")))),
           term = str_replace(term, "haplo_genotype.", "H")) %>%
    filter(!str_detect(term, "(Intercept)|PC|Sex|Age|rare"))
}

#----------#
# Reading and manipulating the association results for illustartion
# 01: Blood biomarkers, 02: Proteins, 03: Metabolites
#results_heatmap <- 

readRDS(rds_file) %>% #select(- haplotype) %>%
    unnest(tidy) #res_to_heat()