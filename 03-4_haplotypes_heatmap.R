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
locus_name  <- gsub("\\d{2}-\\w{3}-\\d{2}_|_haplotypes_association.RDS", "", basename(rds_file))
locus_name

#------------#
# directories
base.dir <- "/home/dghasemisemeskandeh/projects/haploAnalysis/output"
out.plot <- paste0(base.dir, "/plot_heatmaps/", locus_name, "_plot_heatmap_haplotypes_effect.png") #today.date, "_", 

#------------#
# function to install uninstalled required packages
is.installed <- function(package_name){
    is.element(package_name, installed.packages()[,1])
}

# check if package "Haplo.stats" is installed
if (!is.installed("digest")){
   install.packages("digest");
}

#------------#
library("tidyverse")
library("pheatmap")

#------------#
# Function to generate a unique name for each unique combination of variants
change_haplo_name <- function(df) {
  if (all(df$Haplotype == "Ref.")) {
    # Keep "Ref." if the haplotype is identical across traits
    return("Ref.")
  } else {
  # Exclude trait_name and Haplotype columns
  uniq_haplo <- unique(df[- c(1:2)])

  hash <- apply(uniq_haplo, 1, function(x) paste(x, collapse = "_"))
  # Use a consistent name for each group
  haplo_name <- setNames(paste0("H", seq_along(hash)), hash)

  return(haplo_name)
  }
}

# Function to check if haplotype is identical across traits
is_identical_haplotype <- function(df, haplotype) {
  identical_rows <- df %>%
    filter(Haplotype == haplotype) %>%
    select(-trait_name, -Haplotype) %>%
    summarise_all(~all(. == first(.))) %>%
    unlist()
  
  return(all(identical_rows))
}
#----------#
problematic <- c("ALT_GPT", "AST_GOT", "DBP", "SBP", "Pulse_Rate", "FT4", "ALP")
#----------#
# saving haplotype name
#haplo_dict0 <- 
readRDS(rds_file) %>% 
  ungroup() %>%
  select(trait_name, haplotype) %>% 
  unnest(haplotype) %>% 
  select(- hap.freq) %>% 
  filter(Haplotype != "Hrare") %>% 
  #select(trait_name, Haplotype, "chr15.98649165", "chr15.98649166", "chr15.98649359", "chr15.98649374")%>%
  #slice_head(n=72) %>% 
  count(trait_name) %>% print(n=Inf)
  #filter(!trait_name %in% problematic) %>%

quit()

variants <- grep("^chr", names(haplo_dict0), value = TRUE)
variants

haplo_dict <- haplo_dict0 %>% 
  #group_by(trait_name) %>% #count(trait_name) %>% print(n = Inf)
  #change_haplo_name(.)
  mutate(haplo = do.call(paste, c(select(., all_of(variants)), sep = "_"))) %>%
  add_count(haplo, name = "haplo_count") %>%
  filter(haplo_count == max(haplo_count)) %>% #select(trait_name, Haplotype, haplo_count, haplo) %>% print(n =Inf) 
  group_by(trait_name) %>%
  mutate(
    haplo_name = change_haplo_name(.),
    haplo_name = if_else(Haplotype == "Ref.", "Ref.", haplo_name)
  ) %>%
  ungroup() %>%
  select(trait_name, Haplotype, haplo_name)
  
#haplo_dict %>% print(n = Inf)

#table(haplo_dict$trait_name, haplo_dict$Haplotype)
#quit()
#----------#
# preparing results for drawing heatmap
res_to_heat <- function(df){
  
  df %>%
    #select(- haplotype) %>%
    unnest(tidy) %>%
    ungroup() %>%
    mutate(
      #associated = ifelse(p.value <= 0.05, "Yes", "No"),
      term = str_replace(
        term, 
        "(?<=\\.)\\d{1,2}(?!\\d)",
        sprintf("%03d", as.numeric(str_extract(term, "(?<=\\.)\\d{1,2}(?!\\d)")))
        	),
      term = str_replace(term, "haplo_genotype.", "H")
      ) %>%
    filter(!str_detect(term, "(Intercept)|PC|Sex|Age|rare")) %>%
    # reshaping results for pheatmap
    select(trait_name, term, estimate) %>%
    # edited from left_join to right_join to keep only available haplotypes for all traits
    right_join(haplo_dict, by = c("trait_name" = "trait_name", "term" = "Haplotype")) %>%
    filter(term != "Ref.") %>%
    select(- term) %>%
    pivot_wider(names_from = trait_name, values_from = estimate)
}

#----------#
# Reading and manipulating the association results for illustartion
# 01: Blood biomarkers, 02: Proteins, 03: Metabolites
results_heatmap <- readRDS(rds_file) %>% res_to_heat()
results_heatmap

check_haplo <- nrow(results_heatmap) > 1
check_haplo
#----------#
# pheatmap
png(out.plot, units = "in", res = 400, width = 12, height = 6)

pheatmap(results_heatmap[-1],
         #color = hcl.colors(50, "Blue-Red 2"),
         #breaks = seq(-rg, rg, length.out = 100), #rg <- max(abs(results_heatmap[-1]))
         #color = myColor, 
         #breaks = myBreaks,
         labels_row = results_heatmap$haplo_name,
         #display_numbers = results_omics_pval[-1],
         number_color = "gold",
         cluster_cols = check_haplo,
         cluster_rows = check_haplo,
         clustering_method = "ward.D2",
         na_col = "white",
         border_color = NA, 
         #annotation_col = annot_omics,
         #annotation_colors = annot_colors,
         #display_numbers = F,
         fontsize_number = 15,
         fontsize_row = 10,
         fontsize_col = 10,
         angle_col = "270")

dev.off()

#----------#
#sbatch --wrap 'Rscript 03-4_haplotypes_heatmap.R output/result_associations/IGF1R_haplotypes_association.RDS ' -c 2 --mem-per-cpu=32GB -J "03-2_IGF1R.R"