#!/usr/bin/Rscript


# taking variants file as input
args <- commandArgs(trailingOnly = TRUE)

variants_file <- args[1]

# taking the locus name
locus_name  <- gsub("_variants.list", "", basename(variants_file))

#------------------------#
library(dplyr)
library(ggplot2)

#------------------------#
# define the columns names
variants_header = c("CHROM", "POS", "ID", "REF", "ALT", "AF")

# read the varinats_file
df_variants <- data.table::fread(variants_file, header = FALSE, sep = "\t", col.names = variants_header)

#------------------------#
df_variants %>%
  # making allele freq consistent for whole variant to have AF<0.5
  mutate(AF_compl = if_else(AF <= 0.5, AF, 1 - AF)) %>%
  ggplot(aes(AF_compl)) +
  #geom_hline(yintercept = c(5), lty = 2) +
  geom_histogram(bins =50, color = "darkblue", fill = "skyblue") +
  labs(x = paste("Allele Frequency of variants at", locus_name), y = "Number of Variants") +
  theme_classic()


#------------------------#
# saving the histogram
ggsave(paste0("output/26-Oct-23_plot_histo_", locus_name ,".png"), #, 
       last_plot(), width = 9, height = 5.5, dpi = 300, units = "in")
