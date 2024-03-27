#!/usr/bin/Rscript


# The script created on 15/12/2023
# to depic the consequences of the
# variants at each replicated loci.

#----------#
# print time and date
Sys.time()

# date
today.date <- format(Sys.Date(), "%d-%b-%y")

#------------------------#
# taking variants file as input
args <- commandArgs(trailingOnly = TRUE)

annot_file <- args[1]
oplot_file <- args[2]

# taking the locus name
locus_name <- gsub("_annotation.txt", "", basename(annot_file))

#------------------------#

library(dplyr)
library(stringr)
library(ggplot2)
library(ggtext) #, lib.loc = '/tmp/Rtmp7iggeu/downloaded_packages'

#------------------------#
# define the columns names
annot_header <- c("CHROM", "POS", "POS37", "Gene", "Gene_ID", "VEP_annot")

# read the varinats_file
df_annot <- read.delim(annot_file, header = FALSE, sep = "\t", col.names = annot_header)

#------------------------#
# shaping the results for haplotypes plot
prepare_annotation <- function(df) {
  df %>%
  as_tibble() %>%
  mutate(
    annot = str_replace_all(VEP_annot, "_variant", ""),
    annot = str_replace_all(annot, "3_prime_[UTR|utr]", "UTR3'"),
    annot = str_replace_all(annot, "5_prime_UTR|5_prime_utr", "UTR5'"),
    #annot = str_replace_all(annot, "missense&splice_region", "splice region")
  )
}
#------------------------#
# preparing frequency table of the annotations
summarize_annotation <- function(df) {
  
  df %>% 
  count(annot, name = "count_SNP") %>%
  mutate(
    percentage  = round(count_SNP/sum(count_SNP) * 100, 0),
    des_stat    = paste0(count_SNP, " (%", percentage, ")"), 
    annot_ord   = reorder(annot, count_SNP)
  )
}
#------------------------#
# function to draw the bar plot
plot_annotation <- function(df) {

  num_snps <- nrow(df)
  
  df %>%
  summarize_annotation() %>%
  ggplot(aes(x = annot_ord, y = count_SNP)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = "darkblue", fill = "steelblue") +
  geom_text(aes(label = des_stat), nudge_y = 2.5, fontface = 2, angle = 35, size = 3.5, color = "darkgreen") +
  labs(y = paste0("No. of Variants at *", locus_name,"* Locus (Total = ", num_snps, ")")) +
  coord_flip() +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = 2),
        axis.text  = element_text(size = 10, face = 2),
        axis.title.x = ggtext::element_markdown(),
        axis.title.y = element_blank())
}

#------------------------#
# draw bar plot
plt_annot <- df_annot %>% 
  prepare_annotation() %>%
  plot_annotation()
  
#------------------------#
# saving the histogram
ggsave(plt_annot, filename = oplot_file, width = 9.5, height = 5.5, dpi = 300, units = "in")

#----------#
# print time and date
cat("\nAnnotation plot was drawn for", locus_name, "!\n")
Sys.time()
