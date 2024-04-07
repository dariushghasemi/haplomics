#!/usr/bin/Rscript


# The script created on 18/12/2023
# to depic the number of exonic SNPs
# in each gene at the replicated locus.

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
library(ggtext)

#------------------------#
# define the columns names
annot_header <- c("CHROM", "POS", "POS37", "Gene", "Gene_ID", "VEP_annot")

# read the varinats_file
df_annot <- read.delim(annot_file, header = FALSE, sep = "\t", col.names = annot_header)

#------------------------#
# preparing frequency table of the annotations
summarize_annotation <- function(df) {
  
  df %>% 
  count(Gene, name = "count_SNP") %>%
  mutate(
    percentage  = round(count_SNP/sum(count_SNP) * 100, 0),
    des_stat    = paste0(count_SNP, "\n(%", percentage, ")"), 
    annot_ord   = reorder(Gene, count_SNP)
  )
}
#------------------------#
# function to draw the bar plot
plot_genes <- function(df) {

  num_snps <- nrow(df)
  
  df %>%
  summarize_annotation() %>%
  ggplot(aes(x = annot_ord, y = count_SNP)) +
  geom_bar(
	stat = "identity",
	position = position_dodge(),
	width = 0.7,
	color = NA,
	fill = "steelblue"
  ) +
  geom_text(
	aes(label = des_stat), 
	nudge_y = 2.5, 
	fontface = 2, 
	angle = 0, 
	size = 5, 
	vjust = .5, 
	hjust = .2,
	color = "grey20"
  ) +
  labs(y = paste0("Number of instances in each gene at *", locus_name,"* (total = ", num_snps, ")")) +
  coord_flip() +
  theme_classic() +
  theme(
	axis.text.x  = element_text(size = 14, face = 1),
	axis.text.y  = element_text(size = 14, face = 4),
	axis.title.x = ggtext::element_markdown(size = 14, face = 2),
    axis.title.y = element_blank(),
	plot.margin = margin(l = 5, r = 5, t = 5, b = 5, unit = "mm")
  )
}

#------------------------#
# draw bar plot
plt_annot <- df_annot %>% plot_genes()
  
#------------------------#
# saving the histogram
ggsave(plt_annot, filename = oplot_file, width = 11.5, height = 6.5, dpi = 300, units = "in")

#----------#
# print time and date
cat("\nGene plot was drawn for", locus_name, " - done!\n")
Sys.time()
