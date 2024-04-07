#!/usr/bin/Rscript

#----------#
# print time and date
Sys.time()

# date
today.date <- format(Sys.Date(), "%d-%b-%y")

#----------#
# taking variants file as input
args <- commandArgs(trailingOnly = TRUE)
variants_file <- args[1]
hist_plt_file <- args[2]

# taking the locus name
locus_name  <- gsub("_variants.list", "", basename(variants_file))
locus_name

#------------------------#
library(dplyr)
library(ggplot2)

#------------------------#
# define the columns names
variants_header <- c("CHROM", "POS", "ID", "REF", "ALT", "AF")

# read the varinats_file
df_variants <- data.table::fread(variants_file, header = FALSE, sep = "\t", col.names = variants_header)

#------------------------#
af_histo <- df_variants %>%
  mutate(AF_compl = if_else(AF <= 0.5, AF, 1 - AF)) %>% # making allele freq consistent for whole variant to have AF<0.5
  ggplot(aes(AF_compl)) +
  geom_histogram(bins = 50, color = NA, fill = "skyblue") +
  labs(
	x = paste0("Allele frequency of variants at *", locus_name, "*"),
	y = "Number of instances"
  ) +
  theme_classic() +
  theme(
  axis.title = ggtext::element_markdown(size = 14, face = 2),
  axis.text  = element_text(size = 14, face = 1)
  )



#------------------------#
# saving the histogram
ggsave(af_histo, filename = hist_plt_file, width = 9, height = 5.5, dpi = 300, units = "in")
