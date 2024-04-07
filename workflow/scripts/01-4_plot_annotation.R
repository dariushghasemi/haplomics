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
  mutate(annot = str_replace_all(annot, "_", " ")) %>%
  count(annot, name = "count_SNP") %>%
  mutate(
    percentage  = round(count_SNP/sum(count_SNP) * 100, 0),
    des_stat    = paste0(count_SNP, " (%", percentage, ")"), 
    annot_ord   = reorder(annot, desc(count_SNP))
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

draw_annot <- function(df) {
  
  num_snps <- sum(df$count_SNP)
  y_lim <- ceiling(max(df$count_SNP) / 10) * 10
  my_breaks <- round(seq(0, y_lim, length.out = 10), 0)
  
  ggplot(df, aes(x = annot_ord, y = count_SNP), fill = "grey50") +
    geom_bar(
		stat = "identity",
        position = position_dodge(),
        #mapping = aes(x = , y = n),
        show.legend = FALSE,
        width = 0.6,
        fill = "grey40",
        color = NA #"steelblue2"
    ) +
	geom_text(
		aes(label = count_SNP),
		nudge_y = 2.5, 
		fontface = 2, 
		angle = 0, 
		size = 4, 
		color = "grey40", 
		vjust = 0, 
		hjust = .5
	) +
    #jcolors::scale_fill_jcolors(palette = "pal10") +
    scale_y_continuous(breaks = my_breaks, limits = c(0, y_lim)) +
    labs(x = NULL, y = paste0("Number of instances \n at *", locus_name,"*")) + # (Total = ", num_snps, ")\n
    #coord_flip() +
    theme_classic() +
    theme(
      axis.title  = ggtext::element_markdown(size = 13, face = 2),
      axis.text.y = element_text(size = 10, face = 1),      
      axis.text.x = element_text(size = 12, face = 2, angle = -90, vjust = .5, hjust = 0),
      #panel.grid.major.y = element_line(linetype = 'solid', color = "grey70", linewidth = .2),
      plot.margin = margin(l = 2, r = 5, t = 10, b = 5, unit = "mm")
    )
  
}

#------------------------#
# create plot data
plt_annot_df <- df_annot %>% 
  prepare_annotation() %>%
  #plot_annotation()
  summarize_annotation()

plt_annot_df

# draw bar plot
plt_annot <- plt_annot_df %>% draw_annot()

#------------------------#
# saving the histogram
ggsave(plt_annot, filename = oplot_file, width = 7, height = 7, dpi = 300, units = "in")

#----------#
# print time and date
cat("\nAnnotation plot was drawn for", locus_name, " - done!\n")
Sys.time()
