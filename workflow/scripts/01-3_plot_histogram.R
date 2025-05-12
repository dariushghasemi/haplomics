#!/usr/bin/Rscript

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

# print time and date
Sys.time()

#----------#
# Get arguments specified in the sbatch
option_list <- list(
  make_option("--variants", default=NULL, help="List of variants"),
  make_option("--output", default=NULL, help="Output filename")
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


#----------#
# define the columns names
headers <- c("CHROM", "POS", "ID", "REF", "ALT", "AF")

# read the variants file
df_variants <- data.table::fread(opt$variants, header = FALSE, sep = "\t", col.names = headers)

# making allele freq consistent for whole variant to have AF<0.5
df_variants <- df_variants %>% dplyr::mutate(AFc = if_else(AF <= 0.5, AF, 1 - AF)) 

#----------#
# set dynamic breaks for axices
n <- length(df_variants$AFc)
bin_width <- 2 * IQR(df_variants$AFc) / (n^(1/3))
bins <- ceiling((max(df_variants$AFc) - min(df_variants$AFc)) / bin_width)
x_breaks <- pretty(range(df_variants$AFc), n = 7)

# Calculate the histogram counts to determine y-axis breaks
hist_counts <- hist(df_variants$AFc, plot = FALSE, breaks = x_breaks)$counts
y_breaks <- pretty(range(hist_counts), n = 5)
y_lim <- max(y_breaks)

#----------#
# draw the plot
af_histo <- df_variants %>%
  ggplot(aes(AFc)) +
  geom_histogram(bins = 30, color = "skyblue", fill = "steelblue2") +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks, limits = c(0, y_lim)) + 
  labs(
    x = "\nAllele frequency of variants",
    y = "Number of instances"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14, face = 2),
    axis.text  = element_text(size = 14, face = 1)
  )

#----------#
# saving the histogram
ggsave(af_histo, filename = opt$output, width = 9, height = 5.5, dpi = 300, units = "in")
