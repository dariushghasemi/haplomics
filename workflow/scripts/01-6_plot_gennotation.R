
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--annotation", default=NULL, help="Path to file with annotated variants"),
  make_option("--locus",  default=NULL, help="Locus name"),
  make_option("--output", default=NULL, help="Output filename")
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);



# store region name
locus_name <- opt$locus

# read annotated variants
df_annot <- data.table::fread(opt$annotation)


#-----------------------------------------------------#
#-------               Functions             ---------
#-----------------------------------------------------#


# shaping the results for haplotypes plot
prepare_annotation <- function(df) {
  df %>%
    as_tibble() %>%
    dplyr::rename(Gene = gene_symbol, VEP_annot = most_severe_consequence) %>%
    dplyr::mutate(
      annot = str_replace_all(VEP_annot, "_variant", ""),
      annot = str_replace_all(annot, "3_prime_UTR|3_prime_utr", "UTR3'"),
      annot = str_replace_all(annot, "5_prime_UTR|5_prime_utr", "UTR5'"),
      #annot = str_replace_all(annot, "missense&splice_region", "splice region")
    )
}
#------------------------#
# preparing frequency table of the annotations
summarize_annotation <- function(df) {
  
  df %>%
    mutate(annot = str_replace_all(annot, "_", " ")) %>%
    count(annot, Gene, name = "count_SNP") %>%
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
    theme(
      axis.title = element_text(size = 12, face = 2),
      axis.text  = element_text(size = 10, face = 2),
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = element_blank()
      )
}

bar_palette <- c(
  "#5b9aa0", "#4040a1", "#e06377", "#E69F00", "#36486b", "#a2836e", "turquoise2", "gold2",
  "#999999", "#622569",  "steelblue2", "hotpink", "springgreen2","#838060", "#82b74b", "#e0876a"
                 )


draw_annot <- function(df) {
  
  num_snps <- sum(df$count_SNP)
  y_lim <- ceiling(max(df$count_SNP) / 10) * 10
  my_breaks <- pretty(seq(0, y_lim, y_lim/10), n = 10)
  
  ggplot(df, aes(x = count_SNP, y = Gene, fill = annot_ord)) +
    geom_bar(
      stat = "identity",
      #position = position_dodge(),
      #mapping = aes(x = , y = n),
      show.legend = TRUE,
      width = 0.6,
      #fill = "grey40",
      color = NA #"steelblue2"
    ) +
    #jcolors::scale_fill_jcolors(palette = "pal8") +
    #scale_fill_brewer(palette = "Dark2") +
    scale_fill_manual(values = bar_palette, na.value = "white") +
    scale_x_continuous(breaks = my_breaks, limits = c(0, y_lim)) +
    labs(
      y = NULL, 
      x = paste0("\nNumber of variants"), 
      fill = "Annotation"
      ) +
    #coord_flip() +
    theme_linedraw() +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 10, face = 2), 
      axis.title  = element_text(size = 13, face = 2),
      axis.text.x = element_text(size = 12, face = 1),      
      axis.text.y = element_text(size = 12, face = 4, angle = 0, vjust = 0.5, hjust = 0),
      #panel.grid.major.y = element_line(linetype = 'solid', color = "grey70", linewidth = .2),
      plot.margin = margin(l = 20, r = 50, t = 5, b = 5, unit = "mm")
    )
  
}

#-----------------------------------------------------#
#-------         Gene annotation Plot        ---------
#-----------------------------------------------------#

# create plot data
plt_annot_df <- df_annot %>% 
  prepare_annotation() %>%
  #plot_annotation()
  summarize_annotation()


# draw bar plot
plt_annot <- plt_annot_df %>% draw_annot()

#------------------------#
# saving plot
ggsave(plt_annot, filename = opt$output, width = 12, height = 8, dpi = 300, units = "in")

