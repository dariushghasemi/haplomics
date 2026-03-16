
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(data.table)
})

#----------#
# Get arguments specified in the sbatch
option_list <- list(
  make_option("--rds", default=NULL, help="Path to haplo.GLM results"),
  make_option("--annotation", default=NULL, help="Path to annotation file"),
  make_option("--variants", default=NULL, help="Path to variants list"),
  make_option("--char", default=NULL, help="Output haplotype characteristics"),
  make_option("--plotfull", default=NULL, help="Filename of full haplotype plot"),
  make_option("--plottiny", default=NULL, help="Filename of shrinked haplotype plot")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

Sys.time()

#-----------------------------------------------------#
#------                read data                ------
#-----------------------------------------------------#

# read annotated variants
df_annot <- data.table::fread(opt$annotation)

# variants file
names_bim <- c("CHROM", "POS", "ID", "REF", "ALT", "AF")
snps_list <- data.table::fread(opt$variants, header = F, col.names = names_bim)

# association results
rds <- readRDS(opt$rds)


#-----------------------------------------------------#
#------           Required functions          --------
#-----------------------------------------------------#

# Rename haplotypes
rename_haplo <- function(df) {

  huniq <- df %>% select(starts_with("chr")) #omit trait_name & Haplotype
  hash  <- apply(huniq, 1, function(x) paste(x, collapse = "_"))
  hname <- setNames(paste0("H", seq_along(hash)), hash)
  return(hname)
  }

#----------#
# Extract haplotypes from results and harmonizing their names
extract_haplotypes <- function(df) {
  df %>%
  select(trait_name, haplotype) %>%
  slice_head(n = 1) %>% # take an example trait
  unnest(haplotype) %>%
  select(- trait_name) %>%
  filter(Haplotype != "geno.rare") %>%
  mutate(
    haplo_name = rename_haplo(.),
    haplo_name = if_else(Haplotype == "haplo.base", "Ref.", haplo_name),
  ) %>%
  dplyr::select(haplo_name, hap.freq, starts_with("chr")) 
}

#----------#
# shaping the results for haplotypes plot
prepare_annotation <- function(df) {
  mutate(
    df,
    snpid = str_extract(variant, "([0-9]+):([0-9]+)"),
    annot = str_replace_all(most_severe_consequence, "_variant", ""),
    annot = str_replace_all(annot, "3_prime_[UTR|utr]", "UTR3'"),
    annot = str_replace_all(annot, "5_prime_UTR|5_prime_utr", "UTR5'"),
    gene = str_replace(gene_symbol, "^$", "NA")
  ) %>%
    dplyr::select(snpid, variant, gene, annot)
}

#----------#
# changing major/minor format to REF/ALT format
adding_annotation <- function(df) {  
  
  pivot_longer(
    df,
    cols = - c(haplo_name, hap.freq),
    names_to = "SNP",
    values_to = "majorMinor"
    ) %>%
  mutate(
    chr = str_split_fixed(SNP, "_", 4)[,1],
    pos = str_split_fixed(SNP, "_", 4)[,2],
    a1  = str_split_fixed(SNP, "_", 4)[,3],
    a2  = str_split_fixed(SNP, "_", 4)[,4],
    snpid = str_c(str_remove(chr, "chr"), ":", pos),
    ID = str_c(str_remove(chr, "chr"), ":", pos, ":", a2, ":", a1)
  ) %>%       # add variants annotations
  left_join(
    df_annot %>% prepare_annotation(),
    join_by(snpid)
  ) %>%       # add alleles frequencies
  left_join(
    snps_list %>%
      mutate(SNPID = str_c(str_remove(CHROM, "chr"), POS, ALT, REF, sep = ":")) %>%
      dplyr::select(SNPID, REF, ALT, AF),
    join_by(ID == SNPID)
    )
}

#----------#
# Recode minor/major allele to actual allele
labeling_alleles <- function(df) {
  mutate(
    df,
    major_allele = if_else(AF >= 0.5, ALT, REF),
    minor_allele = if_else(AF >= 0.5, REF, ALT),
    allele = case_when(
      majorMinor == 1 ~ minor_allele,
      majorMinor == 2 ~ major_allele,
      TRUE ~ NA_character_
      ),
    # check if a variant is biallelic
    is_snp = nchar(REF) == 1 & nchar(ALT) == 1
    ) %>%
    dplyr::select(haplo_name, hap.freq, SNP, ID, gene, annot, majorMinor, allele, is_snp)
  }


#----------#
dying_alleles <- function(df) {
  
  mutate(
    df,
    N_haplo    = n_distinct(haplo_name),
    N_allele   = n_distinct(allele),
    diallelic  = if_else(N_allele == 1, "", allele),
    .by = ID # fixing polymorphic alleles color across haplotypes
  )
}

#----------#
shrinking_haplotype <- function(df) {
  
  dplyr::filter(
    df,
    N_allele == 2,  # shrinking the plot to colored variants
    is_snp == TRUE  # keep biallelic variants and remove multi-allelic
  )
}

#----------#
# To identify tagging variants, need to first
# add corresponding allele for each variant in Ref.
# haplotype and finally compare each alleles separately.

fading_alleles <- function(df) {

  df_ref <- df %>%
    filter(haplo_name == "Ref.") %>% 
    select(haplo_name, ID, majorMinor) 

  df %>% left_join(
    df_ref, by = "ID",
    suffix = c("", "_ref"),
    relationship = "many-to-many"
    ) %>%
    # alleles differing from reference alleles
    mutate(
      differs_ref = majorMinor != majorMinor_ref | is.na(majorMinor_ref),
      true_allele = ifelse(haplo_name == "Ref." | differs_ref, allele, "")
      )
  }

#----------#
labeling_axises <- function(df) {
  
  mutate(
    df,
    # adding gene and function of the vcariants as a second x-axis
    xlab_annot = factor(paste0(gene, "_", annot), levels = unique(paste0(gene, "_", annot))),
    tlab_snpid = factor(ID, levels = unique(ID), ordered = TRUE),
    # adding haplo frequency to y-axis labels
    hap_labeled = paste0(haplo_name, ", ", round(hap.freq*100, 1), "%"),
    hap_sorted  = fct_reorder(hap_labeled, hap.freq)
    )
}


#----------#
# showing only variants varying across the haplotypes
haplo_plot <- function(df) {
  
  ggplot(df, aes(xlab_annot, hap_sorted)) +
    geom_text(aes(label = true_allele), color = "grey20", size = 3, vjust = .45) +
    geom_hline(yintercept = num_haplo - .5, lty = 1, linewidth = .7, color = "grey50") +
    scale_color_manual(values = c("deepskyblue1", "green1", "magenta1", "#FF3434", "gold1", "grey50", "navyblue", "orange2")) +
    facet_wrap(~ tlab_snpid, scales = "free_x", nrow = 1) +
    scale_y_discrete(position = "right") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 8, angle = -90, vjust = 0.5),
      axis.text.y = element_text(size = 8),
      axis.title = element_blank(),
      panel.border = element_blank(),
      panel.spacing.x = unit(0, "line"),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
      # save more space for x-axis labels
      plot.margin = margin(l = 5, r = 10, t = 2, b = 2, unit = "mm")
      )
}


#-----------------------------------------------------#
#------             Haplotypes plot             ------
#-----------------------------------------------------#

# data for haplotypes plot
data_hap_plt <- rds %>% 
  extract_haplotypes() %>% 
  adding_annotation() %>% 
  labeling_alleles() %>%
  dying_alleles() %>%
  fading_alleles() %>% 
  labeling_axises()

#----------#
# width and height of the plot, also for shrinked plot
num_haplo <- data_hap_plt %>% distinct(haplo_name) %>% nrow()
num_snps  <- data_hap_plt %>% distinct(ID) %>% nrow()
num_snps_shr <- data_hap_plt %>% shrinking_haplotype() %>% distinct(ID) %>% nrow()
num_snps_plt <- data_hap_plt %>% dplyr::filter(N_allele == 2) %>% distinct(ID) %>% nrow()

message(
  "\nNo. haplotypes: ", num_haplo,
  "\nNo. SNPs: ", num_snps,
  "\nNo. varied SNPs: ", num_snps_plt, 
  "\nNo. varied, bi-allelic variants shown on the plot: ", num_snps_shr,
  "\n"
  )


#----------#
# shrinked haplotypes plot
if(num_snps_shr > 0){
  shr_plt <- data_hap_plt %>% shrinking_haplotype() %>% haplo_plot()
  } else {
    # draw full plot if no variant permuted across haplotypes
    shr_plt <- data_hap_plt %>% haplo_plot()
}

# draw haplotypes plot
hap_plt <- data_hap_plt %>% haplo_plot()

#----------#
# Haplotype variants characteristics for HTML report
data_hap_plt %>% 
  pivot_wider(
    id_cols = haplo_name:hap.freq,
    names_from = ID,
    values_from = allele
  ) %>%
  write.csv(opt$char, row.names = F)

#----------#
# save haplotypes plot
ggsave(hap_plt, filename = opt$plotfull, width = num_snps     / 5 + 0.5, height = num_haplo / 2 + 3.0, dpi = 300, units = "in", limitsize = FALSE)
ggsave(shr_plt, filename = opt$plottiny, width = num_snps_shr / 3 + 0.5, height = num_haplo / 2 + 0.0, dpi = 300, units = "in", limitsize = FALSE) #num_snps_shr / 2 - 0.5

#----------#
# print time and date
Sys.time()

