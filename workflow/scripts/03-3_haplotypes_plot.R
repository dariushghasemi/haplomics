
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--rds", default=NULL, help="Path and filename of master coloc table produced by individual traits pre-processing"),
  make_option("--annotation", default=NULL, help="Path to annotation file"),
  make_option("--variants", default=NULL, help="Path to variants list"),
  make_option("--locus", default=20, help="Region name"),
  make_option("--output1", default=NULL, help="Filename of full haplotype plot"),
  make_option("--output2", default=NULL, help="Filename of shrinked haplotype plot")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# taking variants file as input
args <- commandArgs(trailingOnly = TRUE)

# store region name
locus_name <- opt$locus


#-----------------------------------------------------#
#------                read data                ------
#-----------------------------------------------------#

# read annotated variants
df_annot <- data.table::fread(opt$annotation)

# variants file
names_bim <- c("CHROM", "POS", "ID", "REF", "ALT", "AF")
snps_list <- data.table::fread(opt$variants, header = F, col.names = names_bim)

# association results
results <- readRDS(opt$rds)


#-----------------------------------------------------#
#------           Required functions          --------
#-----------------------------------------------------#


# making dictionary to harmonize haplotype name across omics
# Function to generate a unique name for each unique combination of variants
change_haplo_name <- function(df) {
  if (all(df$Haplotype == "Ref.")) {
    # Keep "Ref." if the haplotype is identical across traits
    return("Ref.")
  } else {
  # Exclude trait_name and Haplotype columns
  uniq_haplo <- unique(df[- c(1:2)]) # df %>% select(starts_with("chr"))

  hash <- apply(uniq_haplo, 1, function(x) paste(x, collapse = "_"))
  # Use a consistent name for each group
  haplo_name <- setNames(paste0("H", seq_along(hash)), hash)

  return(haplo_name)
  }
}

#----------#
# Extract haplotypes from results and harmonizing their names
extract_haplotypes <- function(df) {
  df %>%
  ungroup() %>%
  select(trait_name, haplotype) %>%
  # take an example trait for visulization
  slice_head(n = 1) %>%
  unnest(haplotype) %>%
  select(- hap.freq) %>%
  filter(Haplotype != "Hrare") %>%
  group_by(trait_name) %>%
  mutate(
    haplo_name  = change_haplo_name(.),
    haplo_name  = if_else(Haplotype == "Ref.", "Ref.", haplo_name)
  ) %>%
  ungroup() %>%
  relocate(trait_name, Haplotype, haplo_name) %>%
  dplyr::select(- Haplotype) %>%
  rename(Haplotype = haplo_name)
}
#----------#

# harmonizing haplotype names
#mutate(Haplo2 = Haplotype, Haplotype = haplo_factor(Haplo2)) %>%

#----------#
# shaping the results for haplotypes plot
prepare_annotation <- function(df) {
  df %>%
  as_tibble() %>%
  dplyr::rename(
    Gene = gene_symbol,
    VEP_annot = most_severe_consequence
    ) %>%
  dplyr::mutate(
    snpid = str_extract(variant, "([0-9]+):([0-9]+)"),
    annot = str_replace_all(VEP_annot, "_variant", ""),
    annot = str_replace_all(annot, "3_prime_[UTR|utr]", "UTR3'"),
    annot = str_replace_all(annot, "5_prime_UTR|5_prime_utr", "UTR5'"),
    #annot = str_replace_all(annot, "missense&splice_region", "splice region")
    Gene = str_replace(Gene, "^$", "NA")
  ) %>%
    dplyr::select(snpid, variant, annot, Gene)
}

#----------#
# changing major/minor format to REF/ALT format
adding_annotation <- function(df) {  
  
  df %>%
  pivot_longer(
    cols = - c(trait_name, Haplotype),
    names_to = "SNP",
    values_to = "majorMinor"
    ) %>%
  dplyr::mutate(
    snpid = str_extract(SNP, "([0-9]+)_([0-9]+)") %>% str_replace("_", ":")
    ) %>%
  left_join(
    df_annot %>% prepare_annotation(), # add variants annotations
    join_by(snpid)
    ) %>%
  left_join(
    snps_list %>% dplyr::mutate(snpid = paste0(CHROM, ":", POS)), 
    join_by(snpid)
    ) #%>%   # add alleles frequencies
  #mutate(snp_ref_alt  = str_c(snpid, "_", REF, "/", ALT))
}

#----------#
labeling_alleles <- function(df) {
  
  df %>%
    # change the non-aligned columns' name for consistency
    rename(REF_org = REF, ALT_org = ALT, AF_org = AF) %>%
    group_by(snpid) %>%
    # in EPACTS wiki referred -> AC: Total Non-reference Allele Count
    # https://genome.sph.umich.edu/wiki/EPACTS
    mutate(
      # flipping the alleles with ALT allele frequency > 0.50
      REF = case_when(majorMinor == 2 ~ ALT_org, TRUE ~ REF_org), # & AF_org >= 0.5 
      ALT = case_when(majorMinor == 2 ~ REF_org, TRUE ~ ALT_org),
      AF  = case_when(majorMinor == 2 ~ 1 - AF_org, TRUE ~ AF_org),
      aminoAcid = case_when(
        majorMinor == 1 & AF >= 0.5 ~ ALT,
        majorMinor == 1 & AF <  0.5 ~ REF,
        majorMinor == 2 & AF <  0.5 ~ ALT,
        majorMinor == 2 & AF >= 0.5 ~ REF,
        #labeled_Allele == 1 ~ REF, labeled_Allele == 2 ~ ALT
      )
    ) %>% 
    ungroup()
    #chr_pos = str_remove_all(SNP, "[A-T-C-G]+_[A-T-C-G]+"),
}

#----------#
dying_alleles <- function(df) {
  
  df %>%
  # fixing polymorphic alleles color across haplotypes
  group_by(trait_name, snpid) %>%
  mutate(
    N_haplo    = n_distinct(Haplotype),
    N_allele   = n_distinct(aminoAcid),
    diallelic  = if_else(N_allele == 1, "", aminoAcid),
    SNP        = str_replace(SNP, "chr_", ""),
	  SNP        = str_replace(SNP, "_", ":")
  ) %>%
  ungroup()
}

#----------#
shrinking_haplotype <- function(df) {

  df %>%
    # check if a variant is biallelic
    dplyr::mutate(is_snp = nchar(REF) == 1 & nchar(ALT) == 1) %>%
    dplyr::filter(
      N_allele == 2,  # shrinking the plot to colored variants
      is_snp == TRUE  # keep biallelic variants and remove multi-allelic
      ) %>%
    ungroup()
}

#----------#
adding_reference <- function(df) {

  df_ref <- df %>%
  filter(Haplotype == "Ref.") %>% 
  select(Haplotype, snpid, majorMinor) 

  df %>% left_join(df_ref, by = "snpid", suffix = c("", "_ref"), relationship = "many-to-many")
}

#----------#
labeling_axises <- function(df) {
  
  df %>% mutate(
    # adding gene and function of the vcariants as a second x-axis
    xlab_annot  = factor(paste0(Gene, "_", annot), levels = unique(paste0(Gene, "_", annot))),
    #tlab_snpid  = factor(snp_ref_alt, levels = unique(snp_ref_alt), ordered = TRUE)
    tlab_snpid  = factor(SNP, levels = unique(SNP), ordered = TRUE)
)
}

#----------#
fading_alleles <- function(df) {

  df %>% mutate(
    # alleles differing from the alleles in "Reference" haplotypes
    differs_ref = majorMinor != majorMinor_ref | is.na(majorMinor_ref),
    true_allele = ifelse(Haplotype == "Ref." | differs_ref, aminoAcid, "")
  )
}

#----------#
# showing only variants varying across the haplotypes
haplo_plot <- function(df) {
  
  df %>%
    ggplot(aes(xlab_annot, Haplotype)) +
    #geom_point(aes(color = diallelic), size = 4.5, alpha = .75, show.legend = F) +
    geom_text(aes(label = true_allele), color = "grey20", size = 3, vjust = .45) +
    geom_hline(yintercept = num_haplo - .5, lty = 1, linewidth = .7, color = "grey50") +
    scale_color_manual(values = c("deepskyblue1", "green1", "magenta1", "#FF3434", "gold1", "grey50", "navyblue", "orange2")) +
    facet_wrap(~ tlab_snpid, scales = "free_x", nrow = 1) +
    scale_y_discrete(position = "right") +
    labs(x = "", y = "") +
    theme_classic() +
    theme(
      axis.text.x = element_text(face = "bold", angle = -90, vjust = 0.5),
          axis.text.y = element_text(face = "bold"),
          axis.title = element_blank(),
          panel.border = element_blank(), 
          panel.spacing.x = unit(0, "line"),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text.x = element_text(size = 8, face = "bold", angle = 90, vjust = 0.2, hjust = 0.0),
          # save more space for x-axis labels
          plot.margin = margin(l = 5, r = 10, t = 2, b = 2, unit = "mm")
      )
}


#-----------------------------------------------------#
#------             Haplotypes plot             ------
#-----------------------------------------------------#

# data for haplotypes plot
data_hap_plt <- results %>% 
  extract_haplotypes() %>% 
  adding_annotation() %>% 
  labeling_alleles() %>%
  dying_alleles() %>%
  adding_reference() %>% 
  fading_alleles() %>% 
  labeling_axises()

#----------#
# width and height of the plot, also for shrinked plot
num_haplo <- data_hap_plt %>% distinct(Haplotype) %>% nrow()
num_snps  <- data_hap_plt %>% distinct(snpid) %>% nrow()
num_snps_shr <- data_hap_plt %>% shrinking_haplotype() %>% distinct(snpid) %>% nrow()
num_snps_plt <- data_hap_plt %>% dplyr::filter(N_allele == 2) %>% distinct(snpid) %>% nrow()

cat(
  "No. haplotypes:", num_haplo,
  "\nNo. SNPs:", num_snps,
  "\nNo. varied SNPs:", num_snps_plt, 
  "\nNo. varied, bi-allelic variants shown on the plot:", num_snps_shr,
  "\n\n"
  )


#----------#
# shrinked haplotypes plot
shr_plt <- data_hap_plt %>% shrinking_haplotype() %>% haplo_plot()

# draw haplotypes plot
hap_plt <- data_hap_plt %>% haplo_plot()

#----------#
# taking only an example haplotype for characterizing haplo variants
#data_hap_plt %>% filter(Haplotype == "Ref.") %>%
  # dropping joined columns
  #select(- Allele_ref, - Haplotype_ref) %>%
  # creating SNPs ID
  #SNPid = str_c("chr", str_remove(SNP2, "_[A-Z]+.*"))  
  # haplotypes variants characteristics
  #write.csv("20-Nov-23_characteristics_of_haplotype_variants.csv", row.names = F)

#----------#
# save haplotypes plot
ggsave(hap_plt, filename = opt$output1, width = num_snps     / 5 + 0.5, height = num_haplo / 2 + 1.5, dpi = 300, units = "in", limitsize = FALSE)
ggsave(shr_plt, filename = opt$output2, width = num_snps_shr / 3 + 0.5, height = num_haplo / 2 + 4.5, dpi = 300, units = "in", limitsize = FALSE) #num_snps_shr / 2 - 0.5

#----------#
# print time and date
Sys.time()

quit()





#----------#
# shaping the results for haplotypes plot

results_plot <- results %>%
  # changing major/minor format to REF/ALT format
  adding_annotation() %>%
  # bounding to eGFR
  filter(trait_name == "eGFRw") %>%
  # harmonizing haplotype names
  #mutate(Haplo2 = Haplotype, Haplotype = haplo_factor(Haplo2)) %>%
  labeling_alleles() %>%
  # fixing polymorphic alleles color across haplotypes
  group_by(trait_name, snpid) %>%
  mutate(
    N_haplo    = n_distinct(Haplotype),
    N_allele   = n_distinct(aminoAcid),
    diallelic  = if_else(N_allele == 1, "", aminoAcid),
    SNP        = str_replace(SNP, "chr", ""),
	  SNP        = str_replace(SNP, "_", ":"),
    snp_ref_alt  = str_c(snpid, "_", REF, "/", ALT),
  ) %>%
  # shrinking the plot to colored variants
  #filter(N_allele == 2) %>% 
  ungroup() %>%
  # pivot_wider(id_cols = c(trait_name, Haplotype), names_from = SNP, values_from = Allele) %>%
  # mutate(Haplo = str_c(trait_name, "_", Haplotype),
  #        Haplo = str_replace_all(Haplo, "Haplo_", "H")) %>%
  # adding gene and function of the vcariants as a second x-axis
  mutate(annot2  = factor(paste0(Gene, "_", annot), levels = unique(paste0(Gene, "_", annot))),
         snpid2  = factor(snp_ref_alt, levels = unique(snp_ref_alt), ordered = TRUE))



#----------#
ggsave(hap_plt, filename = out.plot, width = num_snps / 3 - 0.5, height = num_haplo + 1.5, dpi = 350, units = "in", limitsize = FALSE)

#sbatch --wrap 'Rscript 03-3_haplotypes_plot.R output/result_associations/IGF1R_haplotypes_association.RDS data/annotation/IGF1R_annotation.txt data/annotation/IGF1R_variants.list' -c 2 --mem-per-cpu=32GB -J "03-3_IGF1R.R"