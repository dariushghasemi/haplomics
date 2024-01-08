#!/usr/bin/Rscript


# The script was created on 14 December 2023
# to automate visualizing reconstructed haplotypes
# for all of the 11 replicated kidney loci
# in two versions: full and shrinked haplotypes.

#----------#

library(tidyverse)

#----------#
# print time and date
Sys.time()

# date
today.date <- format(Sys.Date(), "%d-%b-%y")

#----------#
# taking variants file as input
args <- commandArgs(trailingOnly = TRUE)
rds_file <- args[1]
vep_file <- args[2]
snp_file <- args[3]

# taking locus name
#locus_name  <- gsub("\\d{2}-\\w{3}-\\d{2}_|_association_results_.+.RDS", "", basename(rds_file))
locus_name  <- gsub("\\d{2}-\\w{3}-\\d{2}_|_haplotypes_association.RDS", "", basename(rds_file))
locus_name

#----------#
# directories
out.dir <- "/home/dghasemisemeskandeh/projects/haploAnalysis/output/plot_haplotypes/"
out.plot1 <- paste0(out.dir, locus_name, "_plot_haplotypes.png") #today.date, "_", 
out.plot2 <- paste0(out.dir, locus_name, "_plot_haplotypes_shrinked.png") #today.date, "_", 


#-----------------------------------------------------#
#------                read data                ------
#-----------------------------------------------------#

# annotation data
annotation <- read.delim(
   vep_file,
   col.names  = c("CHROM", "POS", "POS37", "Gene", "Gene_ID", "VEP_annot")
   #colClasses = c("integer", "integer", "character")
)

# variants file
snps_list <- read.delim(
   snp_file,
   col.names  = c("CHROM", "POS", "POS37", "REF", "ALT", "AF")
)

# association results
results <- readRDS(rds_file)
results

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
  select(- Haplotype) %>%
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
  mutate(snpid = str_c(CHROM, ":", POS)) %>%
  mutate(
    annot = str_replace_all(VEP_annot, "_variant", ""),
    annot = str_replace_all(annot, "3_prime_[UTR|utr]", "UTR3'"),
    annot = str_replace_all(annot, "5_prime_UTR|5_prime_utr", "UTR5'"),
    #annot = str_replace_all(annot, "missense&splice_region", "splice region")
  )
}

#----------#
# changing major/minor format to REF/ALT format
adding_annotation <- function(df) {  
  
  df %>%
  pivot_longer(cols = - c(trait_name, Haplotype),
               names_to = "SNP",
               values_to = "majorMinor") %>% #labeled_Allele
  mutate(snpid = str_replace(SNP, "chr", ""),
         snpid = str_replace(snpid, "\\.", ":")) %>% #"chr\\d{1,2}."
  # add alelles frequencies
  inner_join(annotation %>% prepare_annotation(), by = "snpid", relationship = "many-to-many")  %>%
  # add variants annotations
  inner_join(snps_list %>% mutate(snpid = str_c(CHROM, ":", POS)), join_by(CHROM, POS, POS37, snpid), relationship = "many-to-many") %>%
  mutate(snp_ref_alt  = str_c(snpid, "_", REF, "/", ALT)) %>%
  select(- c(CHROM, POS, POS37))
}

#----------#
labeling_alleles <- function(df) {
  
  df %>%
    # change the non-aligned columns' name for consitency with the rest of the script
    rename(REF_org = REF, ALT_org = ALT, AF_org = AF) %>% #SNP = SNP_ID, 
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
    SNP        = str_replace(SNP, "chr", ""),
	  SNP        = str_replace(SNP, "_", ":")
  ) %>%
  ungroup()
}

#----------#
shrinking_haplotype <- function(df) {

  df %>%
  # shrinking the plot to colored variants
  filter(N_allele == 2) %>% 
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
    tlab_snpid  = factor(snp_ref_alt, levels = unique(snp_ref_alt), ordered = TRUE)
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
    #geom_hline(yintercept = 21 - .5, lty = 1, linewidth = .7, color = "grey50") +
    #scale_x_discrete(expand = c(5, 2)) +
    labs(x = "", y = "") +
    theme_classic() +
    theme(axis.text.x = element_text(face = "bold", angle = -35, hjust = 0.01),
          axis.text.y = element_text(face = "bold"),
          axis.title = element_blank(),
          panel.border = element_blank(), 
          panel.spacing.x = unit(0, "line"),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text.x = element_text(size = 8, face = "bold", angle = 90, vjust = 0.2, hjust = 0.0),
          # save more space for x-axis labels
          plot.margin = margin(l = 5, r = 20, t = 2, b = 2, unit = "mm"))
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
num_snps_shr  <- data_hap_plt %>% shrinking_haplotype() %>% distinct(snpid) %>% nrow()

cat("No. haplotypes:", num_haplo,
    "\nNo. SNPs:", num_snps,
    "\nNo. varied SNPs:", num_snps_shr, "\n\n")
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
ggsave(hap_plt, filename = out.plot1, width = num_snps     / 5 + 0.5, height = num_haplo / 2 + 1.5, dpi = 350, units = "in", limitsize = FALSE)
ggsave(shr_plt, filename = out.plot2, width = num_snps_shr / 3 + 0.5, height = num_haplo / 2 - 1.5, dpi = 350, units = "in", limitsize = FALSE) #num_snps_shr / 2 - 0.5

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