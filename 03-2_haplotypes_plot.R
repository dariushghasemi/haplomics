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

# taking the locus name
locus_name  <- gsub("\\d{2}-\\w{3}-\\d{2}_|_association_results_.+.RDS", "", basename(rds_file))
locus_name

#----------#
# directories
base.dir <- "/home/dghasemisemeskandeh/projects/haploAnalysis"
out.plot1 <- paste0(base.dir, "/output/", today.date, "_", locus_name, "_plot_haplotypes.png")
out.plot2 <- paste0(base.dir, "/output/", today.date, "_", locus_name, "_plot_haplotypes_shrinked.png")


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


#-----------------------------------------------------#
#------           Required functions          --------
#-----------------------------------------------------#

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
  as_tibble() %>%
  select(- tidy) %>%
  unnest(haplotype) %>%
  ungroup() %>% #distinct(Haplotype)
  filter(Haplotype != "Hrare") %>%
  pivot_longer(cols = - c(trait_name, Haplotype, hap.freq),
               names_to = "SNP",
               values_to = "majorMinor") %>% #labeled_Allele
  mutate(snpid = str_replace(SNP, "chr", ""),
         snpid = str_replace(snpid, "\\.", ":")) %>% #"chr\\d{1,2}."
  # add alelles frequencies
  inner_join(annotation %>% prepare_annotation(), by = "snpid")  %>%
  # add variants annotations
  inner_join(snps_list %>% mutate(snpid = str_c(CHROM, ":", POS)), join_by(CHROM, POS, POS37, snpid)) %>%
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

  df %>% left_join(df_ref, by = "snpid", suffix = c("", "_ref"))
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
    geom_hline(yintercept = 11 - .5, lty = 1, linewidth = .7, color = "grey50") +
    scale_color_manual(values = c("deepskyblue1", "green1", "magenta1", "#FF3434", "gold1", "grey50", "navyblue", "orange2")) +
    facet_wrap(~ tlab_snpid, scales = "free_x", nrow = 1) +
    geom_hline(yintercept = 21 - .5, lty = 1, linewidth = .7, color = "grey50") +
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

#----------#
# making dictionary to harmonize haplotype name across omics
haplo_factor <- function(x) {
  
  factor(x,
         levels = c("Haplo_030", "Haplo_033", "Haplo_034", "Haplo_065",
                    "Haplo_069", "Haplo_075", "Haplo_077", "Haplo_078",
                    "Haplo_092", "Haplo_106", "Haplo_107", "Haplo_120",
                    "Haplo_124", "Haplo_152", "Haplo_155", "Haplo_174",
                    "Haplo_181", "Haplo_182", "Haplo_208", "Haplo_257", "Haplo_280",
                    "Haplo_324", "Haplo_350", "Haplo_385", "Haplo_394",
                    "Haplo_436", "Haplo_455", "Haplo_505", "Reference"),
         
         labels = c("H1", "H1", "H1", "H2",
                    "H2", "H3", "H3", "H2",
                    "H3", "H4", "H4", "H5",
                    "H4", "H6", "H6", "H7",
                    "H6", "H7", "H7", "H8", "H8",
                    "H8", "H9", "H9", "H10",
                    "H10", "H9", "H10", "Ref."
         ))
}

# harmonizing haplotype names
#mutate(Haplo2 = Haplotype, Haplotype = haplo_factor(Haplo2)) %>%

#-----------------------------------------------------#
#------             Haplotypes plot             ------
#-----------------------------------------------------#

# data for haplotypes plot
data_hap_plt <- results %>%
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
num_haplo_shr <- data_hap_plt %>% shrinking_haplotype() %>% distinct(Haplotype) %>% nrow()
num_snps_shr  <- data_hap_plt %>% shrinking_haplotype() %>% distinct(snpid) %>% nrow()

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
ggsave(hap_plt, filename = out.plot1, width = num_snps / 3 - 0.5, height = num_haplo + 1.5, dpi = 350, units = "in", limitsize = FALSE)
ggsave(shr_plt, filename = out.plot2, width = num_snps_shr / 2 - 0.5, height = num_haplo + 1.5, dpi = 350, units = "in", limitsize = FALSE)

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

