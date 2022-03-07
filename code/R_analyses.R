#########################################################################################
# 2022
# Extraction kit testing - round 2
# Justin Shaffer
# justinparkshaffer@gmail.com
#########################################################################################
#########################################################################################

# Set working directory
#########################################################################################
getwd()
setwd("~/Mycelium/R/2022_extraction_test_round2/")


# Install and load libraries needed for analysis
#########################################################################################

install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")
install.packages("tidyverse")
install.packages("vegan")
install.packages("cluster")
install.packages("dendextend")
install.packages("indicspecies")
install.packages("gdata")
install.packages("reshape2")

library(qiime2R)
library(tidyverse)
library(vegan)
library(cluster)
library(dendextend)
library(indicspecies)
library(plyr)

# For hierarchical clustering
# NOTE: Make sure plyr is OFF
library(vegan)
library(cluster)
library(dendextend)
library(gdata)
library(reshape2)


#######################################################################################################
# Append read count and alpha-diversity vectors to sample metadata
#######################################################################################################

# Read in  metadata
md <- read_tsv("12201_metadata.txt")


# Read in read count data
read_counts_16S <- read_tsv("dna_all_16S_deblur_biom_read_counts.txt")
read_counts_its <- read_tsv("dna_all_ITS_deblur_biom_read_counts.txt")
read_counts_shotgun <- read_tsv("dna_all_shotgun_woltka_wol_biom_read_counts.txt")

# Merge read count data with metadata
md_read_counts_16S <- left_join(md, read_counts_16S, by = "sample_name")
md_read_counts_amplicon <- left_join(md_read_counts_16S, read_counts_its, by = "sample_name")
md_read_counts <- left_join(md_read_counts_amplicon, read_counts_shotgun, by = "sample_name")

# Export new metadata file
write_tsv(md_read_counts, path = "12201_metadata_round2_read_counts.txt")


# Read in alpha-diversity metrics for 16S data
alpha_16S_high_faithspd <- read_qza("dna_bothPS_16S_deblur_biom_lod_noChl_noMit_sepp_gg_hbm_rar12150_alpha_faithspd.qza")
alpha_16S_high_shannon <- read_qza("dna_bothPS_16S_deblur_biom_lod_noChl_noMit_sepp_gg_hbm_rar12150_alpha_shannon.qza")
alpha_16S_high_richness <- read_qza("dna_bothPS_16S_deblur_biom_lod_noChl_noMit_sepp_gg_hbm_rar12150_alpha_richness.qza")

alpha_16S_low_faithspd <- read_qza("dna_bothPS_16S_deblur_biom_lod_noChl_noMit_sepp_gg_lbm_rar3150_alpha_faithspd.qza")
alpha_16S_low_shannon <- read_qza("dna_bothPS_16S_deblur_biom_lod_noChl_noMit_sepp_gg_lbm_rar3150_alpha_shannon.qza")
alpha_16S_low_richness <- read_qza("dna_bothPS_16S_deblur_biom_lod_noChl_noMit_sepp_gg_lbm_rar3150_alpha_richness.qza")


# Read in alpha-diversity metrics for ITS data
alpha_its_high_fishersalpha <- read_qza("dna_bothPS_ITS_deblur_biom_lod_hbm_rar1171_alpha_fishersalpha.qza")
alpha_its_high_shannon <- read_qza("dna_bothPS_ITS_deblur_biom_lod_hbm_rar1171_alpha_shannon.qza")
alpha_its_high_richness <- read_qza("dna_bothPS_ITS_deblur_biom_lod_hbm_rar1171_alpha_richness.qza")

alpha_its_low_fishersalpha <- read_qza("dna_bothPS_ITS_deblur_biom_lod_lbm_rar456_alpha_fishersalpha.qza")
alpha_its_low_shannon <- read_qza("dna_bothPS_ITS_deblur_biom_lod_lbm_rar456_alpha_shannon.qza")
alpha_its_low_richness <- read_qza("dna_bothPS_ITS_deblur_biom_lod_lbm_rar456_alpha_richness.qza")


# Read in alpha-diversity metrics for shotgun metagenomics data
alpha_shotgun_high_faithspd <- read_qza("dna_bothPS_shotgun_woltka_wol_biom_hbm_rar16600_alpha_faithspd.qza")
alpha_shotgun_high_shannon <- read_qza("dna_bothPS_shotgun_woltka_wol_biom_hbm_rar16600_alpha_shannon.qza")
alpha_shotgun_high_richness <- read_qza("dna_bothPS_shotgun_woltka_wol_biom_hbm_rar16600_alpha_richness.qza")

alpha_shotgun_low_faithspd <- read_qza("dna_bothPS_shotgun_woltka_wol_biom_lbm_rar485_alpha_faithspd.qza")
alpha_shotgun_low_shannon <- read_qza("dna_bothPS_shotgun_woltka_wol_biom_lbm_rar485_alpha_shannon.qza")
alpha_shotgun_low_richness <- read_qza("dna_bothPS_shotgun_woltka_wol_biom_lbm_rar485_alpha_richness.qza")


# Concatenate alpha-diversity metrics from low- and high-biomass samples for 16S data
alpha_16S_faithspd <- rbind((rownames_to_column(as.data.frame(alpha_16S_high_faithspd$data), "sample_name")), (rownames_to_column(as.data.frame(alpha_16S_low_faithspd$data), "sample_name")))

alpha_16S_shannon <- rbind((rownames_to_column(as.data.frame(alpha_16S_high_shannon$data), "sample_name")), (rownames_to_column(as.data.frame(alpha_16S_low_shannon$data), "sample_name")))

alpha_16S_richness <- rbind((rownames_to_column(as.data.frame(alpha_16S_high_richness$data), "sample_name")), (rownames_to_column(as.data.frame(alpha_16S_low_richness$data), "sample_name")))


# Concatenate alpha-diversity metrics from low- and high-biomass samples for ITS data
alpha_its_fishersalpha <- rbind((rownames_to_column(as.data.frame(alpha_its_high_fishersalpha$data), "sample_name")), (rownames_to_column(as.data.frame(alpha_its_low_fishersalpha$data), "sample_name")))

alpha_its_shannon <- rbind((rownames_to_column(as.data.frame(alpha_its_high_shannon$data), "sample_name")), (rownames_to_column(as.data.frame(alpha_its_low_shannon$data), "sample_name")))

alpha_its_richness <- rbind((rownames_to_column(as.data.frame(alpha_its_high_richness$data), "sample_name")), (rownames_to_column(as.data.frame(alpha_its_low_richness$data), "sample_name")))


# Concatenate alpha-diversity metrics from low- and high-biomass samples for shotgun metagenomics data
alpha_shotgun_faithspd <- rbind((rownames_to_column(as.data.frame(alpha_shotgun_high_faithspd$data), "sample_name")), (rownames_to_column(as.data.frame(alpha_shotgun_low_faithspd$data), "sample_name")))

alpha_shotgun_shannon <- rbind((rownames_to_column(as.data.frame(alpha_shotgun_high_shannon$data), "sample_name")), (rownames_to_column(as.data.frame(alpha_shotgun_low_shannon$data), "sample_name")))

alpha_shotgun_richness <- rbind((rownames_to_column(as.data.frame(alpha_shotgun_high_richness$data), "sample_name")), (rownames_to_column(as.data.frame(alpha_shotgun_low_richness$data), "sample_name")))


# Combine alpha-diversity metrics with metadata
alpha_16S_faithspd %>%
  right_join(md_read_counts) -> md_read_counts

alpha_16S_shannon %>%
  right_join(md_read_counts) -> md_read_counts

alpha_16S_richness %>%
  right_join(md_read_counts) -> md_read_counts

colnames(md_read_counts)[2] <- "alpha_16S_rar12150_3150_richness"
colnames(md_read_counts)[3] <- "alpha_16S_rar12150_3150_shannon"
colnames(md_read_counts)[4] <- "alpha_16S_rar12150_3150_faithspd"

alpha_its_fishersalpha %>%
  right_join(md_read_counts) -> md_read_counts

alpha_its_shannon %>%
  right_join(md_read_counts) -> md_read_counts

alpha_its_richness %>%
  right_join(md_read_counts) -> md_read_counts

colnames(md_read_counts)[2] <- "alpha_its_rar1171_456_richness"
colnames(md_read_counts)[3] <- "alpha_its_rar1171_456_shannon"
colnames(md_read_counts)[4] <- "alpha_its_rar1171_456_fishersalpha"


alpha_shotgun_faithspd %>%
  right_join(md_read_counts) -> md_read_counts

alpha_shotgun_shannon %>%
  right_join(md_read_counts) -> md_read_counts

alpha_shotgun_richness %>%
  right_join(md_read_counts) -> md_read_counts

colnames(md_read_counts)[2] <- "alpha_shotgun_rar16600_485_richness"
colnames(md_read_counts)[3] <- "alpha_shotgun_rar16600_485_shannon"
colnames(md_read_counts)[4] <- "alpha_shotgun_rar16600_485_faithspd"


# Export new metadata file
write_tsv(md_read_counts, path = "12201_metadata_round2_read_counts_alpha_diversity.txt")


#######################################################################################################
# Scatter plots comparing across all samples types for each kit
#######################################################################################################
library(tidyverse)
library(ggpubr)
library(cowplot)
library(svglite)
library(scales)


# Tidy-up metadata for running correlations and plotting
#######################################################################################################

# Read in metadata
md_long <- read_tsv("12201_metadata_round2_read_counts_alpha_diversity.txt")

# Subset round 1 and round 2 samples and exclude round 1 Zymo
md_long_round2_1 <- subset(md_long, md_long$round != 3)
md_long_round2_2 <- subset(md_long_round2_1, md_long_round2_1$round !=1 |
                             md_long_round2_1$extraction_kit != "Zymo MagBead")
# Subset metadata to include DNA concentrations, read counts, and alpha-diversity metrics
names(md_long_round2_2)
md_long_dna_reads_alpha <- md_long_round2_2[,c(1:10, 12:14, 19, 25:29, 64, 83:85)]

# Create wide form data for dna concentrations
md_wide_dna_reads_alpha <- pivot_wider(md_long_dna_reads_alpha,
                                       id_cols = c(unique_sample_id,
                                                   round,
                                                   biomass_plate,
                                                   sample_type,
                                                   sample_type_1,
                                                   sample_type_2,
                                                   sample_type_3),
                                       names_from = extraction_kit,
                                       values_from = c(dna_conc_ng_ul,
                                                       `16s_deblur_reads`,
                                                       its_deblur_reads,
                                                       shotgun_woltka_reads,
                                                       alpha_shotgun_rar16600_485_richness,
                                                       alpha_shotgun_rar16600_485_shannon,
                                                       alpha_shotgun_rar16600_485_faithspd,
                                                       alpha_its_rar1171_456_richness,
                                                       alpha_its_rar1171_456_shannon,
                                                       alpha_its_rar1171_456_fishersalpha,
                                                       alpha_16S_rar12150_3150_richness,
                                                       alpha_16S_rar12150_3150_shannon,
                                                       alpha_16S_rar12150_3150_faithspd),
                                       values_fn = mean)

# Re-order levels for sample type 2
md_wide_dna_reads_alpha$sample_type_2 <- factor(md_wide_dna_reads_alpha$sample_type_2, levels = c("cat feces", "human feces", "human saliva", "mouse feces", "mouse jejunum tissue", "bare soil",  "rhizosphere soil", "kimchi", "yogurt", "keyboard", "floor", "human female urine", "human male urine", "human foot", "human armpit", "human milk", "seawater", "freshwater", "ZymoBIOMICS Microbial Community Standard I", "ZymoBIOMICS Microbial Community Standard II", "PCR extraction control"))

# Re-name levels for sample type 2
levels(md_wide_dna_reads_alpha$sample_type_2) <- c("Cat feces", "Human feces", "Human saliva", "Mouse feces", "Mouse tissue, jejunum", "Soil, bare",  "Soil, rhizosphere", "Food, kimchi", "Food, yogurt", "Surface, keyboard", "Surface, floor tile", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human milk", "Water, saline", "Water, non-saline", "Mock community I", "Mock community II", "PCR extraction control")

# Make vector for colors for sample type 2
sample_type_2_colors <- c("#F48022", "#A6CFE5", "#8AC752", "#FCBF6D", "#22201F", "#F9F49D", "#B15829", "#6A3F98", "#CAB2D5", "#8AC752", "#1B8D42", "#CAB2D5", "#6A3F98", "#F48021", "#FCC06E", "#22201F", "#1B79B5", "#A6CFE6", "gray60", "#gray60", "#22201F")

# Re-name levels for biomass_plate
md_wide_dna_reads_alpha$biomass_plate <- as.factor(md_wide_dna_reads_alpha$biomass_plate)

# Re-name levels for biomass_plate type 2
levels(md_wide_dna_reads_alpha$biomass_plate) <- c("High biomass", "Low biomass", "not provided")

# Subset for round 1 vs. round 2
md_wide_round1 <- subset(md_wide_dna_reads_alpha, md_wide_dna_reads_alpha$round == 1)
md_wide_round2 <- subset(md_wide_dna_reads_alpha, md_wide_dna_reads_alpha$round == 2)

# Exclude controls
md_wide_round1_clean <- subset(md_wide_round1, md_wide_round1$sample_type != 'control blank' &
                                 md_wide_round1$sample_type != 'mock community')
md_wide_round2_clean <- subset(md_wide_round2, md_wide_round2$sample_type != 'control blank' &
                                 md_wide_round2$sample_type != 'mock community')


# Plot DNA concentrations
#######################################################################################################

# Exclude negative values
md_wide_round1_clean <- subset(md_wide_round1_clean, md_wide_round1_clean$dna_conc_ng_ul_PowerSoil > 0 &
                                 md_wide_round1_clean$dna_conc_ng_ul_Norgen > 0)
md_wide_round2_clean <- subset(md_wide_round2_clean, 
                               md_wide_round2_clean$`dna_conc_ng_ul_PowerSoil Pro` > 0 |
                               md_wide_round2_clean$`dna_conc_ng_ul_MagMAX Microbiome` > 0 |
                               md_wide_round2_clean$`dna_conc_ng_ul_NucleoMag Food` > 0 |
                               md_wide_round2_clean$`dna_conc_ng_ul_Zymo MagBead` > 0)

# Test for normality
ggqqplot(md_wide_round1_clean$`dna_conc_ng_ul_PowerSoil`)
ggqqplot(md_wide_round1_clean$`dna_conc_ng_ul_PowerSoil Pro`)
ggqqplot(md_wide_round1_clean$`dna_conc_ng_ul_Norgen`)
ggqqplot(md_wide_round2_clean$`dna_conc_ng_ul_MagMAX Microbiome`)
ggqqplot(md_wide_round2_clean$`dna_conc_ng_ul_NucleoMag Food`)
ggqqplot(md_wide_round2_clean$`dna_conc_ng_ul_Zymo MagBead`)
## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
powersoil_pro_dna_corr_kendall <- cor.test(x = md_wide_round1_clean$`dna_conc_ng_ul_PowerSoil`, y = md_wide_round1_clean$`dna_conc_ng_ul_PowerSoil Pro`, method = "kendall")

norgen_dna_corr_kendall <- cor.test(x = md_wide_round1_clean$`dna_conc_ng_ul_PowerSoil`, y = md_wide_round1_clean$`dna_conc_ng_ul_Norgen`, method = "kendall")

magmax_dna_corr_kendall <- cor.test(x = md_wide_round2_clean$`dna_conc_ng_ul_PowerSoil`, y = md_wide_round2_clean$`dna_conc_ng_ul_MagMAX Microbiome`, method = "kendall")

nucleomag_dna_corr_kendall <- cor.test(x = md_wide_round2_clean$`dna_conc_ng_ul_PowerSoil`, y = md_wide_round2_clean$`dna_conc_ng_ul_NucleoMag Food`, method = "kendall")

zymo_dna_corr_kendall <- cor.test(x = md_wide_round2_clean$`dna_conc_ng_ul_PowerSoil`, y = md_wide_round2_clean$`dna_conc_ng_ul_Zymo MagBead`, method = "kendall")


# View correlations
powersoil_pro_dna_corr_kendall
norgen_dna_corr_kendall
magmax_dna_corr_kendall
nucleomag_dna_corr_kendall
zymo_dna_corr_kendall


# PowerSoil Pro
dna_powersoil_pro <- ggplot(md_wide_round1_clean, aes(x = `dna_conc_ng_ul_PowerSoil`,
                                                      y = `dna_conc_ng_ul_PowerSoil Pro`)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\nDNA yield (ng/µL)") +
  ylab("PowerSoil Pro\nDNA yield (ng/µL)") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_log10() +
  scale_y_log10()


# Norgen
dna_norgen <- ggplot(md_wide_round1_clean, aes(x = `dna_conc_ng_ul_PowerSoil` + 1,
                                               y = `dna_conc_ng_ul_Norgen` + 1)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\nDNA yield (ng/µL)") +
  ylab("Norgen\nDNA yield (ng/µL)") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_log10() +
  scale_y_log10()


# MagMAX Microbiome
dna_magmax <- ggplot(md_wide_round2_clean, aes(x = `dna_conc_ng_ul_PowerSoil`,
                                               y = `dna_conc_ng_ul_MagMAX Microbiome`)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\nDNA yield (ng/µL)") +
  ylab("MagMAX Microbiome\nDNA yield (ng/µL)") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_log10() +
  scale_y_log10()


# NucleoMag
dna_nucleomag <- ggplot(md_wide_round2_clean, aes(x = `dna_conc_ng_ul_PowerSoil`,
                                                  y = `dna_conc_ng_ul_NucleoMag Food`)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\nDNA yield (ng/µL)") +
  ylab("NucleoMag Food\nDNA yield (ng/µL)") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_log10() +
  scale_y_log10()


# NucleoMag
dna_zymo <- ggplot(md_wide_round2_clean, aes(x = `dna_conc_ng_ul_PowerSoil`,
                                                  y = `dna_conc_ng_ul_Zymo MagBead`)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\nDNA yield (ng/µL)") +
  ylab("Zymo MagBead\nDNA yield (ng/µL)") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_log10() +
  scale_y_log10()


# Cowplot - all kits
figure_corr_dna <- plot_grid(dna_powersoil_pro,
                             dna_norgen,
                             dna_magmax,
                             dna_nucleomag,
                             dna_zymo,
                             labels = c("A", "B", "C", "D", "E"),
                             label_size = 24,
                             label_fontfamily = 'sans',
                             ncol = 3)

figure_corr_dna

save_plot('figure_scatter_dna_11x7.svg',
          figure_corr_dna,
          base_height = 7,
          base_width = 11)


# 16S deblur reads
#######################################################################################################

# Test for normality
ggqqplot(md_wide_round1_clean$`16s_deblur_reads_PowerSoil`)
ggqqplot(md_wide_round1_clean$`16s_deblur_reads_PowerSoil Pro`)
ggqqplot(md_wide_round1_clean$`16s_deblur_reads_Norgen`)
ggqqplot(md_wide_round2_clean$`16s_deblur_reads_MagMAX Microbiome`)
ggqqplot(md_wide_round2_clean$`16s_deblur_reads_NucleoMag Food`)
ggqqplot(md_wide_round2_clean$`16s_deblur_reads_Zymo MagBead`)
## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
powersoil_pro_reads_16S_corr_kendall <- cor.test(x = md_wide_round1_clean$`16s_deblur_reads_PowerSoil`, y = md_wide_round1_clean$`16s_deblur_reads_PowerSoil Pro`, method = "kendall")

norgen_reads_16S_corr_kendall <- cor.test(x = md_wide_round1_clean$`16s_deblur_reads_PowerSoil`, y = md_wide_round1_clean$`16s_deblur_reads_Norgen`, method = "kendall")

magmax_reads_16S_corr_kendall <- cor.test(x = md_wide_round2_clean$`16s_deblur_reads_PowerSoil`, y = md_wide_round2_clean$`16s_deblur_reads_MagMAX Microbiome`, method = "kendall")

nucleomag_reads_16S_corr_kendall <- cor.test(x = md_wide_round2_clean$`16s_deblur_reads_PowerSoil`, y = md_wide_round2_clean$`16s_deblur_reads_NucleoMag Food`, method = "kendall")

zymo_reads_16S_corr_kendall <- cor.test(x = md_wide_round2_clean$`16s_deblur_reads_PowerSoil`, y = md_wide_round2_clean$`16s_deblur_reads_Zymo MagBead`, method = "kendall")


# View correlations
powersoil_pro_reads_16S_corr_kendall
norgen_reads_16S_corr_kendall
magmax_reads_16S_corr_kendall
nucleomag_reads_16S_corr_kendall
zymo_reads_16S_corr_kendall


# Scatterplot - PowerSoil Pro
reads_16S_powersoil_pro <- ggplot(md_wide_round1_clean, aes(x = `16s_deblur_reads_PowerSoil`/1000,
                                 y = `16s_deblur_reads_PowerSoil Pro`/1000)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\n16S Deblur reads (x1000)") +
  ylab("PowerSoil Pro\n16S Deblur reads (x1000)") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type")


# Scatterplot - Norgen
reads_16S_norgen <- ggplot(md_wide_round1_clean, aes(x = `16s_deblur_reads_PowerSoil`/1000,
                                 y = `16s_deblur_reads_Norgen`/1000)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\n16S Deblur reads (x1000)") +
  ylab("Norgen\n16S Deblur reads (x1000)") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type")


# Scatterplot - MagMAX Microbiome
reads_16S_magmax <- ggplot(md_wide_round2_clean, aes(x = `16s_deblur_reads_PowerSoil`/1000,
                                 y = `16s_deblur_reads_MagMAX Microbiome`/1000)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\n16S Deblur reads (x1000)") +
  ylab("MagMAX Microbiome\n16S Deblur reads (x1000)") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type")


# Scatterplot - Nucleomag Food
reads_16S_nucleomag <- ggplot(md_wide_round2_clean, aes(x = `16s_deblur_reads_PowerSoil`/1000,
                                 y = `16s_deblur_reads_NucleoMag Food`/1000)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\n16S Deblur reads (x1000)") +
  ylab("NucleoMag Food\n16S Deblur reads (x1000)") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type")


# Scatterplot - Zymo MagBead
reads_16S_zymo <- ggplot(md_wide_round2_clean, aes(x = `16s_deblur_reads_PowerSoil`/1000,
                                 y = `16s_deblur_reads_Zymo MagBead`/1000)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\n16S Deblur reads (x1000)") +
  ylab("Zymo MagBead\n16S Deblur reads (x1000)") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type")


# Cowplot - all kits
figure_corr_reads_16S <- plot_grid(reads_16S_powersoil_pro,
                                   reads_16S_norgen,
                                   reads_16S_magmax,
                                   reads_16S_nucleomag,
                                   reads_16S_zymo,
                                   labels = c("A", "B", "C", "D", "E"),
                                   label_size = 18,
                                   label_fontfamily = 'sans',
                                   ncol = 3)

figure_corr_reads_16S

save_plot('figure_scatter_reads_16S_11x7.svg',
          figure_corr_reads_16S,
          base_height = 7,
          base_width = 11)

save_plot('figure_scatter_reads_16S_11x7.png',
          figure_corr_reads_16S,
          base_height = 7,
          base_width = 11)

save_plot('figure_scatter_reads_16S_11x7.tiff',
          figure_corr_reads_16S,
          base_height = 7,
          base_width = 11)


# ITS deblur reads
#######################################################################################################

# Test for normality
ggqqplot(md_wide_round1_clean$`its_deblur_reads_PowerSoil`)
ggqqplot(md_wide_round1_clean$`its_deblur_reads_PowerSoil Pro`)
ggqqplot(md_wide_round1_clean$`its_deblur_reads_Norgen`)
ggqqplot(md_wide_round2_clean$`its_deblur_reads_MagMAX Microbiome`)
ggqqplot(md_wide_round2_clean$`its_deblur_reads_NucleoMag Food`)
ggqqplot(md_wide_round2_clean$`its_deblur_reads_Zymo MagBead`)
## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
powersoil_pro_reads_ITS_corr_kendall <- cor.test(x = md_wide_round1_clean$`its_deblur_reads_PowerSoil`, y = md_wide_round1_clean$`its_deblur_reads_PowerSoil Pro`, method = "kendall")

norgen_reads_ITS_corr_kendall <- cor.test(x = md_wide_round1_clean$`its_deblur_reads_PowerSoil`, y = md_wide_round1_clean$`its_deblur_reads_Norgen`, method = "kendall")

magmax_reads_ITS_corr_kendall <- cor.test(x = md_wide_round2_clean$`its_deblur_reads_PowerSoil`, y = md_wide_round2_clean$`its_deblur_reads_MagMAX Microbiome`, method = "kendall")

nucleomag_reads_ITS_corr_kendall <- cor.test(x = md_wide_round2_clean$`its_deblur_reads_PowerSoil`, y = md_wide_round2_clean$`its_deblur_reads_NucleoMag Food`, method = "kendall")

zymo_reads_ITS_corr_kendall <- cor.test(x = md_wide_round2_clean$`its_deblur_reads_PowerSoil`, y = md_wide_round2_clean$`its_deblur_reads_Zymo MagBead`, method = "kendall")


# View correlations
powersoil_pro_reads_ITS_corr_kendall
norgen_reads_ITS_corr_kendall
magmax_reads_ITS_corr_kendall
nucleomag_reads_ITS_corr_kendall
zymo_reads_ITS_corr_kendall


# Scatterplot - PowerSoil Pro
reads_ITS_powersoil_pro <- ggplot(md_wide_round1_clean, aes(x = `its_deblur_reads_PowerSoil`,
                                                            y = `its_deblur_reads_PowerSoil Pro`)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\nITS Deblur reads") +
  ylab("PowerSoil Pro\nITS Deblur reads") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))


# Scatterplot - Norgen
reads_ITS_norgen <- ggplot(md_wide_round1_clean, aes(x = `its_deblur_reads_PowerSoil`,
                                                     y = `its_deblur_reads_Norgen`)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\nITS Deblur reads") +
  ylab("Norgen\nITS Deblur reads") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))


# Scatterplot - MagMAX Microbiome
reads_ITS_magmax <- ggplot(md_wide_round2_clean, aes(x = `its_deblur_reads_PowerSoil`,
                                                     y = `its_deblur_reads_MagMAX Microbiome`)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\nITS Deblur reads") +
  ylab("MagMAX Microbiome\nITS Deblur reads") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))


# Scatterplot - Nucleomag Food
reads_ITS_nucleomag <- ggplot(md_wide_round2_clean, aes(x = `its_deblur_reads_PowerSoil`,
                                                        y = `its_deblur_reads_NucleoMag Food`)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\nITS Deblur reads") +
  ylab("NucleoMag Food\nITS Deblur reads") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))


# Scatterplot - Zymo MagBead
reads_ITS_zymo <- ggplot(md_wide_round2_clean, aes(x = `its_deblur_reads_PowerSoil`,
                                                   y = `its_deblur_reads_Zymo MagBead`)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\nITS Deblur reads") +
  ylab("Zymo MagBead\nITS Deblur reads") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))


# Cowplot - all kits
figure_corr_reads_ITS <- plot_grid(reads_ITS_powersoil_pro,
                                   reads_ITS_norgen,
                                   reads_ITS_magmax,
                                   reads_ITS_nucleomag,
                                   reads_ITS_zymo,
                                   labels = c("A", "B", "C", "D", "E"),
                                   label_size = 18,
                                   label_fontfamily = 'sans',
                                   ncol = 3)

figure_corr_reads_ITS

save_plot('figure_scatter_reads_ITS_11x7.svg',
          figure_corr_reads_ITS,
          base_height = 7,
          base_width = 11)

save_plot('figure_scatter_reads_ITS_11x7.png',
          figure_corr_reads_ITS,
          base_height = 7,
          base_width = 11)

save_plot('figure_scatter_reads_ITS_11x7.tiff',
          figure_corr_reads_ITS,
          base_height = 7,
          base_width = 11)


# Plot Shotgun host- and quality-filtered reads - Scatterplots comparing the same sample between two kits
#######################################################################################################

# Test for normality
ggqqplot(md_wide_round1_clean$`shotgun_woltka_reads_PowerSoil`)
ggqqplot(md_wide_round1_clean$`shotgun_woltka_reads_PowerSoil Pro`)
ggqqplot(md_wide_round1_clean$`shotgun_woltka_reads_Norgen`)
ggqqplot(md_wide_round2_clean$`shotgun_woltka_reads_MagMAX Microbiome`)
ggqqplot(md_wide_round2_clean$`shotgun_woltka_reads_NucleoMag Food`)
ggqqplot(md_wide_round2_clean$`shotgun_woltka_reads_Zymo MagBead`)
## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
powersoil_pro_reads_shotgun_corr_kendall <- cor.test(x = md_wide_round1_clean$`shotgun_woltka_reads_PowerSoil`, y = md_wide_round1_clean$`shotgun_woltka_reads_PowerSoil Pro`, method = "kendall")

norgen_reads_shotgun_corr_kendall <- cor.test(x = md_wide_round1_clean$`shotgun_woltka_reads_PowerSoil`, y = md_wide_round1_clean$`shotgun_woltka_reads_Norgen`, method = "kendall")

magmax_reads_shotgun_corr_kendall <- cor.test(x = md_wide_round2_clean$`shotgun_woltka_reads_PowerSoil`, y = md_wide_round2_clean$`shotgun_woltka_reads_MagMAX Microbiome`, method = "kendall")

nucleomag_reads_shotgun_corr_kendall <- cor.test(x = md_wide_round2_clean$`shotgun_woltka_reads_PowerSoil`, y = md_wide_round2_clean$`shotgun_woltka_reads_NucleoMag Food`, method = "kendall")

zymo_reads_shotgun_corr_kendall <- cor.test(x = md_wide_round2_clean$`shotgun_woltka_reads_PowerSoil`, y = md_wide_round2_clean$`shotgun_woltka_reads_Zymo MagBead`, method = "kendall")


# View correlations
powersoil_pro_reads_shotgun_corr_kendall
norgen_reads_shotgun_corr_kendall
magmax_reads_shotgun_corr_kendall
nucleomag_reads_shotgun_corr_kendall
zymo_reads_shotgun_corr_kendall


# Scatterplot - PowerSoil Pro
reads_shotgun_powersoil_pro <- ggplot(md_wide_round1_clean, aes(x = `shotgun_woltka_reads_PowerSoil`,
                                                                y = `shotgun_woltka_reads_PowerSoil Pro`)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\nMetagenomic reads") +
  ylab("PowerSoil Pro\nMetagenomic reads") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))



# Scatterplot - Norgen
reads_shotgun_norgen <- ggplot(md_wide_round1_clean, aes(x = `shotgun_woltka_reads_PowerSoil`,
                                                         y = `shotgun_woltka_reads_Norgen`)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\nMetagenomic reads") +
  ylab("Norgen\nMetagenomic reads") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))


# Scatterplot - MagMAX Microbiome
reads_shotgun_magmax <- ggplot(md_wide_round2_clean, aes(x = `shotgun_woltka_reads_PowerSoil`,
                                                         y = `shotgun_woltka_reads_MagMAX Microbiome`)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\nMetagenomic reads") +
  ylab("MagMAX Microbiome\nMetagenomic reads") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))


# Scatterplot - Nucleomag Food
reads_shotgun_nucleomag <- ggplot(md_wide_round2_clean, aes(x = `shotgun_woltka_reads_PowerSoil`,
                                                            y = `shotgun_woltka_reads_NucleoMag Food`)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\nMetagenomic reads") +
  ylab("NucleoMag Food\nMetagenomic reads") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))


# Scatterplot - Zymo MagBead
reads_shotgun_zymo <- ggplot(md_wide_round2_clean, aes(x = `shotgun_woltka_reads_PowerSoil`,
                                                       y = `shotgun_woltka_reads_Zymo MagBead`)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("PowerSoil\nMetagenomic reads") +
  ylab("Zymo MagBead\nMetagenomic reads") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type_2,
                 shape = biomass_plate)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))


# Cowplot - all kits
figure_corr_reads_shotgun <- plot_grid(reads_shotgun_powersoil_pro,
                                       reads_shotgun_norgen,
                                       reads_shotgun_magmax,
                                       reads_shotgun_nucleomag,
                                       reads_shotgun_zymo,
                                       labels = c("A", "B", "C", "D", "E"),
                                       label_size = 18,
                                       label_fontfamily = 'sans',
                                       ncol = 3)

figure_corr_reads_shotgun

save_plot('figure_scatter_reads_shotgun_11x7.svg',
          figure_corr_reads_shotgun,
          base_height = 7,
          base_width = 11)

save_plot('figure_scatter_reads_shotgun_11x7.png',
          figure_corr_reads_shotgun,
          base_height = 7,
          base_width = 11)

save_plot('figure_scatter_reads_shotgun_11x7.tiff',
          figure_corr_reads_shotgun,
          base_height = 7,
          base_width = 11)


# Cowplot - all read counts
#######################################################################################################
figure_corr_reads_plot <- plot_grid(reads_16S_powersoil_pro,
                                       reads_16S_norgen,
                                       reads_16S_magmax,
                                       reads_16S_nucleomag,
                                       reads_16S_zymo,
                                       reads_ITS_powersoil_pro,
                                       reads_ITS_norgen,
                                       reads_ITS_magmax,
                                       reads_ITS_nucleomag,
                                       reads_ITS_zymo,
                                       reads_shotgun_powersoil_pro,
                                       reads_shotgun_norgen,
                                       reads_shotgun_magmax,
                                       reads_shotgun_nucleomag,
                                       reads_shotgun_zymo,
                                       labels = c("A", "B", "C", "D", "E",
                                                  "F", "G", "H", "I", "J",
                                                  "K", "L", "M", "N", "O"),
                                       label_size = 24,
                                       label_fontfamily = 'sans',
                                       ncol = 5)

figure_corr_reads_plot

save_plot('figure_scatter_reads_20x12.tiff',
          figure_corr_reads_plot,
          base_width = 20,
          base_height = 12)

save_plot('figure_scatter_reads_20x12.svg',
          figure_corr_reads_plot,
          base_width = 20,
          base_height = 12)


#######################################################################################################
# Plot Alpha-diversity
#######################################################################################################
library(tidyverse)
library(ggpubr)
library(cowplot)
library(svglite)

# Read in metadata
md_long <- read_tsv("12201_metadata_round2_read_counts_alpha_diversity.txt")


# Subset round 1 and round 2 samples and exclude round 1 Zymo
md_long_round2_1 <- subset(md_long, md_long$round != 3)
md_long_round2_2 <- subset(md_long_round2_1, md_long_round2_1$round !=1 |
                             md_long_round2_1$extraction_kit != "Zymo MagBead")

# Subset metadata to include DNA concentrations, read counts, and alpha-diversity metrics
names(md_long_round2_2)
md_long_dna_reads_alpha <- md_long_round2_2[,c(1:10, 12:14, 19, 25:29, 64, 83:85)]


# Create wide form data for dna concentrations
md_wide_dna_reads_alpha <- pivot_wider(md_long_dna_reads_alpha,
                                       id_cols = c(unique_sample_id,
                                                   round,
                                                   biomass_plate,
                                                   sample_type,
                                                   sample_type_1,
                                                   sample_type_2,
                                                   sample_type_3),
                                       names_from = extraction_kit,
                                       values_from = c(dna_conc_ng_ul,
                                                       `16s_deblur_reads`,
                                                       its_deblur_reads,
                                                       shotgun_woltka_reads,
                                                       alpha_shotgun_rar16600_485_richness,
                                                       alpha_shotgun_rar16600_485_shannon,
                                                       alpha_shotgun_rar16600_485_faithspd,
                                                       alpha_its_rar1171_456_richness,
                                                       alpha_its_rar1171_456_shannon,
                                                       alpha_its_rar1171_456_fishersalpha,
                                                       alpha_16S_rar12150_3150_richness,
                                                       alpha_16S_rar12150_3150_shannon,
                                                       alpha_16S_rar12150_3150_faithspd),
                                       values_fn = mean)


# Re-order levels for sample type 2
md_wide_dna_reads_alpha$sample_type_2 <- factor(md_wide_dna_reads_alpha$sample_type_2, levels = c("cat feces", "human feces", "human saliva", "mouse feces", "mouse jejunum tissue", "bare soil",  "rhizosphere soil", "kimchi", "yogurt", "keyboard", "floor", "human female urine", "human male urine", "human foot", "human armpit", "human milk", "seawater", "freshwater", "ZymoBIOMICS Microbial Community Standard I", "ZymoBIOMICS Microbial Community Standard II", "PCR extraction control"))


# Re-name levels for sample type 2
levels(md_wide_dna_reads_alpha$sample_type_2) <- c("Cat feces", "Human feces", "Human saliva", "Mouse feces", "Mouse tissue, jejunum", "Soil, bare",  "Soil, rhizosphere", "Food, kimchi", "Food, yogurt", "Surface, keyboard", "Surface, floor tile", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human milk", "Water, saline", "Water, non-saline", "Mock community I", "Mock community II", "PCR extraction control")


# Make vector for colors for sample type 2
sample_type_2_colors <- c("#F48022", "#A6CFE5", "#8AC752", "#FCBF6D", "#22201F", "#F9F49D", "#B15829", "#6A3F98", "#CAB2D5", "#8AC752", "#1B8D42", "#CAB2D5", "#6A3F98", "#F48021", "#FCC06E", "#22201F", "#1B79B5", "#A6CFE6", "gray60", "#gray60", "#22201F")


# Re-name levels for biomass_plate
md_wide_dna_reads_alpha$biomass_plate <- as.factor(md_wide_dna_reads_alpha$biomass_plate)


# Re-name levels for biomass_plate type 2
levels(md_wide_dna_reads_alpha$biomass_plate) <- c("High biomass", "Low biomass", "not provided")


# Subset for round 1 vs. round 2
md_wide_round1 <- subset(md_wide_dna_reads_alpha, md_wide_dna_reads_alpha$round == 1)
md_wide_round2 <- subset(md_wide_dna_reads_alpha, md_wide_dna_reads_alpha$round == 2)


# Exclude controls
md_wide_round1_clean <- subset(md_wide_round1, md_wide_round1$sample_type != 'control blank' &
                                 md_wide_round1$sample_type != 'mock community')
md_wide_round2_clean <- subset(md_wide_round2, md_wide_round2$sample_type != 'control blank' &
                                 md_wide_round2$sample_type != 'mock community')


# 16S
#######################################################################################################

# Test for normality
ggqqplot(md_wide_round1_clean$`alpha_16S_rar12150_3150_faithspd_PowerSoil`)
ggqqplot(md_wide_round1_clean$`alpha_16S_rar12150_3150_faithspd_PowerSoil Pro`)
ggqqplot(md_wide_round1_clean$`alpha_16S_rar12150_3150_faithspd_Norgen`)
ggqqplot(md_wide_round2_clean$`alpha_16S_rar12150_3150_faithspd_MagMAX Microbiome`)
ggqqplot(md_wide_round2_clean$`alpha_16S_rar12150_3150_faithspd_NucleoMag Food`)
ggqqplot(md_wide_round2_clean$`alpha_16S_rar12150_3150_faithspd_Zymo MagBead`)
## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
powersoil_pro_alpha_16S_corr_kendall <- cor.test(x = md_wide_round1_clean$`alpha_16S_rar12150_3150_faithspd_PowerSoil`, y = md_wide_round1_clean$`alpha_16S_rar12150_3150_faithspd_PowerSoil Pro`, method = "kendall")

norgen_alpha_16S_corr_kendall <- cor.test(x = md_wide_round1_clean$`alpha_16S_rar12150_3150_faithspd_PowerSoil`, y = md_wide_round1_clean$`alpha_16S_rar12150_3150_faithspd_Norgen`, method = "kendall")

magmax_alpha_16S_corr_kendall <- cor.test(x = md_wide_round2_clean$`alpha_16S_rar12150_3150_faithspd_PowerSoil`, y = md_wide_round2_clean$`alpha_16S_rar12150_3150_faithspd_MagMAX Microbiome`, method = "kendall")

nucleomag_alpha_16S_corr_kendall <- cor.test(x = md_wide_round2_clean$`alpha_16S_rar12150_3150_faithspd_PowerSoil`, y = md_wide_round2_clean$`alpha_16S_rar12150_3150_faithspd_NucleoMag Food`, method = "kendall")

zymo_alpha_16S_corr_kendall <- cor.test(x = md_wide_round2_clean$`alpha_16S_rar12150_3150_faithspd_PowerSoil`, y = md_wide_round2_clean$`alpha_16S_rar12150_3150_faithspd_Zymo MagBead`, method = "kendall")


# View correlations
powersoil_pro_alpha_16S_corr_kendall
norgen_alpha_16S_corr_kendall
magmax_alpha_16S_corr_kendall
nucleomag_alpha_16S_corr_kendall
zymo_alpha_16S_corr_kendall


# Scatterplot - PowerSoil Pro - with legend for adding to combined figure post-production
ggplot(md_wide_round1_clean, aes(x = `alpha_16S_rar12150_3150_faithspd_PowerSoil`,
                                 y = `alpha_16S_rar12150_3150_faithspd_PowerSoil Pro`)) +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\n16S Faith's PD") +
  ylab("PowerSoil Pro\n16S Faith's PD") +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,21,21,21,21,21,
                                                     22,22,22,22,22,22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type")


# Scatterplot - PowerSoil Pro
alpha_16S_powersoil_pro <- ggplot(md_wide_round1_clean, aes(x = `alpha_16S_rar12150_3150_faithspd_PowerSoil`,
                                 y = `alpha_16S_rar12150_3150_faithspd_PowerSoil Pro`)) +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\n16S Faith's PD") +
  ylab("PowerSoil Pro\n16S Faith's PD") +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22))


# Scatterplot - Norgen
alpha_16S_norgen <- ggplot(md_wide_round1_clean, aes(x = `alpha_16S_rar12150_3150_faithspd_PowerSoil`,
                                 y = `alpha_16S_rar12150_3150_faithspd_Norgen`)) +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\n16S Faith's PD") +
  ylab("Norgen\n16S Faith's PD") +

  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22))


# Scatterplot - MagMAX Microbiome
alpha_16S_magmax <- ggplot(md_wide_round2_clean, aes(x = `alpha_16S_rar12150_3150_faithspd_PowerSoil`,
                                 y = `alpha_16S_rar12150_3150_faithspd_MagMAX Microbiome`)) +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\n16S Faith's PD") +
  ylab("MagMAX Microbiome\n16S Faith's PD") +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22))


# Scatterplot - Nucleomag Food
alpha_16S_nucleomag <- ggplot(md_wide_round2_clean, aes(x = `alpha_16S_rar12150_3150_faithspd_PowerSoil`,
                                 y = `alpha_16S_rar12150_3150_faithspd_NucleoMag Food`)) +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\n16S Faith's PD") +
  ylab("NucleoMag Food\n16S Faith's PD") +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22))


# Scatterplot - Zymo MagBead
alpha_16S_zymo <- ggplot(md_wide_round2_clean, aes(x = `alpha_16S_rar12150_3150_faithspd_PowerSoil`,
                                 y = `alpha_16S_rar12150_3150_faithspd_Zymo MagBead`)) +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\n16S Faith's PD") +
  ylab("Zymo MagBead\n16S Faith's PD") +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22))


# Cowplot - all kits
figure_corr_alpha_16S <- plot_grid(alpha_16S_powersoil_pro,
                                   alpha_16S_norgen,
                                   alpha_16S_magmax,
                                   alpha_16S_nucleomag,
                                   alpha_16S_zymo,
                                   labels = c("A", "B", "C", "D", "E"),
                                   label_size = 18,
                                   label_fontfamily = 'sans',
                                   ncol = 3)

figure_corr_alpha_16S

save_plot('figure_scatter_alpha_16S_11x7.svg',
          figure_corr_alpha_16S,
          base_height = 7,
          base_width = 11)

save_plot('figure_scatter_alpha_16S_11x7.png',
          figure_corr_alpha_16S,
          base_height = 7,
          base_width = 11)

save_plot('figure_scatter_alpha_16S_11x7.tiff',
          figure_corr_alpha_16S,
          base_height = 7,
          base_width = 11)


# ITS
#######################################################################################################

# Test for normality
ggqqplot(md_wide_round1_clean$`alpha_its_rar1171_456_fishersalpha_PowerSoil`)
ggqqplot(md_wide_round1_clean$`alpha_its_rar1171_456_fishersalpha_PowerSoil Pro`)
ggqqplot(md_wide_round1_clean$`alpha_its_rar1171_456_fishersalpha_Norgen`)
ggqqplot(md_wide_round2_clean$`alpha_its_rar1171_456_fishersalpha_MagMAX Microbiome`)
ggqqplot(md_wide_round2_clean$`alpha_its_rar1171_456_fishersalpha_NucleoMag Food`)
ggqqplot(md_wide_round2_clean$`alpha_its_rar1171_456_fishersalpha_Zymo MagBead`)
## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
powersoil_pro_alpha_its_corr_kendall <- cor.test(x = md_wide_round1_clean$`alpha_its_rar1171_456_fishersalpha_PowerSoil`, y = md_wide_round1_clean$`alpha_its_rar1171_456_fishersalpha_PowerSoil Pro`, method = "kendall")

norgen_alpha_its_corr_kendall <- cor.test(x = md_wide_round1_clean$`alpha_its_rar1171_456_fishersalpha_PowerSoil`, y = md_wide_round1_clean$`alpha_its_rar1171_456_fishersalpha_Norgen`, method = "kendall")

magmax_alpha_its_corr_kendall <- cor.test(x = md_wide_round2_clean$`alpha_its_rar1171_456_fishersalpha_PowerSoil`, y = md_wide_round2_clean$`alpha_its_rar1171_456_fishersalpha_MagMAX Microbiome`, method = "kendall")

nucleomag_alpha_its_corr_kendall <- cor.test(x = md_wide_round2_clean$`alpha_its_rar1171_456_fishersalpha_PowerSoil`, y = md_wide_round2_clean$`alpha_its_rar1171_456_fishersalpha_NucleoMag Food`, method = "kendall")

zymo_alpha_its_corr_kendall <- cor.test(x = md_wide_round2_clean$`alpha_its_rar1171_456_fishersalpha_PowerSoil`, y = md_wide_round2_clean$`alpha_its_rar1171_456_fishersalpha_Zymo MagBead`, method = "kendall")


# View correlations
powersoil_pro_alpha_its_corr_kendall
norgen_alpha_its_corr_kendall
magmax_alpha_its_corr_kendall
nucleomag_alpha_its_corr_kendall
zymo_alpha_its_corr_kendall


# Scatterplot - PowerSoil Pro
alpha_its_powersoil_pro <- ggplot(md_wide_round1_clean, aes(x = `alpha_its_rar1171_456_fishersalpha_PowerSoil`,
                                                            y = `alpha_its_rar1171_456_fishersalpha_PowerSoil Pro`)) +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\nITS Fisher's alpha") +
  ylab("PowerSoil Pro\nITS Fisher's alpha") +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22))


# Scatterplot - Norgen
alpha_its_norgen <- ggplot(md_wide_round1_clean, aes(x = `alpha_its_rar1171_456_fishersalpha_PowerSoil`,
                                                     y = `alpha_its_rar1171_456_fishersalpha_Norgen`)) +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\nITS Fisher's alpha") +
  ylab("Norgen\nITS Fisher's alpha") +
  
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22))


# Scatterplot - MagMAX Microbiome
alpha_its_magmax <- ggplot(md_wide_round2_clean, aes(x = `alpha_its_rar1171_456_fishersalpha_PowerSoil`,
                                                     y = `alpha_its_rar1171_456_fishersalpha_MagMAX Microbiome`)) +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\nITS Fisher's alpha") +
  ylab("MagMAX Microbiome\nITS Fisher's alpha") +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22))


# Scatterplot - Nucleomag Food
alpha_its_nucleomag <- ggplot(md_wide_round2_clean, aes(x = `alpha_its_rar1171_456_fishersalpha_PowerSoil`,
                                                        y = `alpha_its_rar1171_456_fishersalpha_NucleoMag Food`)) +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\nITS Fisher's alpha") +
  ylab("NucleoMag Food\nITS Fisher's alpha") +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22))


# Scatterplot - Zymo MagBead
alpha_its_zymo <- ggplot(md_wide_round2_clean, aes(x = `alpha_its_rar1171_456_fishersalpha_PowerSoil`,
                                                   y = `alpha_its_rar1171_456_fishersalpha_Zymo MagBead`)) +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\nITS Fisher's alpha") +
  ylab("Zymo MagBead\nITS Fisher's alpha") +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22))


# Cowplot - all kits
figure_corr_alpha_its <- plot_grid(alpha_its_powersoil_pro,
                                   alpha_its_norgen,
                                   alpha_its_magmax,
                                   alpha_its_nucleomag,
                                   alpha_its_zymo,
                                   labels = c("A", "B", "C", "D", "E"),
                                   label_size = 18,
                                   label_fontfamily = 'sans',
                                   ncol = 3)

figure_corr_alpha_its

save_plot('figure_scatter_alpha_its_11x7.svg',
          figure_corr_alpha_its,
          base_height = 7,
          base_width = 11)

save_plot('figure_scatter_alpha_its_11x7.png',
          figure_corr_alpha_its,
          base_height = 7,
          base_width = 11)

save_plot('figure_scatter_alpha_its_11x7.tiff',
          figure_corr_alpha_its,
          base_height = 7,
          base_width = 11)


# Shotgun
#######################################################################################################

# Test for normality
ggqqplot(md_wide_round1_clean$`alpha_shotgun_rar16600_485_faithspd_PowerSoil`)
ggqqplot(md_wide_round1_clean$`alpha_shotgun_rar16600_485_faithspd_PowerSoil Pro`)
ggqqplot(md_wide_round1_clean$`alpha_shotgun_rar16600_485_faithspd_Norgen`)
ggqqplot(md_wide_round2_clean$`alpha_shotgun_rar16600_485_faithspd_MagMAX Microbiome`)
ggqqplot(md_wide_round2_clean$`alpha_shotgun_rar16600_485_faithspd_NucleoMag Food`)
ggqqplot(md_wide_round2_clean$`alpha_shotgun_rar16600_485_faithspd_Zymo MagBead`)
## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
powersoil_pro_alpha_shotgun_corr_kendall <- cor.test(x = md_wide_round1_clean$`alpha_shotgun_rar16600_485_faithspd_PowerSoil`, y = md_wide_round1_clean$`alpha_shotgun_rar16600_485_faithspd_PowerSoil Pro`, method = "kendall")

norgen_alpha_shotgun_corr_kendall <- cor.test(x = md_wide_round1_clean$`alpha_shotgun_rar16600_485_faithspd_PowerSoil`, y = md_wide_round1_clean$`alpha_shotgun_rar16600_485_faithspd_Norgen`, method = "kendall")

magmax_alpha_shotgun_corr_kendall <- cor.test(x = md_wide_round2_clean$`alpha_shotgun_rar16600_485_faithspd_PowerSoil`, y = md_wide_round2_clean$`alpha_shotgun_rar16600_485_faithspd_MagMAX Microbiome`, method = "kendall")

nucleomag_alpha_shotgun_corr_kendall <- cor.test(x = md_wide_round2_clean$`alpha_shotgun_rar16600_485_faithspd_PowerSoil`, y = md_wide_round2_clean$`alpha_shotgun_rar16600_485_faithspd_NucleoMag Food`, method = "kendall")

zymo_alpha_shotgun_corr_kendall <- cor.test(x = md_wide_round2_clean$`alpha_shotgun_rar16600_485_faithspd_PowerSoil`, y = md_wide_round2_clean$`alpha_shotgun_rar16600_485_faithspd_Zymo MagBead`, method = "kendall")


# View correlations
powersoil_pro_alpha_shotgun_corr_kendall
norgen_alpha_shotgun_corr_kendall
magmax_alpha_shotgun_corr_kendall
nucleomag_alpha_shotgun_corr_kendall
zymo_alpha_shotgun_corr_kendall


# Scatterplot - PowerSoil Pro
alpha_shotgun_powersoil_pro <- ggplot(md_wide_round1_clean, aes(x = `alpha_shotgun_rar16600_485_faithspd_PowerSoil`,
                                                                y = `alpha_shotgun_rar16600_485_faithspd_PowerSoil Pro`)) +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\nShotgun Faith's PD") +
  ylab("PowerSoil Pro\nShotgun Faith's PD") +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22))


# Scatterplot - Norgen
alpha_shotgun_norgen <- ggplot(md_wide_round1_clean, aes(x = `alpha_shotgun_rar16600_485_faithspd_PowerSoil`,
                                                         y = `alpha_shotgun_rar16600_485_faithspd_Norgen`)) +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\nShotgun Faith's PD") +
  ylab("Norgen\nShotgun Faith's PD") +
  
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22))


# Scatterplot - MagMAX Microbiome
alpha_shotgun_magmax <- ggplot(md_wide_round2_clean, aes(x = `alpha_shotgun_rar16600_485_faithspd_PowerSoil`,
                                                         y = `alpha_shotgun_rar16600_485_faithspd_MagMAX Microbiome`)) +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\nShotgun Faith's PD") +
  ylab("MagMAX Microbiome\nShotgun Faith's PD") +
  xlim(c(0,80)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22))


# Scatterplot - Nucleomag Food
alpha_shotgun_nucleomag <- ggplot(md_wide_round2_clean, aes(x = `alpha_shotgun_rar16600_485_faithspd_PowerSoil`,
                                                            y = `alpha_shotgun_rar16600_485_faithspd_NucleoMag Food`)) +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\nShotgun Faith's PD") +
  ylab("NucleoMag Food\nShotgun Faith's PD") +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22))


# Scatterplot - Zymo MagBead
alpha_shotgun_zymo <- ggplot(md_wide_round2_clean, aes(x = `alpha_shotgun_rar16600_485_faithspd_PowerSoil`,
                                                       y = `alpha_shotgun_rar16600_485_faithspd_Zymo MagBead`)) +
  geom_smooth(color = "gray60",
              size = 0.5,
              method = "lm") +
  geom_abline(linetype = "dashed",
              color = "gray60",
              size = 0.5) +
  geom_point(aes(shape = biomass_plate,
                 fill = sample_type_2),
             size = 2) +
  xlab("PowerSoil\nShotgun Faith's PD") +
  ylab("Zymo MagBead\nShotgun Faith's PD") +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        legend.position = "none") +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22))


# Cowplot - all shotgun
figure_corr_alpha_shotgun <- plot_grid(alpha_shotgun_powersoil_pro,
                                       alpha_shotgun_norgen,
                                       alpha_shotgun_magmax,
                                       alpha_shotgun_nucleomag,
                                       alpha_shotgun_zymo,
                                       labels = c("A", "B", "C", "D", "E"),
                                       label_size = 18,
                                       label_fontfamily = 'sans',
                                       ncol = 3)

figure_corr_alpha_shotgun

save_plot('figure_scatter_alpha_shotgun_11x7.svg',
          figure_corr_alpha_shotgun,
          base_height = 7,
          base_width = 11)

save_plot('figure_scatter_alpha_shotgun_11x7.png',
          figure_corr_alpha_shotgun,
          base_height = 7,
          base_width = 11)

save_plot('figure_scatter_alpha_shotgun_11x7.tiff',
          figure_corr_alpha_shotgun,
          base_height = 7,
          base_width = 11)


# Cowplot - all read counts
#######################################################################################################
figure_corr_alpha_plot <- plot_grid(alpha_16S_powersoil_pro,
                                    alpha_16S_norgen,
                                    alpha_16S_magmax,
                                    alpha_16S_nucleomag,
                                    alpha_16S_zymo,
                                    alpha_its_powersoil_pro,
                                    alpha_its_norgen,
                                    alpha_its_magmax,
                                    alpha_its_nucleomag,
                                    alpha_its_zymo,
                                    alpha_shotgun_powersoil_pro,
                                    alpha_shotgun_norgen,
                                    alpha_shotgun_magmax,
                                    alpha_shotgun_nucleomag,
                                    alpha_shotgun_zymo,
                                    labels = c("A", "B", "C", "D", "E",
                                               "F", "G", "H", "I", "J",
                                               "K", "L", "M", "N", "O"),
                                    label_size = 24,
                                    label_fontfamily = 'sans',
                                    ncol = 5)

figure_corr_alpha_plot

save_plot('figure_scatter_alpha_20x12.tiff',
          figure_corr_alpha_plot,
          base_width = 20,
          base_height = 12)

save_plot('figure_scatter_alpha_20x12.svg',
          figure_corr_alpha_plot,
          base_width = 20,
          base_height = 12)


#######################################################################################################
# Boxplots
#######################################################################################################

# Tidy-up metadata for plotting
#######################################################################################################

# Read in sample metadata
md_long <- read_tsv("12201_metadata_round2_read_counts_alpha_diversity.txt")

# Re-order levels for sample biomass
md_long$biomass_plate <- factor(md_long$biomass_plate, levels = c("high", "low", "not provided"))

# Re-order levels for extraction kits
md_long$extraction_kit_round <- factor(md_long$extraction_kit_round, levels = c("PowerSoil r1", "PowerSoil Pro", "MagMax Beta", "Norgen", "Homebrew", "PowerSoil r2", "MagMAX Microbiome","NucleoMag Food", "Zymo MagBead"))

# Re-name levels for extraction kits
levels(md_long$extraction_kit_round) = c("PowerSoil r1", "PowerSoil Pro", "MagMax Beta", "Norgen", "Homebrew", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")

# Re-order levels for sample type 2
md_long$sample_type_2 <- factor(md_long$sample_type_2, levels = c("keyboard", "floor", "kimchi", "yogurt", "seawater", "freshwater", "bare soil", "rhizosphere soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human saliva", "human milk", "ZymoBIOMICS Microbial Community Standard I", "ZymoBIOMICS Microbial Community Standard II", "PCR extraction control"))

# Re-name levels for sample type 2
levels(md_long$sample_type_2) <- c("Surface, keyboard", "Surface, floor tile", "Food, kimchi", "Food, yogurt", "Water, saline", "Water, non-saline", "Soil, bare", "Soil, rhizosphere", "Mouse tissue, jejunum ", "Mouse feces", "Cat feces", "Human feces", "Human urine, female", "Human  urine, male", "Human skin, foot", "Human skin, armpit", "Human saliva", "Human milk", "mock community I", "mock community II", "PCR extraction control")

# Subset data to remove round 3
md_long_clean_0 <- subset(md_long,
                          md_long$round != 3)

# Subset data to remove round 1 Zymo
md_long_clean_1 <- subset(md_long_clean_0,
                          md_long_clean_0$extraction_kit_round != 'Zymo MagBead' |
                            md_long_clean_0$round != 1)

# Subset further to exclude positive and negative controls and samples with zero-or-negative DNA concentrations
md_long_clean <- subset(md_long_clean_1,
                        md_long_clean_1$sample_type != 'control blank' &
                          md_long_clean_1$sample_type != 'mock community' &
                          md_long_clean_1$dna_conc_ng_ul > 0)

# Subset data to remove round 1 MagMAX Beta
md_long_clean_seven_kits <- subset(md_long_clean,
                                   md_long_clean$extraction_kit_round != 'MagMax Beta')

# Subset further to exclude Homebrew
md_long_clean_six_kits <- subset(md_long_clean_seven_kits,
                                 md_long_clean_seven_kits$extraction_kit_round != 'Homebrew')


# Plot DNA concentrations - Boxplots comparing averages across kits for each sample type
#######################################################################################################

# Boxplot - All kits
## x = extraction_kit_round
## y = dna_conc_ng_ul
## facets = sample_type_2
ggplot(md_long_clean, aes(x = extraction_kit_round, y = as.numeric(dna_conc_ng_ul))) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit_round)) +
  geom_jitter(, width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, scales = "free_y") +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction kit") +
  ylab("DNA yield (ng/µL)") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 5.5, color = "gray60") +
  scale_fill_manual(values = c("blue", "steelblue1", "green", "chocolate4", "gray60", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.background = element_rect(fill="white"))


# Boxplot - All kits
## x = extraction_kit_round
## y = log10(dna_conc_ng_ul)
## facets = sample_type_2
ggplot(md_long_clean, aes(x = extraction_kit_round, y = as.numeric(log10(dna_conc_ng_ul + 1)))) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit_round)) +
  geom_jitter(, width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, scales = "free_y") +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction kit") +
  ylab("log(DNA yield [ng/µL] + 1)") +
  geom_hline(yintercept = log(1.1), linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 5.5, color = "gray60") +
  scale_fill_manual(values = c("blue", "steelblue1", "green", "chocolate4", "gray60", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.background = element_rect(fill="white"))


# Boxplot - Six kits for manuscript
## x = extraction_kit_round
## y = dna_conc_ng_ul
## facets = sample_type_2
ggplot(md_long_clean_six_kits, aes(x = extraction_kit_round, y = as.numeric(dna_conc_ng_ul))) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit_round)) +
  geom_jitter(, width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, 
             scales = "free_y",
             ncol = 6) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction kit") +
  ylab("\n\nDNA yield (ng/µL)") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 3.5, color = "gray60") +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.background = element_rect(fill="white"))


# Boxplot - Six kits for manuscript (Figure SX)
## x = extraction_kit_round
## y = log(dna_conc_ng_ul)
## facets = sample_type_2
ggplot(md_long_clean_six_kits, aes(x = extraction_kit_round, y = as.numeric(log10(dna_conc_ng_ul + 1)))) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit_round)) +
  geom_jitter(, width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, 
             scales = "free_y",
             ncol = 6) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction kit") +
  ylab("\nlog(DNA yield [ng/µL] + 1)") +
  geom_hline(yintercept = log(1.1), linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 3.5, color = "gray60") +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.background = element_rect(fill="white"))


# Plot 16S read counts - Boxplots comparing averages across kits for each sample type
#######################################################################################################

# Boxplot - All kits
## x = extraction_kit_round
## y = 16S_deblur_reads
## facets = sample_type_2
ggplot(md_long_clean, aes(x = extraction_kit_round, y = (`16s_deblur_reads`/1000))) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit_round)) +
  geom_jitter(, width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, 
             scales = "free_y",
             ncol = 6) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction kit") +
  ylab("Quality-filtered 16S reads (x1000)") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 5.5, color = "gray60") +
  scale_fill_manual(values = c("blue", "steelblue1", "green", "chocolate4", "gray60", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.background =element_rect(fill="white"))


# Boxplot - Six kits for manuscript (Figure SX)
## x = extraction_kit_round
## y = 16S_deblur_reads
## facets = sample_type_2
ggplot(md_long_clean_six_kits, aes(x = extraction_kit_round, y = (`16s_deblur_reads`/1000))) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit_round)) +
  geom_jitter(, width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, 
             scales = "free_y",
             ncol = 6) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction kit") +
  ylab("Quality-filtered 16S reads (x1000)") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 3.5, color = "gray60") +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.background =element_rect(fill="white"))


# Plot ITS read-counts - Boxplots comparing averages across kits for each sample type
#######################################################################################################

# Boxplot - Including Homebrew kit
## x = extraction_kit_round
## y = its_deblur_reads
## facets = sample_type_2
ggplot(md_long_clean, aes(x = extraction_kit_round, y = (its_deblur_reads))) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit_round)) +
  geom_jitter(, width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, scales = "free_y") +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction kit") +
  ylab("Quality-filtered ITS reads") +
  geom_vline(xintercept = 5.5, color = "gray60") +
  scale_fill_manual(values = c("blue", "steelblue1", "green", "chocolate4", "gray60", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.background =element_rect(fill="white"))


# Boxplot - Six kits for manuscript (Figure SX)
## x = extraction_kit_round
## y = its_deblur_reads
## facets = sample_type_2
ggplot(md_long_clean_six_kits, aes(x = extraction_kit_round, y = (its_deblur_reads/1000))) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit_round)) +
  geom_jitter(, width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, 
             scales = "free_y",
             ncol = 6) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction kit") +
  ylab("Quality-filtered ITS reads (x1000)") +
  geom_vline(xintercept = 3.5, color = "gray60") +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.background =element_rect(fill="white"))


# Plot Shotgun metagenomic read-counts - Boxplots comparing averages across kits for each sample type
#######################################################################################################

# Boxplot - All kits
## x = extraction_kit_round
## y = shotgun_woltka_reads
## facets = sample_type_2
ggplot(md_long_clean, aes(x = extraction_kit_round, y = (shotgun_woltka_reads/1000000))) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit_round)) +
  geom_jitter(, width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, scales = "free_y") +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 5.5, color = "gray60") +
  scale_fill_manual(values = c("blue", "steelblue1", "green", "chocolate4", "gray60", "blue", "green4", "white", "yellow")) +
  xlab("extraction kit") +
  ylab("Host- and quality-filtered metagenomic reads (x1E6)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.background =element_rect(fill="white"))


# Boxplot - Six kits for manuscript (Figure SX)
## x = extraction_kit_round
## y = shotgun_woltka_reads
## facets = sample_type_2
ggplot(md_long_clean_six_kits, aes(x = extraction_kit_round, y = (shotgun_woltka_reads/1000000))) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit_round)) +
  geom_jitter(, width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, 
             scales = "free_y",
             ncol = 6) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 3.5, color = "gray60") +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  xlab("extraction kit") +
  ylab("Host- and quality-filtered\nmetagenomic reads (x1E6)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.background =element_rect(fill="white"))


# Cowplot - DNA yields and read counts
#######################################################################################################

# DNA yield plot
dna_plot <- ggplot(md_long_clean_six_kits, aes(x = extraction_kit_round, y = as.numeric(log10(dna_conc_ng_ul + 1)))) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit_round)) +
  geom_jitter(, width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, 
             scales = "free_y",
             ncol = 6) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction kit") +
  ylab("\nlog(DNA yield [ng/µL] + 1)") +
  geom_hline(yintercept = log(1.1), linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 3.5, color = "gray60") +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.background = element_rect(fill="white"))


# 16S deblur reads plot
reads_16S_plot <- ggplot(md_long_clean_six_kits, aes(x = extraction_kit_round, y = (`16s_deblur_reads`/1000))) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit_round)) +
  geom_jitter(, width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, 
             scales = "free_y",
             ncol = 6) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction kit") +
  ylab("Quality-filtered 16S reads (x1000)") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 3.5, color = "gray60") +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.background =element_rect(fill="white"))


# ITS deblur reads plot
reads_ITS_plot <- ggplot(md_long_clean_six_kits, aes(x = extraction_kit_round, y = (its_deblur_reads/1000))) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit_round)) +
  geom_jitter(, width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, 
             scales = "free_y",
             ncol = 6) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction kit") +
  ylab("Quality-filtered ITS reads (x1000)") +
  geom_vline(xintercept = 3.5, color = "gray60") +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.background =element_rect(fill="white"))


# Shotgun reads plot
reads_shotgun_plot <- ggplot(md_long_clean_six_kits, aes(x = extraction_kit_round, y = (shotgun_woltka_reads/1000000))) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit_round)) +
  geom_jitter(, width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, 
             scales = "free_y",
             ncol = 6) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 3.5, color = "gray60") +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  xlab("extraction kit") +
  ylab("Host- and quality-filtered\nmetagenomic reads (x1E6)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.background =element_rect(fill="white"))


# Cowplot - DNA yield and 16S reads
figure_dna_yield_16s_reads <- plot_grid(dna_plot,
                                        reads_16S_plot,
                                        labels = c("A", "B"),
                                        label_size = 24,
                                        label_fontfamily = 'sans',
                                        ncol = 1)

figure_dna_yield_16s_reads

save_plot('figure_boxplot_combined_dna_yield_16S_read_count_11x9.5.svg',
          figure_dna_yield_16s_reads,
          base_width = 11,
          base_height = 9.5)


# Cowplot - DNA yield and all read counts
figure_dna_yield_its_reads <- plot_grid(reads_ITS_plot,
                                    reads_shotgun_plot,
                                    labels = c("A", "B"),
                                    label_size = 24,
                                    label_fontfamily = 'sans',
                                    ncol = 1)

figure_dna_yield_its_reads

save_plot('figure_boxplot_combined_ITS_and_shotgun_read_counts_11x9.5.svg',
          figure_dna_yield_its_reads,
          base_width = 11,
          base_height = 9.5)


# Cowplot - DNA yield and 16S reads
figure_dna_yield_16s_reads <- plot_grid(dna_plot,
                                        reads_16S_plot,
                                        labels = c("A", "B"),
                                        label_size = 18,
                                        label_fontfamily = 'sans',
                                        ncol = 1)

figure_dna_yield_16s_reads

save_plot('figure_boxplot_combined_dna_yield_16S_read_count_11x9.5.svg',
          figure_dna_yield_16s_reads,
          base_width = 11,
          base_height = 9.5)


# Cowplot - DNA yield and all read counts
figure_dna_yield_reads <- plot_grid(dna_plot,
                                        reads_16S_plot,
                                        reads_ITS_plot,
                                        reads_shotgun_plot,
                                        labels = c("A", "B", "C", "D"),
                                        label_size = 32,
                                        label_fontfamily = 'sans',
                                        ncol = 1)

figure_dna_yield_reads

save_plot('figure_boxplot_combined_dna_yield_read_counts_17x22.tiff',
          figure_dna_yield_reads,
          base_width = 17,
          base_height = 22)


# Scatterplot - DNA yield and 16S read counts
#########################################################################################
ggplot(md_long_clean_six_kits, aes(x = dna_conc_ng_ul, y = (`16s_deblur_reads`))) +
  geom_point() +
  xlab("DNA yield (ng/µL)") +
  ylab("16S Deblur reads") +
  geom_smooth(method = 'lm') +
  theme_bw() +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        text = element_text(size = 14))


#########################################################################################
# Well-to-well contamination - heatmap visualization
#########################################################################################
library(viridis)

# Read in data
all_w2w <- read.csv("synthetic_plasmid_results_plotting.csv", sep = ",")
all_w2w$extraction_kit <- factor(all_w2w$extraction_kit, levels = c("PowerSoil r1",
                                   "PowerSoil Pro",
                                   "MagMax Beta",
                                   "Norgen",
                                   "Homebrew",
                                   "PowerSoil r2",
                                   "MagMAX Microbiome",
                                   "NucleoMag Food",
                                   "Zymo MagBead"))


# Plot - All kits
ggplot(all_w2w, aes(x = column, y = reorder(row, desc(row)), fill = plasmid_reads_percent)) +
  geom_tile() +
  facet_wrap(~extraction_kit,
             ncol = 5) +
  scale_fill_viridis(name = "Spike-in reads (%)", discrete = F) +
  scale_x_discrete(limits = c(1:12)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5),
        strip.background =element_rect(fill="white"))


# Plot - Six kits
all_w2w_six_kits <- subset(all_w2w, all_w2w$extraction_kit != "Homebrew" &
                             all_w2w$extraction_kit != "MagMax Beta")
ggplot(all_w2w_six_kits, aes(x = column, y = reorder(row, desc(row)), fill = plasmid_reads_percent)) +
  geom_tile() +
  facet_wrap(~extraction_kit,
             ncol = 7) +
  scale_fill_viridis(name = "Spike-in reads (%)", discrete = F) +
  scale_x_discrete(limits = c(1:12)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5),
        strip.background =element_rect(fill="white"))


# Well-to-well contamination - boxplot
#########################################################################################
install.packages("PMCMRplus")
library(tidyverse)
library(PMCMRplus)

# Read in data
all_w2w <- read.csv("synthetic_plasmid_results_plotting.csv", sep = ",")

# Subset data to include only sink wells
all_w2w_sinkOnly <- subset(all_w2w, all_w2w$well_type == "sink")

# Re-order levels for extraction kits
all_w2w_sinkOnly$extraction_kit <- factor(all_w2w_sinkOnly$extraction_kit, levels = c("PowerSoil r1", "PowerSoil Pro", "MagMax Beta", "Norgen",  "Homebrew", "PowerSoil r2", "MagMAX Microbiome", "NucleoMag Food", "Zymo MagBead"))

# Re-name levels for extraction kits
levels(all_w2w_sinkOnly$extraction_kit) = c("PowerSoil r1", "PowerSoil Pro", "MagMax Beta", "Norgen", "Homebrew", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")

# Subset to include only extraction kits for manuscript
all_w2w_sinkOnly_six_kits <- subset(all_w2w_sinkOnly,
                                    all_w2w_sinkOnly$extraction_kit != "MagMax Beta" &
                                      all_w2w_sinkOnly$extraction_kit != "Homebrew")

# boxplot - including all kits
ggplot(all_w2w_sinkOnly, aes(x = extraction_kit, y = plasmid_reads_log10)) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction kit") +
  ylab("log(spike-in read count in sink well)") +
  ylim(c(-0.01,5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "black", "chocolate4", "gray60", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 6),
        strip.background =element_rect(fill="white"))

# boxplot - including only six kits
ggplot(all_w2w_sinkOnly_six_kits, aes(x = extraction_kit, y = plasmid_reads_log10)) +
  #facet_wrap(~round,
  #           ncol = 1,
  #           scales = "free") +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_kit)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction kit") +
  ylab("log(spike-in read count in sink well)") +
  geom_vline(xintercept = 3.5,
        color = "black",
        size = 0.25) +
  ylim(c(-0.01,5.5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14),
        strip.background =element_rect(fill="white"))


# Kruskal-Wallis tests
all_w2w_sinkOnly_six_kits_r1 <- subset(all_w2w_sinkOnly_six_kits,
                                       all_w2w_sinkOnly_six_kits$round == "Round 1")
all_w2w_sinkOnly_six_kits_r2 <- subset(all_w2w_sinkOnly_six_kits,
                                       all_w2w_sinkOnly_six_kits$round == "Round 2")

kruskal.test(all_w2w_sinkOnly_six_kits_r1$plasmid_reads ~ all_w2w_sinkOnly_six_kits_r1$extraction_kit)
kruskal.test(all_w2w_sinkOnly_six_kits_r2$plasmid_reads ~ all_w2w_sinkOnly_six_kits_r2$extraction_kit)


# Post-hoc Dunn's test
kwAllPairsDunnTest(all_w2w_sinkOnly_six_kits_r1$plasmid_reads, 
                          as.factor(all_w2w_sinkOnly_six_kits_r1$extraction_kit), 
                          p.adjust.method = "hochberg")
kwAllPairsDunnTest(all_w2w_sinkOnly_six_kits_r2$plasmid_reads, 
                   as.factor(all_w2w_sinkOnly_six_kits_r2$extraction_kit), 
                   p.adjust.method = "hochberg")


#########################################################################################
# Feature abundance correlations
#########################################################################################
library(tidyverse)
library(ggpubr)
library(cowplot)
library(svglite)

# 16S
#########################################################################################

# Read in long-format files generated with python notebook
powersoil_r2_16S_genera_corr_data <- read_tsv('table_correlation_16S_genera_PowerSoil_r2_swap.txt')
powersoil_pro_16S_genera_corr_data <- read_tsv('table_correlation_16S_genera_PowerSoil_Pro.txt')
magmax_16S_genera_corr_data <- read_tsv('table_correlation_16S_genera_MagMAX_Microbiome.txt')
nucleomag_16S_genera_corr_data <- read_tsv('table_correlation_16S_genera_NucleoMag_Food.txt')
norgen_16S_genera_corr_data <- read_tsv('table_correlation_16S_genera_Norgen.txt')
zymo_16S_genera_corr_data <- read_tsv('table_correlation_16S_genera_Zymo_MagBead_swap.txt')

# Test for normality
ggqqplot(log10(powersoil_r2_16S_genera_corr_data$`PowerSoil r2`))
## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
powersoil_r2_16S_genera_corr_kendall <- cor.test(x = powersoil_r2_16S_genera_corr_data$`PowerSoil r2`,
         y = powersoil_r2_16S_genera_corr_data$`PowerSoil r1`,
         method = "kendall")
powersoil_pro_16S_genera_corr_kendall <- cor.test(x = powersoil_pro_16S_genera_corr_data$`PowerSoil Pro`,
                                                  y = powersoil_pro_16S_genera_corr_data$`PowerSoil r1`,
                                                 method = "kendall")
magmax_16S_genera_corr_kendall <- cor.test(x = magmax_16S_genera_corr_data$`MagMAX Microbiome`,
                                                  y = magmax_16S_genera_corr_data$`PowerSoil r2`,
                                                  method = "kendall")
nucleomag_16S_genera_corr_kendall <- cor.test(x = nucleomag_16S_genera_corr_data$`NucleoMag Food`,
                                           y = nucleomag_16S_genera_corr_data$`PowerSoil r2`,
                                           method = "kendall")
norgen_16S_genera_corr_kendall <- cor.test(x = norgen_16S_genera_corr_data$Norgen,
                                              y = norgen_16S_genera_corr_data$`PowerSoil r1`,
                                              method = "kendall")
zymo_16S_genera_corr_kendall <- cor.test(x = zymo_16S_genera_corr_data$`Zymo MagBead`,
                                           y = zymo_16S_genera_corr_data$`PowerSoil r2`,
                                           method = "kendall")

# View correlations
powersoil_r2_16S_genera_corr_kendall
powersoil_pro_16S_genera_corr_kendall
magmax_16S_genera_corr_kendall
nucleomag_16S_genera_corr_kendall
norgen_16S_genera_corr_kendall
zymo_16S_genera_corr_kendall


# Plot correlations
powersoil_corr_16S <- ggplot(data = powersoil_r2_16S_genera_corr_data,
       aes(x = `PowerSoil r1`,
           y = `PowerSoil r2`)) +
  geom_point(fill = 'blue',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()

powersoil_pro_corr_16S <- ggplot(data = powersoil_pro_16S_genera_corr_data,
                             aes(x = `PowerSoil r1`,
                                 y = `PowerSoil Pro`)) +
  geom_point(fill = 'steelblue1',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("PowerSoil") +
  theme(text = element_text(size = 14))

magmax_corr_16S <- ggplot(data = magmax_16S_genera_corr_data,
                                 aes(x = `PowerSoil r2`,
                                     y = `MagMAX Microbiome`)) +
  geom_point(fill = 'green4',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("PowerSoil") +
  theme(text = element_text(size = 14))

nucleomag_corr_16S <- ggplot(data = nucleomag_16S_genera_corr_data,
                          aes(x = `PowerSoil r2`,
                              y = `NucleoMag Food`)) +
  geom_point(fill = 'white',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("PowerSoil") +
  theme(text = element_text(size = 14))

norgen_corr_16S <- ggplot(data = norgen_16S_genera_corr_data,
                             aes(x = `PowerSoil r1`,
                                 y = `Norgen`)) +
  geom_point(fill = 'chocolate4',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("PowerSoil") +
  theme(text = element_text(size = 14))

zymo_corr_16S <- ggplot(data = zymo_16S_genera_corr_data,
                          aes(x = `PowerSoil r2`,
                              y = `Zymo MagBead`)) +
  geom_point(fill = 'yellow',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("PowerSoil") +
  theme(text = element_text(size = 14))

figure_corr_16S <- plot_grid(powersoil_pro_corr_16S,
          norgen_corr_16S,
          magmax_corr_16S,
          nucleomag_corr_16S,
          zymo_corr_16S,
          powersoil_corr_16S,
          labels = c("A", "B", "C", "D", "E", "F"),
          label_size = 18,
          label_fontfamily = 'sans',
          ncol = 3)

figure_corr_16S

save_plot('figure_scatter_corr_16S_genera_10x6.svg',
          figure_corr_16S,
          base_height = 6,
          base_width = 10)

save_plot('figure_scatter_corr_16S_genera_10x6.png',
          figure_corr_16S,
          base_height = 6,
          base_width = 10)

save_plot('figure_scatter_corr_16S_genera_10x6.tiff',
          figure_corr_16S,
          base_height = 6,
          base_width = 10)


# ITS
#########################################################################################

# Read in long-format files generated with python notebook
powersoil_r2_ITS_genera_corr_data <- read_tsv('table_correlation_ITS_genera_PowerSoil_r2_swap.txt')
powersoil_pro_ITS_genera_corr_data <- read_tsv('table_correlation_ITS_genera_PowerSoil_Pro.txt')
magmax_ITS_genera_corr_data <- read_tsv('table_correlation_ITS_genera_MagMAX_Microbiome.txt')
nucleomag_ITS_genera_corr_data <- read_tsv('table_correlation_ITS_genera_NucleoMag_Food.txt')
norgen_ITS_genera_corr_data <- read_tsv('table_correlation_ITS_genera_Norgen.txt')
zymo_ITS_genera_corr_data <- read_tsv('table_correlation_ITS_genera_Zymo_MagBead_swap.txt')

# Test for normality
ggqqplot(log10(powersoil_r2_ITS_genera_corr_data$`PowerSoil r2`))
## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
powersoil_r2_ITS_genera_corr_kendall <- cor.test(x = powersoil_r2_ITS_genera_corr_data$`PowerSoil r2`,
                                                 y = powersoil_r2_ITS_genera_corr_data$`PowerSoil r1`,
                                                 method = "kendall")
powersoil_pro_ITS_genera_corr_kendall <- cor.test(x = powersoil_pro_ITS_genera_corr_data$`PowerSoil Pro`,
                                                  y = powersoil_pro_ITS_genera_corr_data$`PowerSoil r1`,
                                                  method = "kendall")
magmax_ITS_genera_corr_kendall <- cor.test(x = magmax_ITS_genera_corr_data$`MagMAX Microbiome`,
                                           y = magmax_ITS_genera_corr_data$`PowerSoil r2`,
                                           method = "kendall")
nucleomag_ITS_genera_corr_kendall <- cor.test(x = nucleomag_ITS_genera_corr_data$`NucleoMag Food`,
                                              y = nucleomag_ITS_genera_corr_data$`PowerSoil r2`,
                                              method = "kendall")
norgen_ITS_genera_corr_kendall <- cor.test(x = norgen_ITS_genera_corr_data$Norgen,
                                           y = norgen_ITS_genera_corr_data$`PowerSoil r1`,
                                           method = "kendall")
zymo_ITS_genera_corr_kendall <- cor.test(x = zymo_ITS_genera_corr_data$`Zymo MagBead`,
                                         y = zymo_ITS_genera_corr_data$`PowerSoil r2`,
                                         method = "kendall")

# View correlations
powersoil_r2_ITS_genera_corr_kendall
powersoil_pro_ITS_genera_corr_kendall
magmax_ITS_genera_corr_kendall
nucleomag_ITS_genera_corr_kendall
norgen_ITS_genera_corr_kendall
zymo_ITS_genera_corr_kendall


# Plot correlations
powersoil_corr_ITS <- ggplot(data = powersoil_r2_ITS_genera_corr_data,
                             aes(x = `PowerSoil r1`,
                                 y = `PowerSoil r2`)) +
  geom_point(fill = 'blue',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()

powersoil_pro_corr_ITS <- ggplot(data = powersoil_pro_ITS_genera_corr_data,
                                 aes(x = `PowerSoil r1`,
                                     y = `PowerSoil Pro`)) +
  geom_point(fill = 'steelblue1',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("PowerSoil") +
  theme(text = element_text(size = 14))

magmax_corr_ITS <- ggplot(data = magmax_ITS_genera_corr_data,
                          aes(x = `PowerSoil r2`,
                              y = `MagMAX Microbiome`)) +
  geom_point(fill = 'green4',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("PowerSoil") +
  theme(text = element_text(size = 14))

nucleomag_corr_ITS <- ggplot(data = nucleomag_ITS_genera_corr_data,
                             aes(x = `PowerSoil r2`,
                                 y = `NucleoMag Food`)) +
  geom_point(fill = 'white',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("PowerSoil") +
  theme(text = element_text(size = 14))

norgen_corr_ITS <- ggplot(data = norgen_ITS_genera_corr_data,
                          aes(x = `PowerSoil r1`,
                              y = `Norgen`)) +
  geom_point(fill = 'chocolate4',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("PowerSoil") +
  theme(text = element_text(size = 14))

zymo_corr_ITS <- ggplot(data = zymo_ITS_genera_corr_data,
                        aes(x = `PowerSoil r2`,
                            y = `Zymo MagBead`)) +
  geom_point(fill = 'yellow',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("PowerSoil") +
  theme(text = element_text(size = 14))

figure_corr_ITS <- plot_grid(powersoil_pro_corr_ITS,
                             norgen_corr_ITS,
                             magmax_corr_ITS,
                             nucleomag_corr_ITS,
                             zymo_corr_ITS,
                             powersoil_corr_ITS,
                             labels = c("A", "B", "C", "D", "E", "F"),
                             label_size = 18,
                             label_fontfamily = 'sans',
                             ncol = 3)

figure_corr_ITS

save_plot('figure_scatter_corr_ITS_genera_10x6.svg',
          figure_corr_ITS,
          base_height = 6,
          base_width = 10)

save_plot('figure_scatter_corr_ITS_genera_10x6.png',
          figure_corr_ITS,
          base_height = 6,
          base_width = 10)

save_plot('figure_scatter_corr_ITS_genera_10x6.tiff',
          figure_corr_ITS,
          base_height = 6,
          base_width = 10)


# Shotgun
#########################################################################################

# Read in long-format files generated with python notebook
powersoil_r2_shotgun_genera_corr_data <- read_tsv('table_correlation_shotgun_genera_PowerSoil_r2_swap.txt')
powersoil_pro_shotgun_genera_corr_data <- read_tsv('table_correlation_shotgun_genera_PowerSoil_Pro.txt')
magmax_shotgun_genera_corr_data <- read_tsv('table_correlation_shotgun_genera_MagMAX_Microbiome.txt')
nucleomag_shotgun_genera_corr_data <- read_tsv('table_correlation_shotgun_genera_NucleoMag_Food.txt')
norgen_shotgun_genera_corr_data <- read_tsv('table_correlation_shotgun_genera_Norgen.txt')
zymo_shotgun_genera_corr_data <- read_tsv('table_correlation_shotgun_genera_Zymo_MagBead_swap.txt')

# Test for normality
ggqqplot(log10(powersoil_r2_shotgun_genera_corr_data$`PowerSoil r2`))
## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
powersoil_r2_shotgun_genera_corr_kendall <- cor.test(x = powersoil_r2_shotgun_genera_corr_data$`PowerSoil r2`,
                                                     y = powersoil_r2_shotgun_genera_corr_data$`PowerSoil r1`,
                                                     method = "kendall")
powersoil_pro_shotgun_genera_corr_kendall <- cor.test(x = powersoil_pro_shotgun_genera_corr_data$`PowerSoil Pro`,
                                                      y = powersoil_pro_shotgun_genera_corr_data$`PowerSoil r1`,
                                                      method = "kendall")
magmax_shotgun_genera_corr_kendall <- cor.test(x = magmax_shotgun_genera_corr_data$`MagMAX Microbiome`,
                                               y = magmax_shotgun_genera_corr_data$`PowerSoil r2`,
                                               method = "kendall")
nucleomag_shotgun_genera_corr_kendall <- cor.test(x = nucleomag_shotgun_genera_corr_data$`NucleoMag Food`,
                                                  y = nucleomag_shotgun_genera_corr_data$`PowerSoil r2`,
                                                  method = "kendall")
norgen_shotgun_genera_corr_kendall <- cor.test(x = norgen_shotgun_genera_corr_data$Norgen,
                                               y = norgen_shotgun_genera_corr_data$`PowerSoil r1`,
                                               method = "kendall")
zymo_shotgun_genera_corr_kendall <- cor.test(x = zymo_shotgun_genera_corr_data$`Zymo MagBead`,
                                             y = zymo_shotgun_genera_corr_data$`PowerSoil r2`,
                                             method = "kendall")

# View correlations
powersoil_r2_shotgun_genera_corr_kendall
powersoil_pro_shotgun_genera_corr_kendall
magmax_shotgun_genera_corr_kendall
nucleomag_shotgun_genera_corr_kendall
norgen_shotgun_genera_corr_kendall
zymo_shotgun_genera_corr_kendall


# Plot correlations
powersoil_corr_shotgun <- ggplot(data = powersoil_r2_shotgun_genera_corr_data,
                                 aes(x = `PowerSoil r1`,
                                     y = `PowerSoil r2`)) +
  geom_point(fill = 'blue',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()

powersoil_pro_corr_shotgun <- ggplot(data = powersoil_pro_shotgun_genera_corr_data,
                                     aes(x = `PowerSoil r1`,
                                         y = `PowerSoil Pro`)) +
  geom_point(fill = 'steelblue1',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("PowerSoil") +
  theme(text = element_text(size = 14))

magmax_corr_shotgun <- ggplot(data = magmax_shotgun_genera_corr_data,
                              aes(x = `PowerSoil r2`,
                                  y = `MagMAX Microbiome`)) +
  geom_point(fill = 'green4',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("PowerSoil") +
  theme(text = element_text(size = 14))

nucleomag_corr_shotgun <- ggplot(data = nucleomag_shotgun_genera_corr_data,
                                 aes(x = `PowerSoil r2`,
                                     y = `NucleoMag Food`)) +
  geom_point(fill = 'white',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("PowerSoil") +
  theme(text = element_text(size = 14))

norgen_corr_shotgun <- ggplot(data = norgen_shotgun_genera_corr_data,
                              aes(x = `PowerSoil r1`,
                                  y = `Norgen`)) +
  geom_point(fill = 'chocolate4',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("PowerSoil") +
  theme(text = element_text(size = 14))

zymo_corr_shotgun <- ggplot(data = zymo_shotgun_genera_corr_data,
                            aes(x = `PowerSoil r2`,
                                y = `Zymo MagBead`)) +
  geom_point(fill = 'yellow',
             shape = 21) +
  geom_smooth(method = 'lm',
              color = 'blue',
              size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("PowerSoil") +
  theme(text = element_text(size = 14))

figure_corr_shotgun <- plot_grid(powersoil_pro_corr_shotgun,
                                 norgen_corr_shotgun,
                                 magmax_corr_shotgun,
                                 nucleomag_corr_shotgun,
                                 zymo_corr_shotgun,
                                 powersoil_corr_shotgun,
                                 labels = c("A", "B", "C", "D", "E", "F"),
                                 label_size = 18,
                                 label_fontfamily = 'sans',
                                 ncol = 3)

figure_corr_shotgun

save_plot('figure_scatter_corr_shotgun_genera_10x6.svg',
          figure_corr_shotgun,
          base_height = 6,
          base_width = 10)

save_plot('figure_scatter_corr_shotgun_genera_10x6.png',
          figure_corr_shotgun,
          base_height = 6,
          base_width = 10)

save_plot('figure_scatter_corr_shotgun_genera_10x6.tiff',
          figure_corr_shotgun,
          base_height = 6,
          base_width = 10)


# Cowplot - all data layers
#########################################################################################
figure_corr_abundances <- plot_grid(powersoil_pro_corr_16S,
                                 norgen_corr_16S,
                                 magmax_corr_16S,
                                 nucleomag_corr_16S,
                                 zymo_corr_16S,
                                 powersoil_pro_corr_ITS,
                                 norgen_corr_ITS,
                                 magmax_corr_ITS,
                                 nucleomag_corr_ITS,
                                 zymo_corr_ITS,
                                 powersoil_pro_corr_shotgun,
                                 norgen_corr_shotgun,
                                 magmax_corr_shotgun,
                                 nucleomag_corr_shotgun,
                                 zymo_corr_shotgun,
                                 labels = c("A", "B", "C", "D", "E",
                                            "F", "G", "H", "I", "J",
                                            "K", "L", "M", "N", "O"),
                                 label_size = 24,
                                 label_fontfamily = 'sans',
                                 ncol = 5,
                                 hjust = -0.25)

save_plot('figure_scatter_corr_abundances_20x10.tiff',
          figure_corr_abundances,
          base_width = 20,
          base_height = 10)

save_plot('figure_scatter_corr_abundances_20x10.svg',
          figure_corr_abundances,
          base_width = 20,
          base_height = 10)


#########################################################################################
# Technical replicate distance analysis
#########################################################################################

# Read in data
tech_rep_16S_hbm_jaccard <- read_tsv("data_tech_reps_16S_hbm_jaccard.txt")
tech_rep_16S_hbm_rpca <- read_tsv("data_tech_reps_16S_hbm_rpca.txt")
tech_rep_16S_hbm_unifrac <- read_tsv("data_tech_reps_16S_hbm_unifrac.txt")
tech_rep_16S_hbm_wunifrac <- read_tsv("data_tech_reps_16S_hbm_wunifrac.txt")
tech_rep_16S_lbm_jaccard <- read_tsv("data_tech_reps_16S_lbm_jaccard.txt")
tech_rep_16S_lbm_rpca <- read_tsv("data_tech_reps_16S_lbm_rpca.txt")
tech_rep_16S_lbm_unifrac <- read_tsv("data_tech_reps_16S_lbm_unifrac.txt")
tech_rep_16S_lbm_wunifrac <- read_tsv("data_tech_reps_16S_lbm_wunifrac.txt")
tech_rep_its_hbm_jaccard <- read_tsv("data_tech_reps_its_hbm_jaccard.txt")
tech_rep_its_hbm_rpca <- read_tsv("data_tech_reps_its_hbm_rpca.txt")
tech_rep_its_lbm_jaccard <- read_tsv("data_tech_reps_its_lbm_jaccard.txt")
tech_rep_its_lbm_rpca <- read_tsv("data_tech_reps_its_lbm_rpca.txt")
tech_rep_shotgun_hbm_jaccard <- read_tsv("data_tech_reps_shotgun_hbm_jaccard.txt")
tech_rep_shotgun_hbm_rpca <- read_tsv("data_tech_reps_shotgun_hbm_rpca.txt")
tech_rep_shotgun_hbm_unifrac <- read_tsv("data_tech_reps_shotgun_hbm_unifrac.txt")
tech_rep_shotgun_hbm_wunifrac <- read_tsv("data_tech_reps_shotgun_hbm_wunifrac.txt")
tech_rep_shotgun_lbm_jaccard <- read_tsv("data_tech_reps_shotgun_lbm_jaccard.txt")
tech_rep_shotgun_lbm_rpca <- read_tsv("data_tech_reps_shotgun_lbm_rpca.txt")
tech_rep_shotgun_lbm_unifrac <- read_tsv("data_tech_reps_shotgun_lbm_unifrac.txt")
tech_rep_shotgun_lbm_wunifrac <- read_tsv("data_tech_reps_shotgun_lbm_wunifrac.txt")

# Merge data
tech_rep_16S_jaccard <- rbind(tech_rep_16S_hbm_jaccard, tech_rep_16S_lbm_jaccard)
tech_rep_16S_rpca <- rbind(tech_rep_16S_hbm_rpca, tech_rep_16S_lbm_rpca)
tech_rep_16S_unifrac <- rbind(tech_rep_16S_hbm_unifrac, tech_rep_16S_lbm_unifrac)
tech_rep_16S_wunifrac <- rbind(tech_rep_16S_hbm_wunifrac, tech_rep_16S_lbm_wunifrac)
tech_rep_its_jaccard <- rbind(tech_rep_its_hbm_jaccard, tech_rep_its_lbm_jaccard)
tech_rep_its_rpca <- rbind(tech_rep_its_hbm_rpca, tech_rep_its_lbm_rpca)
tech_rep_shotgun_jaccard <- rbind(tech_rep_shotgun_hbm_jaccard, tech_rep_shotgun_lbm_jaccard)
tech_rep_shotgun_rpca <- rbind(tech_rep_shotgun_hbm_rpca, tech_rep_shotgun_lbm_rpca)
tech_rep_shotgun_unifrac <- rbind(tech_rep_shotgun_hbm_unifrac, tech_rep_shotgun_lbm_unifrac)
tech_rep_shotgun_wunifrac <- rbind(tech_rep_shotgun_hbm_wunifrac, tech_rep_shotgun_lbm_wunifrac)

# Re-order levels for extraction kits
tech_rep_16S_jaccard$sample1_extraction_kit_round <- factor(tech_rep_16S_jaccard$sample1_extraction_kit_round, levels = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX Microbiome","NucleoMag Food", "Zymo MagBead"))
tech_rep_16S_rpca$sample1_extraction_kit_round <- factor(tech_rep_16S_rpca$sample1_extraction_kit_round, levels = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX Microbiome","NucleoMag Food", "Zymo MagBead"))
tech_rep_16S_unifrac$sample1_extraction_kit_round <- factor(tech_rep_16S_unifrac$sample1_extraction_kit_round, levels = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX Microbiome","NucleoMag Food", "Zymo MagBead"))
tech_rep_16S_wunifrac$sample1_extraction_kit_round <- factor(tech_rep_16S_wunifrac$sample1_extraction_kit_round, levels = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX Microbiome","NucleoMag Food", "Zymo MagBead"))
tech_rep_its_jaccard$sample1_extraction_kit_round <- factor(tech_rep_its_jaccard$sample1_extraction_kit_round, levels = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX Microbiome","NucleoMag Food", "Zymo MagBead"))
tech_rep_its_rpca$sample1_extraction_kit_round <- factor(tech_rep_its_rpca$sample1_extraction_kit_round, levels = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX Microbiome","NucleoMag Food", "Zymo MagBead"))
tech_rep_shotgun_jaccard$sample1_extraction_kit_round <- factor(tech_rep_shotgun_jaccard$sample1_extraction_kit_round, levels = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX Microbiome","NucleoMag Food", "Zymo MagBead"))
tech_rep_shotgun_rpca$sample1_extraction_kit_round <- factor(tech_rep_shotgun_rpca$sample1_extraction_kit_round, levels = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX Microbiome","NucleoMag Food", "Zymo MagBead"))
tech_rep_shotgun_unifrac$sample1_extraction_kit_round <- factor(tech_rep_shotgun_unifrac$sample1_extraction_kit_round, levels = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX Microbiome","NucleoMag Food", "Zymo MagBead"))
tech_rep_shotgun_wunifrac$sample1_extraction_kit_round <- factor(tech_rep_shotgun_wunifrac$sample1_extraction_kit_round, levels = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX Microbiome","NucleoMag Food", "Zymo MagBead"))


# Re-name levels for extraction kits
levels(tech_rep_16S_jaccard$sample1_extraction_kit_round) = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")
levels(tech_rep_16S_rpca$sample1_extraction_kit_round) = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")
levels(tech_rep_16S_unifrac$sample1_extraction_kit_round) = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")
levels(tech_rep_16S_wunifrac$sample1_extraction_kit_round) = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")
levels(tech_rep_its_jaccard$sample1_extraction_kit_round) = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")
levels(tech_rep_its_rpca$sample1_extraction_kit_round) = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")
levels(tech_rep_shotgun_jaccard$sample1_extraction_kit_round) = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")
levels(tech_rep_shotgun_rpca$sample1_extraction_kit_round) = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")
levels(tech_rep_shotgun_unifrac$sample1_extraction_kit_round) = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")
levels(tech_rep_shotgun_wunifrac$sample1_extraction_kit_round) = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")


# Re-order levels for sample type 2
tech_rep_16S_jaccard$sample_type_2 <- factor(tech_rep_16S_jaccard$sample_type_2, levels = c("keyboard", "floor", "kimchi", "yogurt", "seawater", "freshwater", "bare soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human saliva", "human milk"))
tech_rep_16S_rpca$sample_type_2 <- factor(tech_rep_16S_rpca$sample_type_2, levels = c("keyboard", "floor", "kimchi", "yogurt", "seawater", "freshwater", "bare soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human saliva", "human milk"))
tech_rep_16S_unifrac$sample_type_2 <- factor(tech_rep_16S_unifrac$sample_type_2, levels = c("keyboard", "floor", "kimchi", "yogurt", "seawater", "freshwater", "bare soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human saliva", "human milk"))
tech_rep_16S_wunifrac$sample_type_2 <- factor(tech_rep_16S_wunifrac$sample_type_2, levels = c("keyboard", "floor", "kimchi", "yogurt", "seawater", "freshwater", "bare soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human saliva", "human milk"))
tech_rep_its_jaccard$sample_type_2 <- factor(tech_rep_its_jaccard$sample_type_2, levels = c("keyboard", "floor", "kimchi", "yogurt", "seawater", "freshwater", "bare soil", "rhizosphere soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human saliva", "human milk"))
tech_rep_its_rpca$sample_type_2 <- factor(tech_rep_its_rpca$sample_type_2, levels = c("keyboard", "floor", "kimchi", "yogurt", "seawater", "freshwater", "bare soil", "rhizosphere soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human saliva", "human milk"))
tech_rep_shotgun_jaccard$sample_type_2 <- factor(tech_rep_shotgun_jaccard$sample_type_2, levels = c("keyboard", "floor", "kimchi", "yogurt", "seawater", "freshwater", "bare soil", "rhizosphere soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human saliva", "human milk"))
tech_rep_shotgun_rpca$sample_type_2 <- factor(tech_rep_shotgun_rpca$sample_type_2, levels = c("keyboard", "floor", "kimchi", "yogurt", "seawater", "freshwater", "bare soil", "rhizosphere soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human saliva", "human milk"))
tech_rep_shotgun_unifrac$sample_type_2 <- factor(tech_rep_shotgun_unifrac$sample_type_2, levels = c("keyboard", "floor", "kimchi", "yogurt", "seawater", "freshwater", "bare soil", "rhizosphere soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human saliva", "human milk"))
tech_rep_shotgun_wunifrac$sample_type_2 <- factor(tech_rep_shotgun_wunifrac$sample_type_2, levels = c("keyboard", "floor", "kimchi", "yogurt", "seawater", "freshwater", "bare soil", "rhizosphere soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human saliva", "human milk"))


# Re-name levels for sample type 2
levels(tech_rep_16S_jaccard$sample_type_2) <- c("Surface, keyboard", "Surface, floor tile", "Food, kimchi", "Food, yogurt", "Water, saline", "Water, non-saline", "Soil, bare", "Mouse tissue, jejunum ", "Mouse feces", "Cat feces", "Human feces", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human saliva", "Human milk")
levels(tech_rep_16S_rpca$sample_type_2) <- c("Surface, keyboard", "Surface, floor tile", "Food, kimchi", "Food, yogurt", "Water, saline", "Water, non-saline", "Soil, bare", "Mouse tissue, jejunum ", "Mouse feces", "Cat feces", "Human feces", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human saliva", "Human milk")
levels(tech_rep_16S_unifrac$sample_type_2) <- c("Surface, keyboard", "Surface, floor tile", "Food, kimchi", "Food, yogurt", "Water, saline", "Water, non-saline", "Soil, bare", "Mouse tissue, jejunum ", "Mouse feces", "Cat feces", "Human feces", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human saliva", "Human milk")
levels(tech_rep_16S_wunifrac$sample_type_2) <- c("Surface, keyboard", "Surface, floor tile", "Food, kimchi", "Food, yogurt", "Water, saline", "Water, non-saline", "Soil, bare", "Mouse tissue, jejunum ", "Mouse feces", "Cat feces", "Human feces", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human saliva", "Human milk")
levels(tech_rep_its_jaccard$sample_type_2) <- c("Surface, keyboard", "Surface, floor tile", "Food, kimchi", "Food, yogurt", "Water, saline", "Water, non-saline", "Soil, bare", "Soil, rhizosphere", "Mouse tissue, jejunum ", "Mouse feces", "Cat feces", "Human feces", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human saliva", "Human milk")
levels(tech_rep_its_rpca$sample_type_2) <- c("Surface, keyboard", "Surface, floor tile", "Food, kimchi", "Food, yogurt", "Water, saline", "Water, non-saline", "Soil, bare", "Soil, rhizosphere", "Mouse tissue, jejunum ", "Mouse feces", "Cat feces", "Human feces", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human saliva", "Human milk")
levels(tech_rep_shotgun_jaccard$sample_type_2) <- c("Surface, keyboard", "Surface, floor tile", "Food, kimchi", "Food, yogurt", "Water, saline", "Water, non-saline", "Soil, bare", "Soil, rhizosphere", "Mouse tissue, jejunum ", "Mouse feces", "Cat feces", "Human feces", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human saliva", "Human milk")
levels(tech_rep_shotgun_rpca$sample_type_2) <- c("Surface, keyboard", "Surface, floor tile", "Food, kimchi", "Food, yogurt", "Water, saline", "Water, non-saline", "Soil, bare", "Soil, rhizosphere", "Mouse tissue, jejunum ", "Mouse feces", "Cat feces", "Human feces", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human saliva", "Human milk")
levels(tech_rep_shotgun_unifrac$sample_type_2) <- c("Surface, keyboard", "Surface, floor tile", "Food, kimchi", "Food, yogurt", "Water, saline", "Water, non-saline", "Soil, bare", "Soil, rhizosphere", "Mouse tissue, jejunum ", "Mouse feces", "Cat feces", "Human feces", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human saliva", "Human milk")
levels(tech_rep_shotgun_wunifrac$sample_type_2) <- c("Surface, keyboard", "Surface, floor tile", "Food, kimchi", "Food, yogurt", "Water, saline", "Water, non-saline", "Soil, bare", "Soil, rhizosphere", "Mouse tissue, jejunum ", "Mouse feces", "Cat feces", "Human feces", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human saliva", "Human milk")


# Plot data
boxplot_16S_jaccard <- ggplot(tech_rep_16S_jaccard, aes(x = sample1_extraction_kit_round, y = value)) +
  facet_wrap(~sample_type_2) +
  geom_boxplot(outlier.alpha = 0, aes(fill = sample1_extraction_kit_round)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  ylab("Within sample\nJaccard distance") +
  geom_vline(xintercept = 3.5,
             color = "black",
             size = 0.25) +
  ylim(c(-0.01,5.5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.background =element_rect(fill="white"))

boxplot_16S_rpca <- ggplot(tech_rep_16S_rpca, aes(x = sample1_extraction_kit_round, y = value)) +
  facet_wrap(~sample_type_2,
             ncol = 5) +
  geom_boxplot(outlier.alpha = 0, aes(fill = sample1_extraction_kit_round)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  ylab("Within sample\nRobust Aitchison distance") +
  geom_vline(xintercept = 3.5,
             color = "black",
             size = 0.25) +
  ylim(c(-0.01,5.5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.background =element_rect(fill="white"))

boxplot_16S_unifrac <- ggplot(tech_rep_16S_unifrac, aes(x = sample1_extraction_kit_round, y = value)) +
  facet_wrap(~sample_type_2) +
  geom_boxplot(outlier.alpha = 0, aes(fill = sample1_extraction_kit_round)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  ylab("Within sample\nunweighted UniFrac distance") +
  geom_vline(xintercept = 3.5,
             color = "black",
             size = 0.25) +
  ylim(c(-0.01,5.5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.background =element_rect(fill="white"))

boxplot_16S_wunifrac <- ggplot(tech_rep_16S_wunifrac, aes(x = sample1_extraction_kit_round, y = value)) +
  facet_wrap(~sample_type_2) +
  geom_boxplot(outlier.alpha = 0, aes(fill = sample1_extraction_kit_round)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  ylab("Within sample\nweighted UniFrac distance") +
  geom_vline(xintercept = 3.5,
             color = "black",
             size = 0.25) +
  ylim(c(-0.01,5.5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.background =element_rect(fill="white"))

boxplot_its_jaccard <- ggplot(tech_rep_its_jaccard, aes(x = sample1_extraction_kit_round, y = value)) +
  facet_wrap(~sample_type_2) +
  geom_boxplot(outlier.alpha = 0, aes(fill = sample1_extraction_kit_round)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  ylab("Within sample\nJaccard distance") +
  geom_vline(xintercept = 3.5,
             color = "black",
             size = 0.25) +
  ylim(c(-0.01,5.5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.background =element_rect(fill="white"))

boxplot_its_rpca <- ggplot(tech_rep_its_rpca, aes(x = sample1_extraction_kit_round, y = value)) +
  facet_wrap(~sample_type_2) +
  geom_boxplot(outlier.alpha = 0, aes(fill = sample1_extraction_kit_round)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  ylab("Within sample\nRobust Aitchison distance") +
  geom_vline(xintercept = 3.5,
             color = "black",
             size = 0.25) +
  ylim(c(-0.01,5.5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.background =element_rect(fill="white"))

boxplot_shotgun_jaccard <- ggplot(tech_rep_shotgun_jaccard, aes(x = sample1_extraction_kit_round, y = value)) +
  facet_wrap(~sample_type_2) +
  geom_boxplot(outlier.alpha = 0, aes(fill = sample1_extraction_kit_round)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  ylab("Within sample\nJaccard distance") +
  geom_vline(xintercept = 3.5,
             color = "black",
             size = 0.25) +
  ylim(c(-0.01,5.5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.background =element_rect(fill="white"))

boxplot_shotgun_rpca <- ggplot(tech_rep_shotgun_rpca, aes(x = sample1_extraction_kit_round, y = value)) +
  facet_wrap(~sample_type_2) +
  geom_boxplot(outlier.alpha = 0, aes(fill = sample1_extraction_kit_round)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  ylab("Within sample\nRobust Aitchison distance") +
  geom_vline(xintercept = 3.5,
             color = "black",
             size = 0.25) +
  ylim(c(-0.01,5.5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.background =element_rect(fill="white"))

boxplot_shotgun_unifrac <- ggplot(tech_rep_shotgun_unifrac, aes(x = sample1_extraction_kit_round, y = value)) +
  facet_wrap(~sample_type_2) +
  geom_boxplot(outlier.alpha = 0, aes(fill = sample1_extraction_kit_round)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  ylab("Within sample\nunweighted UniFrac distance") +
  geom_vline(xintercept = 3.5,
             color = "black",
             size = 0.25) +
  ylim(c(-0.01,5.5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.background =element_rect(fill="white"))

boxplot_shotgun_wunifrac <- ggplot(tech_rep_shotgun_wunifrac, aes(x = sample1_extraction_kit_round, y = value)) +
  facet_wrap(~sample_type_2) +
  geom_boxplot(outlier.alpha = 0, aes(fill = sample1_extraction_kit_round)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  ylab("Within sample\nweighted UniFrac distance") +
  geom_vline(xintercept = 3.5,
             color = "black",
             size = 0.25) +
  ylim(c(-0.01,5.5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.background =element_rect(fill="white"))


# Cowplots
#########################################################################################
library(tidyverse)
library(ggpubr)
library(cowplot)
library(svglite)
library(scales)

figure_tech_rep_16S <- plot_grid(boxplot_16S_jaccard,
                                 boxplot_16S_rpca,
                                 boxplot_16S_unifrac,
                                 boxplot_16S_wunifrac,
                                    labels = c("A", "B", "C", "D"),
                                    label_size = 24,
                                    label_fontfamily = 'sans',
                                    ncol = 2,
                                    hjust = -0.25)

figure_tech_rep_16S

save_plot('figure_tech_rep_16S_20x15.tiff',
          figure_tech_rep_16S,
          base_width = 20,
          base_height = 15)

save_plot('figure_tech_rep_16S_20x15.svg',
          figure_tech_rep_16S,
          base_width = 20,
          base_height = 15)


figure_tech_rep_its <- plot_grid(boxplot_its_jaccard,
                                 boxplot_its_rpca,
                                 labels = c("A", "B"),
                                 label_size = 24,
                                 label_fontfamily = 'sans',
                                 ncol = 2,
                                 hjust = -0.25)

figure_tech_rep_its

save_plot('figure_tech_rep_its_20x7.5.tiff',
          figure_tech_rep_its,
          base_width = 20,
          base_height = 7.5)
save_plot('figure_tech_rep_its_20x7.5.svg',
          figure_tech_rep_its,
          base_width = 20,
          base_height = 7.5)


figure_tech_rep_shotgun <- plot_grid(boxplot_shotgun_jaccard,
                                 boxplot_shotgun_rpca,
                                 boxplot_shotgun_unifrac,
                                 boxplot_shotgun_wunifrac,
                                 labels = c("A", "B", "C", "D"),
                                 label_size = 24,
                                 label_fontfamily = 'sans',
                                 ncol = 2,
                                 hjust = -0.25)

figure_tech_rep_shotgun

save_plot('figure_tech_rep_shotgun_20x15.tiff',
          figure_tech_rep_shotgun,
          base_width = 20,
          base_height = 15)

save_plot('figure_tech_rep_shotgun_20x15.svg',
          figure_tech_rep_shotgun,
          base_width = 20,
          base_height = 15)

