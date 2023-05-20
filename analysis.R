library(dplyr)
library(IntegrationFeatureHeatmap)
library(ggplot2)
library(factoextra)
library(FactoMineR)
library(reshape2)
library(vegan)
options(scipen = 999)

chromosome_lengths <- readr::read_csv('hg38_chromo_lengths.csv')

sites1 <- readr::read_tsv('/home/ubuntu/projects/Weinberger_lat/intSites.tsv')
sites1$posid <- paste0(sites1$chromosome, sites1$strand, sites1$position)
sites1 <- subset(sites1, reads > 1)

sites2 <- readr::read_tsv('ADA.tsv')
sites2$subject <- paste('ADA', sites2$subject)
sites2$externalSampleID <- as.character(sites2$externalSampleID)
sites3 <- readr::read_tsv('SCID_GTSP5545.tsv')
sites3$subject <- paste('SCID', sites3$subject)
sites3$externalSampleID <- as.character(sites3$externalSampleID)
sitesX <- bind_rows(sites2, sites3)

# Prepare experimental sites for heatmap generation.
s1 <- select(sites1, chromosome, position, subject, internalSampleID)
names(s1) <- c('seqname', 'start', 'subject', 'sample')
s1$end <- s1$start
s1$width <- 1
s1$strand <- '*'
s1$mid <- s1$start
s1$type <- 'insertion'
s1$heatmap_group <- s1$subject
s1 <- select(s1, -subject, -sample)


# Prepare reference sites for heatmap generation.
s2 <- select(sitesX, chromosome, position, subject, internalSampleID)
names(s2) <- c('seqname', 'start', 'subject', 'sample')
s2$end <- s2$start
s2$width <- 1
s2$strand <- '*'
s2$mid <- s2$start
s2$type <- 'insertion'
s2$heatmap_group <- s2$subject
s2 <- select(s2, -subject, -sample)

cleanUpMapNames <- function(x){
  names(x) <- sub('\\.rds', '', names(x))
  names(x) <- sub('\\.RData', '', names(x))
  names(x) <- sub('hg38\\.', '', names(x))
  names(x) <- sub('\\.1000$', ' 1K', names(x))
  names(x) <- sub('\\.10000$', ' 10K', names(x))
  names(x) <- sub('\\.15000$', ' 15K', names(x))
  names(x) <- sub('\\.20000$', ' 20K', names(x))
  names(x) <- sub('\\.25000$', ' 25K', names(x))
  x
}

s <- bind_rows(s1, s2)

# generate random-matched dataframes
random_match_df <- IntegrationFeatureHeatmap::aavenger_sites_random_match(
  aavenger_sites = s,
  chromosome_lengths = chromosome_lengths,
  random_seed_value = 10,
  match_row_number_modifier = 3
)

combined_df <- rbind(s, random_match_df)


# Genomic heatmap.
#-------------------------------------------------------------------------------

# test each integration site for overlap in each feature at each given overlap 
trackPath <- '/home/ubuntu/data/heatMapTracks/genomicMarkers'

combined_overlap_test_results_genomic_ranges <- cleanUpMapNames(IntegrationFeatureHeatmap::test_for_overlaps(
  matched_aavenger_sites_df = combined_df,
  list_of_feature_files = list.files(trackPath, pattern = "\\.(rds|RData)$", full.names = TRUE),
  overlap_ranges_to_test = c(1000, 10000)
))

df_roc <- format_hot_roc_result(hotroc_compare_insertion_to_match(combined_overlap_test_results_genomic_ranges))
df_roc_pvals <- df_roc %>% dplyr::select(feature_window, heatmap_group, p_value, sig_p_value)
df_roc <- df_roc %>% dplyr::select(feature_window, heatmap_group,  ROC_value)

df_roc$heatmap_group <- gsub('\\.', '_', df_roc$heatmap_group)
df_roc_pvals$heatmap_group <- gsub('\\.', '_', df_roc_pvals$heatmap_group)

df_roc$heatmap_group <- factor(df_roc$heatmap_group, 
                               levels = c('C1_Top', 'C2_Top', 'C3_Top', 'C4_Top', 'LF1_Top', 'LF2_Top', 'LF3_Top', 'LF4_Top', 
                                          'C1_Bottom', 'C2_Bottom', 'C3_Bottom', 'C4_Bottom', 'LF1_Bottom', 'LF2_Bottom', 
                                          'LF3_Bottom', 'LF4_Bottom', 'ADA_p402', 'ADA_p404', 'ADA_p410', 'SCID_pP10'))
heatmap <- 
  ggplot() + 
  geom_tile(data = df_roc, aes(x = heatmap_group, y = feature_window, fill = ROC_value), color='black', linewidth = 0.75) + 
  geom_text(data = df_roc_pvals, aes(x = heatmap_group, y = feature_window, label = sig_p_value)) + 
  theme_classic() + 
  theme(axis.title = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12)) + 
  scale_fill_gradientn(colours = c("blue", "grey90", "red"), 
                       na.value = "transparent", 
                       breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = c(0, 0.25, 0.5, 0.75, 1), 
                       limits = c(0, 1)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  scale_x_discrete(expand = c(0, 0)) + 
  labs(fill = "ROC") 
  
ggsave(heatmap, file = 'heatMap1.png', units = 'in', height = 8, width = 9)



# Epigenomic heatmap
#-------------------------------------------------------------------------------
trackPath <- '/home/ubuntu/data/heatMapTracks/epiGenetricMarkers'

combined_overlap_test_results_genomic_ranges <- cleanUpMapNames(IntegrationFeatureHeatmap::test_for_overlaps(
  matched_aavenger_sites_df = combined_df,
  list_of_feature_files = list.files(trackPath, pattern = "\\.(rds|RData)$", full.names = TRUE),
  overlap_ranges_to_test = c(10000)
))

df_roc <- format_hot_roc_result(hotroc_compare_insertion_to_match(combined_overlap_test_results_genomic_ranges))
df_roc_pvals <- df_roc %>% dplyr::select(feature_window, heatmap_group, p_value, sig_p_value)
df_roc <- df_roc %>% dplyr::select(feature_window, heatmap_group,  ROC_value)

df_roc$heatmap_group <- gsub('\\.', '_', df_roc$heatmap_group)
df_roc_pvals$heatmap_group <- gsub('\\.', '_', df_roc_pvals$heatmap_group)

df_roc$heatmap_group <- factor(df_roc$heatmap_group, 
                               levels = c('C1_Top', 'C2_Top', 'C3_Top', 'C4_Top', 'LF1_Top', 'LF2_Top', 'LF3_Top', 'LF4_Top', 
                                          'C1_Bottom', 'C2_Bottom', 'C3_Bottom', 'C4_Bottom', 'LF1_Bottom', 'LF2_Bottom', 
                                          'LF3_Bottom', 'LF4_Bottom', 'ADA_p402', 'ADA_p404', 'ADA_p410', 'SCID_pP10'))
heatmap <- 
  ggplot() + 
  geom_tile(data = df_roc, aes(x = heatmap_group, y = feature_window, fill = ROC_value), color='black', linewidth = 0.75) + 
  geom_text(data = df_roc_pvals, aes(x = heatmap_group, y = feature_window, label = sig_p_value)) + 
  theme_classic() + 
  theme(axis.title = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12)) + 
  scale_fill_gradientn(colours = c("blue", "grey90", "red"), 
                       na.value = "transparent", 
                       breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = c(0, 0.25, 0.5, 0.75, 1), 
                       limits = c(0, 1)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  scale_x_discrete(expand = c(0, 0)) + 
  labs(fill = "ROC") 

ggsave(heatmap, file = 'heatMap2.png', units = 'in', height = 12, width = 8.5)


# Create ROCs for PCA
#-------------------------------------------------------------------------------
trackPath <- '/home/ubuntu/data/heatMapTracks/allMarkers'

combined_overlap_test_results_genomic_ranges <- cleanUpMapNames(IntegrationFeatureHeatmap::test_for_overlaps(
  matched_aavenger_sites_df = combined_df,
  list_of_feature_files = list.files(trackPath, pattern = "\\.(rds|RData)$", full.names = TRUE),
  overlap_ranges_to_test = c(1000, 10000)
))

df_roc <- format_hot_roc_result(hotroc_compare_insertion_to_match(combined_overlap_test_results_genomic_ranges))
df_roc_pvals <- df_roc %>% dplyr::select(feature_window, heatmap_group, p_value, sig_p_value)
df_roc <- df_roc %>% dplyr::select(feature_window, heatmap_group,  ROC_value)

df_roc$heatmap_group <- gsub('\\.', '_', df_roc$heatmap_group)
df_roc_pvals$heatmap_group <- gsub('\\.', '_', df_roc_pvals$heatmap_group)

df_roc$heatmap_group <- factor(df_roc$heatmap_group, 
                               levels = c('C1_Top', 'C2_Top', 'C3_Top', 'C4_Top', 'LF1_Top', 'LF2_Top', 'LF3_Top', 'LF4_Top', 
                                          'C1_Bottom', 'C2_Bottom', 'C3_Bottom', 'C4_Bottom', 'LF1_Bottom', 'LF2_Bottom', 
                                          'LF3_Bottom', 'LF4_Bottom', 'ADA_p402', 'ADA_p404', 'ADA_p410', 'SCID_pP10'))


# Create HIV/ADA/SCID PCA plot
#--------------------------------------------------------------------------------
o <- dcast(df_roc, heatmap_group~feature_window, value.var = 'ROC_value')

o$g1 <- ifelse(grepl('^C', o$heatmap_group), 'Control', 'LatFree')
o$g1 <- ifelse(grepl('ADA|SCID', o$heatmap_group), 'gamma', o$g1)

o$g2 <- ifelse(grepl('Top', o$heatmap_group), 'Top', 'Bottom')
o$g2 <- ifelse(grepl('ADA|SCID', o$heatmap_group), 'gamma', o$g2)

o <- relocate(o, g1, g2, .after = 'heatmap_group')
o$g1 <- factor(o$g1)
o$g2 <- factor(o$g2)

pca.data <- PCA(o[,c(-1,-2,-3)], scale.unit = TRUE, graph = FALSE)

fviz_pca_ind(pca.data,
             geom.ind =  "point",
             col.ind = o$g1, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses.
             legend.title = "Groups",
             mean.point = FALSE,
             pointsize = 3) + labs(title ="")


# Plot without gammas
#-------------------------------------------------------------------------------

o <- dcast(df_roc[!grepl('ADA|SCID', df_roc$heatmap_group),], heatmap_group~feature_window, value.var = 'ROC_value')

o$g1 <- ifelse(grepl('^C', o$heatmap_group), 'Control', 'LatFree')
o$g2 <- ifelse(grepl('Top', o$heatmap_group), 'Top', 'Bottom')

o <- relocate(o, g1, g2, .after = 'heatmap_group')
o$g1 <- factor(o$g1)
o$g2 <- factor(o$g2)

pca.data <- PCA(o[,c(-1,-2,-3)], scale.unit = TRUE, graph = FALSE)

fviz_pca_ind(pca.data,
             geom.ind =  "point",
             col.ind = o$g1, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses.
             legend.title = "Groups",
             mean.point = FALSE,
             pointsize = 3) + labs(title ="")


# permutational MANOVA
#-------------------------------------------------------------------------------
set.seed(1)
adonis2(select(o, -heatmap_group, -g1, -g2) ~ o$g1, method="bray",perm=999)
adonis2(select(o, -heatmap_group, -g1, -g2) ~ o$g2, method="bray",perm=999)

