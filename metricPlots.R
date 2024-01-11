####script for Salvador et al. Analysis of nuclear dysmorphia
####created by Madeline Melzer 20230731, updated by Madeline E Melzer 20231221
####adapted from 2023Feb24_NuclearMetricPlots

library(ggplot2)
library(tidyverse)
library(reshape2)
library(svglite)
library(ggrepel)
library(umap)
library(Seurat)
library(ggfortify)
library(FactoMineR)
library(factoextra)

set.seed(23)

## defining the directories
dataDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/output_metrics/"
plotDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/plots/20231221_finalPlots/"

## loading the data
data = as_tibble(read.table(paste0(dataDirectory, "20231221_allAges.csv"), stringsAsFactors=F, header = T, sep = ","))

## preprocessing
control = data %>% select(-centroid.0,-centroid.1)%>% mutate(concavity_ct.log = log(concavity_ct+1), area_norm = 3.14*(major_axis_length*minor_axis_length)/4, dArea = area_norm-area)

## area cutoffs
lower_bound =  5509.125
upper_bound =  17812.125
controlCombined = control %>% filter(area >= lower_bound & area <= upper_bound & location != "") #added lower and upper bounds from 1.5 IQR area, calcualted in metrics jupyter notebook

############################################################################################
############################## Dimensional Plotting ########################################
############################################################################################

controlCombined <- controlCombined %>% mutate(ID=row_number())

controlCombined.metrics = controlCombined %>% select(eccentricity,
                                                     perimeter,
                                                     dArea,
                                                     area,
                                                     major_axis_length,
                                                     minor_axis_length,
                                                     r_avg,
                                                     rms_avg,
                                                     rms_subsampled,
                                                     x_intercepts,
                                                     concavity_ct,
                                                     circularity,
                                                     solidity,
                                                     equivalent_diameter_area,
                                                     area_bbox,
                                                     area_convex,
                                                     area_filled,
                                                     extent,
                                                     feret_diameter_max,
                                                     compactness,
                                                     roundness,
                                                     convexity,
                                                     aspect_ratio,
                                                     concavity_ct.log,
                                                     perimeter_crofton) %>% mutate(ID = row_number())

controlCombined.meta = controlCombined %>% select(file,
                                                     label,
                                                     age,
                                                     sex,
                                                     animal,
                                                     location,
                                                     eccentricity,
                                                     perimeter,
                                                     dArea,
                                                     area,
                                                     major_axis_length,
                                                     minor_axis_length,
                                                     r_avg,
                                                     rms_avg,
                                                     rms_subsampled,
                                                     x_intercepts,
                                                     concavity_ct,
                                                     circularity,
                                                     solidity,
                                                     equivalent_diameter_area,
                                                     area_bbox,
                                                     area_convex,
                                                     area_filled,
                                                     extent,
                                                     feret_diameter_max,
                                                     compactness,
                                                     roundness,
                                                     convexity,
                                                     aspect_ratio,
                                                     concavity_ct.log,
                                                     perimeter_crofton) %>% mutate(ID = row_number())

# getting rid of id column so the UMAP does not include this identifier
controlCombined.metrics.noID = controlCombined.metrics %>% select(-ID)
controlCombined.metrics.noID = na.omit(controlCombined.metrics.noID)

# PCA and saving PCA
controlCombined.metrics.pca = prcomp(controlCombined.metrics.noID, scale = TRUE)
controlCombined.metrics.pca$rotation = -1*controlCombined.metrics.pca$rotation
#controlCombined.pca.meta = as.data.frame(controlCombined.metrics.pca$x) %>% inner_join(as.data.frame(controlCombined.meta), by = character(0))
#write.csv(controlCombined.pca.meta, "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/output_metrics/UMAPs/20230918_age01_pca.csv")

# elbow plot for pca above, if desired, to help pick number of PCs
explained_variance <- (controlCombined.metrics.pca$sdev)^2 / sum((controlCombined.metrics.pca$sdev)^2)
pca_data <- data.frame(PC = 1:length(explained_variance), Explained_Variance = explained_variance)

# Create the elbow plot
ggplot(pca_data, aes(x = PC, y = Explained_Variance)) +
  geom_line() +
  geom_point() +
  theme_classic(base_size = 18) +
  labs(x = "Principal Component", y = "Explained Variance")
#plot
#ggsave(plot, file = paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/plots/", '20230918_explainedVariance.png'), width = 6, height = 4)

# UMAP based on PCA
controlCombined.pca.umap = umap(controlCombined.metrics.pca$x[,1:11])
controlCombined.umap.df = controlCombined.pca.umap$layout %>%
  as.data.frame() %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  #rename(UMAP1="X1", UMAP2="X2") %>%
  mutate(ID = row_number()) %>%
  inner_join(controlCombined.meta, by = "ID")
write.csv(controlCombined.umap.df, "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/output_metrics/UMAPs/20231221_umap_allAges.csv")


### dot plotting UMAP
controlCombined.umap.df.plot = controlCombined.umap.df %>% filter(age == 24)
plot = ggplot(controlCombined.umap.df.plot, aes(x = UMAP1, y = UMAP2, color = as.factor(location))) +
  geom_point(size = 3) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Arch" = "#6b6ca3", "Desc" = "#c93f55"))
plot
ggsave(plot, file = paste0(plotDirectory, '52wk_dots.svg'), width = 6, height = 4)
ggsave(plot, file = paste0(plotDirectory, '52wk_dots.png'), width = 6, height = 4)

### density plot UMAP

#controlCombined.umap.df.plot = controlCombined.umap.df %>% filter(age == 52)
plot = ggplot() +
  geom_point(data = controlCombined.umap.df.plot, aes(x = UMAP1, y = UMAP2), size = 1, color = "grey93") +
  geom_density_2d(data = controlCombined.umap.df.plot, aes(x = UMAP1, y = UMAP2, color = location), contour_var = "ndensity") +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Arch" = "#6b6ca3", "Desc" = "#c93f55"))
plot
ggsave(plot, file = paste0(plotDirectory, '52wk_density.svg'), width = 6, height = 4)
ggsave(plot, file = paste0(plotDirectory, '52wk_density.png'), width = 6, height = 4)


#### counting number of cells for each age
age_counts <- controlCombined.umap.df %>%
  group_by(age, location) %>%
  summarise(count = n())
print(age_counts)


### plotting just the PCs
controlCombined.metrics.pca.df <- data.frame(controlCombined.metrics.pca$x, ID = seq_len(nrow(controlCombined.metrics.pca$x)))
controlCombined.metrics.pca.df <- inner_join(controlCombined.metrics.pca.df, controlCombined.meta, by = "ID")

plot = ggplot(controlCombined.metrics.pca.df, aes(x = PC1, y = PC2, color = location)) +
  geom_point() +
  theme_classic(base_size = 18) +
  #theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_color_manual(values = c(Arch = "#762A83", Desc = "#DFE725FF")) +
  labs(x = "PC1", y = "PC2")
#ggsave(plot, file = paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/plots/", '20230919_pca_location.svg'), width = 6, height = 4)
#ggsave(plot, file = paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/plots/", '20230919_pca_location.png'), width = 6, height = 4)
plot