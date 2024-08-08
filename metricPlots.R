####script for Salvador et al. Analysis of nuclear dysmorphia
####created by Madeline Melzer 20230731, updated by Madeline E Melzer 20240806
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
library(vegan)

set.seed(23)

### defining the directories
getwd()
dataDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/vimentin/results/metricTables/"
saveDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/vimentin/results/umapTables/"
plotDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/plots/"

### loading the data
raw_data = as_tibble(read.table(paste0(dataDirectory, "20240807_allAges.csv"), stringsAsFactors=F, header = T, sep = ","))

### applying area threshold
lower_bound =  5113
upper_bound =  18825
filtered_data = raw_data %>% filter(area >= lower_bound & area <= upper_bound & location != "" & !is.na(pointRadiusRMSDifference)) #added lower and upper bounds from 1.5 IQR area, cross check with metrics.ipynb result

# subset data to get the same number of cells for all 

sampled_data <- filtered_data %>%
  group_by(age, location, genotype) %>%
  slice_sample(n = 250, replace = FALSE) %>%
  ungroup()

### processing data for dimensional reduction with PCA and UMAP
data <- sampled_data %>% filter(genotype == "KO") %>% mutate(ID=row_number()) #the ID is to link the UMAP data with the metadata


metrics = data %>% select('area',
                          'perimeter',
                          'major_axis_length',
                          'minor_axis_length',
                          'eccentricity',
                          'solidity',
                          'equivalent_diameter_area',
                          'area_bbox',
                          'area_convex',
                          'extent',
                          'feret_diameter_max',
                          'perimeter_crofton',
                          'moments_hu.0',
                          'moments_hu.1',
                          'moments_hu.2',
                          'moments_hu.3',
                          'moments_hu.4',
                          'moments_hu.5',
                          'moments_hu.6',
                          'circularity',
                          'aspectRatio',
                          'flattening',
                          'averagePointRadius',
                          'pointRadiusExtremeDifference',
                          'pointRadiusRMSDifference',
                          'areaEquivalentEllipse',
                          'perimeterEquivalentEllipse',
                          'dArea',
                          'dPerimeter',
                          'dAreaConvex',
                          'xInterceptCount',
                          'concavityCount',
                          'ID')
                          
  
metadata = data %>% select('file',
                           'label',
                           'age',
                           'sex',
                           'animal',
                           'location',
                           'genotype', ## for Vimentin Knockout vs. WT only!
                           'area',
                           'perimeter',
                           'major_axis_length',
                           'minor_axis_length',
                           'eccentricity',
                           'solidity',
                           'equivalent_diameter_area',
                           'area_bbox',
                           'area_convex',
                           'extent',
                           'feret_diameter_max',
                           'perimeter_crofton',
                           'moments_hu.0',
                           'moments_hu.1',
                           'moments_hu.2',
                           'moments_hu.3',
                           'moments_hu.4',
                           'moments_hu.5',
                           'moments_hu.6',
                           'circularity',
                           'aspectRatio',
                           'flattening',
                           'averagePointRadius',
                           'pointRadiusExtremeDifference',
                           'pointRadiusRMSDifference',
                           'areaEquivalentEllipse',
                           'perimeterEquivalentEllipse',
                           'dArea',
                           'dPerimeter',
                           'dAreaConvex',
                           'xInterceptCount',
                           'concavityCount',
                           'ID')

metrics.noID = metrics %>% select(-ID)
metrics.noID = na.omit(metrics.noID)

metrics.pca = prcomp(metrics.noID, scale = TRUE)
metrics.pca$rotation = -1*metrics.pca$rotation

## if you want to save the PCA results with the metadata:
pca_scores <- as.data.frame(metrics.pca$x)
metadata_subset <- metadata %>% 
  select(file, label, age, sex, animal, location, ID) # Select only necessary metadata columns
data_pca_meta <- cbind(pca_scores, metadata_subset)# Combine PCA scores with metadata subset
write.csv(data_pca_meta, paste0(saveDirectory, "20240807_pca_vimentin_KO.csv"), row.names = FALSE)
#data_pca_meta = read.csv(paste0(saveDirectory, "20240807_pca_vimentin.csv"))


### looking at PCAs to see how much variance is explained
metrics.pca.variances <- metrics.pca$sdev^2
cumulative_variances <- cumsum(metrics.pca.variances) / sum(metrics.pca.variances) # Compute cumulative variance

plot(metrics.pca$sdev^2, type = "b", xlab = "Principal Component", ylab = "Eigenvalue")
plot(cumulative_variances, type = "b", xlab = "Number of Components", ylab = "Cumulative Proportion of Variance Explained")
screeplot(metrics.pca, type = "line", npcs = 20)


### choosing number of PCAs to use for UMAP based on screeplot for variance
data.pca.umap = umap(metrics.pca$x[,1:11])

data.umap.df = data.pca.umap$layout %>%
  as.data.frame() %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  #rename(UMAP1="X1", UMAP2="X2") %>%
  mutate(ID = row_number()) %>%
  inner_join(metadata, by = "ID")
write.csv(data.umap.df, paste0(saveDirectory, "20240807_umap_allMetrics_vimentin_KO.csv"))
#data.umap.df = read.csv(paste0(saveDirectory, "20240807_umap_allMetrics_vimentin.csv"))


### dot plots, UMAP
data.umap.df.plot = data.umap.df %>% filter(age == 8)
plot = ggplot(data.umap.df.plot, aes(x = UMAP1, y = UMAP2, color = location)) +
  geom_point(size = 2, shape = 16) +
  theme_void() +
  theme(legend.position = "none")+
  scale_color_manual(values = c("Arch" = "#6b6ca3", "Desc" = "#c93f55"))
plot
ggsave(plot, file = paste0(plotDirectory, 'dot_KO_8wk.svg'), width = 6, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'dot_KO_8wk.png'), width = 6, height = 4)

### density plots, UMAP
plot = ggplot() +
  geom_point(data = data.umap.df, aes(x = UMAP1, y = UMAP2), size = 1, shape = 16, color = "grey93") +
  geom_density_2d(data = data.umap.df.plot, aes(x = UMAP1, y = UMAP2, color = location), contour_var = "ndensity") +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Arch" = "#6b6ca3", "Desc" = "#c93f55"))
plot
ggsave(plot, file = paste0(plotDirectory, 'density_KO_8wk.svg'), width = 6, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'density_KO_8wk.png'), width = 6, height = 4)


### dot plots, UMAP, structured for genotype
data.umap.df.plot = data.umap.df %>% filter(age == 8)
plot = ggplot(data.umap.df.plot, aes(x = UMAP1, y = UMAP2, color = genotype)) +
  geom_point(size = 2, shape = 16) +
  theme_void() +
  theme(legend.position = "none")+
  scale_color_manual(values = c("WT" = "#6b6ca3", "KO" = "#c93f55"))
plot
#ggsave(plot, file = paste0(plotDirectory, 'dot_52wk.svg'), width = 6, height = 4)
#ggsave(plot, file = paste0(plotDirectory, 'dot_52wk.png'), width = 6, height = 4)

### density plots, UMAP, structured for genotype
plot = ggplot() +
  geom_point(data = data.umap.df, aes(x = UMAP1, y = UMAP2), size = 1, shape = 16, color = "grey93") +
  geom_density_2d(data = data.umap.df.plot, aes(x = UMAP1, y = UMAP2, color = genotype), contour_var = "ndensity") +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("WT" = "#6b6ca3", "KO" = "#c93f55"))
plot
#ggsave(plot, file = paste0(plotDirectory, 'density_52wk.svg'), width = 6, height = 4)
#ggsave(plot, file = paste0(plotDirectory, 'density_52wk.png'), width = 6, height = 4)




### counting number of cells for each age, aortic location (500)
age_counts <- filtered_data %>%
  group_by(age, location, genotype) %>%
  #group_by(location) %>%
  summarise(count = n())
print(age_counts)





