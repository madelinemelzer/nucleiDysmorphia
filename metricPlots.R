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
dataDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/results/metricTables/"
saveDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/results/umapTables/"
plotDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/plots/"

### loading the data
raw_data = as_tibble(read.table(paste0(dataDirectory, "20240806_allAges.csv"), stringsAsFactors=F, header = T, sep = ","))

### applying area threshold
lower_bound =  5619.0
upper_bound =  18115.0
filtered_data = raw_data %>% filter(area >= lower_bound & area <= upper_bound & location != "" & !is.na(pointRadiusRMSDifference)) #added lower and upper bounds from 1.5 IQR area, cross check with metrics.ipynb result

# subset data to get the same number of cells for all 

sampled_data <- filtered_data %>%
  group_by(age, location) %>%
  slice_sample(n = 500, replace = FALSE) %>%
  ungroup()

### processing data for dimensional reduction with PCA and UMAP
data <- sampled_data %>% mutate(ID=row_number()) #the ID is to link the UMAP data with the metadata

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
#write.csv(data_pca_meta, paste0(saveDirectory, "20240806_pca.csv"), row.names = FALSE)
data_pca_meta = read.csv(paste0(saveDirectory, "20240806_pca.csv"))


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
#write.csv(data.umap.df, paste0(saveDirectory, "20240806_umap_allMetrics.csv"))
data.umap.df = read.csv(paste0(saveDirectory, "20240806_umap_allMetrics.csv"))


### dot plots, UMAP
data.umap.df.plot = data.umap.df %>% filter(age == 52)
plot = ggplot(data.umap.df.plot, aes(x = UMAP1, y = UMAP2, color = location)) +
  geom_point(size = 2, shape = 16) +
  theme_void() +
  theme(legend.position = "none")+
  scale_color_manual(values = c("Arch" = "#6b6ca3", "Desc" = "#c93f55"))
plot
ggsave(plot, file = paste0(plotDirectory, 'dot_52wk.svg'), width = 6, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'dot_52wk.png'), width = 6, height = 4)

### density plots, UMAP
plot = ggplot() +
  geom_point(data = data.umap.df, aes(x = UMAP1, y = UMAP2), size = 1, shape = 16, color = "grey93") +
  geom_density_2d(data = data.umap.df.plot, aes(x = UMAP1, y = UMAP2, color = location), contour_var = "ndensity") +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Arch" = "#6b6ca3", "Desc" = "#c93f55"))
plot
ggsave(plot, file = paste0(plotDirectory, 'density_52wk.svg'), width = 6, height = 4)
ggsave(plot, file = paste0(plotDirectory, 'density_52wk.png'), width = 6, height = 4)


### testing to see if groups are separate


# Function to perform PERMANOVA for one age group
perform_permanova <- function(age_group, data) {
  age_data <- data %>% filter(age == age_group)
  pca_coords <- age_data %>% select(starts_with("PC")) %>% select(1:11)  # Select first 11 PCs
  
  # Print some diagnostic information
  cat("Processing age group:", age_group, "\n")
  cat("Number of samples per location:\n")
  print(table(age_data$location))
  
  permanova_result <- adonis2(pca_coords ~ location, data = age_data, method = "euclidean")
  
  # Print the full PERMANOVA result
  print(permanova_result)
  
  return(data.frame(
    Age = age_group,
    F_value = permanova_result$F[1],
    R2 = permanova_result$R2[1],
    p_value = permanova_result$`Pr(>F)`[1]
  ))
}

# Apply PERMANOVA to each age group
age_groups <- unique(data_pca_meta$age)
permanova_results <- do.call(rbind, lapply(age_groups, perform_permanova, data = data_pca_meta))

# View results
print(permanova_results)



### ANOISM
perform_anosim <- function(age_group, data) {
  age_data <- data %>% filter(age == age_group)
  pca_coords <- age_data %>% select(starts_with("PC"))
  
  anosim_result <- anosim(pca_coords, age_data$location)
  
  return(data.frame(
    Age = age_group,
    R_statistic = anosim_result$statistic,
    p_value = anosim_result$signif
  ))
}

anosim_results <- do.call(rbind, lapply(age_groups, perform_anosim, data = data_pca_meta))
print(anosim_results)


### MORAN'S I
library(ape)
library(dplyr)

library(fields)  # for rdist function

perform_morans_i <- function(age_group, data) {
  age_data <- data %>% filter(age == age_group)
  pca_coords <- age_data %>% select(starts_with("PC")) %>% select(1:11)
  
  # Create a distance matrix based on location
  # 0 for same location, 1 for different location
  location_dist <- as.matrix(dist(as.numeric(as.factor(age_data$location)), method = "euclidean"))
  location_dist[location_dist > 0] <- 1
  location_dist <- 1 - location_dist  # Invert so that same location = 1, different = 0
  
  # Calculate distances between samples in PC space
  pc_dist <- as.matrix(dist(pca_coords))
  
  # Calculate Moran's I
  n <- nrow(pc_dist)
  W <- sum(location_dist)
  z <- pc_dist - mean(pc_dist)
  
  moran_i <- (n / W) * (sum(location_dist * z) / sum(z^2))
  
  # Calculate expected Moran's I
  expected_i <- -1 / (n - 1)
  
  # Calculate variance of Moran's I
  s2 <- n^2 * sum(location_dist^2) / (W^2 * (n-1)) - expected_i^2
  
  # Calculate z-score
  z_score <- (moran_i - expected_i) / sqrt(s2)
  
  # Calculate p-value
  p_value <- 2 * pnorm(-abs(z_score))
  
  return(data.frame(
    Age = age_group,
    Observed = moran_i,
    Expected = expected_i,
    Z_score = z_score,
    P_value = p_value
  ))
}

# Apply Moran's I to each age group
age_groups <- unique(data_pca_meta$age)
morans_i_results <- do.call(rbind, lapply(age_groups, perform_morans_i, data = data_pca_meta))

# View results
print(morans_i_results)








library(pairwiseAdonis)
pairwise_results <- lapply(age_groups, function(age_group) {
  age_data <- data_pca_meta %>% filter(age == age_group)
  pca_coords <- age_data %>% select(starts_with("PC"))
  
  pw_adonis <- pairwise.adonis(pca_coords, age_data$location)
  return(data.frame(Age = age_group, pw_adonis))
})

pairwise_results_df <- do.call(rbind, pairwise_results)
print(pairwise_results_df)






### plotting PC plot to check the spreads of the locations
pc_plot <- ggplot(data_pca_meta, aes(x = PC1, y = PC2, color = location)) +
  geom_point(alpha = 0.7) +  # Add some transparency to the points
  facet_wrap(~ age, scales = "free") +  # Create separate plots for each age
  theme_minimal() +
  labs(title = "PC1 vs PC2 by Age and Location",
       x = "PC1",
       y = "PC2",
       color = "Location") +
  theme(legend.position = "bottom")
pc_plot






### counting number of cells for each age, aortic location (500)
age_counts <- data.umap.df %>%
  group_by(age, location) %>%
  #group_by(location) %>%
  summarise(count = n())
print(age_counts)





