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
#library(plotly)


set.seed(23)

##defining the directories
getwd()


dataDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/output_metrics/"

plotDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/plots/20231221_finalPlots/"



##loading the data


data = as_tibble(read.table(paste0(dataDirectory, "20231221_allAges.csv"), stringsAsFactors=F, header = T, sep = ","))



#preprocessing


#removing extraneous columns, adding circularity

#artificial = data %>% select(-centroid.0,-centroid.1) %>% mutate(circularity = 4*3.14*area/(perimeter^2), type = "artificial", concavity_ct.log = log(concavity_ct+1), area_norm = 3.14*(major_axis_length*minor_axis_length)/2, dArea = area_norm-area)


#control = data %>% select(-centroid.0,-centroid.1, -"Unnamed..0")%>% mutate(concavity_ct.log = log(concavity_ct+1), area_norm = 3.14*(major_axis_length*minor_axis_length)/2, dArea = area_norm-area)
control = data %>% select(-centroid.0,-centroid.1)%>% mutate(concavity_ct.log = log(concavity_ct+1), area_norm = 3.14*(major_axis_length*minor_axis_length)/2, dArea = area_norm-area)

#allCombined.summary = allCombined %>% group_by(age, genotype, location) %>% summarise(sd = sd(metric), metric = mean(metric))
#allCombined.summary = allCombined %>% group_by(age, genotype, location) %>% summarise(sd = sd(rms_subsampled,na.rm = TRUE), rms_subsampled = mean(rms_subsampled,na.rm = TRUE))
#allCombined.summary = allCombined %>% group_by(age, genotype, location) %>% summarise(sd = sd(rms_avg, na.rm = TRUE), rms_avg = mean(rms_avg, na.rm = TRUE))


###########################################################################################################
#######################control hand-picked elliptical and weird nuclei ####################################
###########################################################################################################

#area cutoffs
lower_bound =  5509.125
upper_bound =  17812.125



controlCombined = control %>% filter(area >= lower_bound & area <= upper_bound & location != "") #added lower and upper bounds from 1.5 IQR area, calcualted in metrics jupyter notebook

metric = controlCombined$eccentricity
controlCombined.summary = controlCombined %>% group_by(location) %>% summarise(sd = sd(eccentricity, na.rm = TRUE), metric = mean(eccentricity, na.rm = TRUE))
#write.csv(controlCombined, "G:\\.shortcut-targets-by-id\\1VrZRYnU2bIYFYRrrUbM4YJeezJ09nRFb\\Melzer_nucleiDysmorphia\\data\\controlCombined.csv", row.names = FALSE)
plot1 = ggplot(controlCombined, aes(x = location, y = metric, color = location)) +
  geom_point(shape = 16, alpha = 0.6, size = 1.5, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5)) +
  #geom_jitter() +
  geom_violin(alpha = 0, size = 0.8, color = "black") +
  #geom_text_repel(aes(x=location, y= metric, label = label), color = "red", max.overlaps = 40) +
  theme_classic(base_size = 18) +
  geom_pointrange(aes(ymin = metric - sd, ymax = metric + sd), size = 1, position = position_dodge(width = 0.5), data = controlCombined.summary, color = "black") +
  theme(legend.position = "none", strip.background = element_blank(), strip.text.x = element_blank(), axis.title.x = element_blank(),  axis.title.y = element_blank()) +  # Remove y axis title
  scale_color_manual(values = c(Arch = "#762A83", Desc = "#DFE725FF"))
  #scale_color_manual(values = c(Arch = "#008B8B", Desc = "#E6E600"))
plot1
ggsave(plot1, file = paste0(plotDirectory, 'eccentricity.svg'), width = , height = 4)
ggsave(plot1, file = paste0(plotDirectory, 'eccentricity.png'), width = 4, height = 4)














#scatter plot
ggplot(controlCombined, aes(x = circularity, y = concavity_ct.log)) +
  geom_point(aes(color = type), shape = 16, size = 1) +
  theme_classic((base_size = 18)) +
  labs(x = "circularity", y = "concavity_ct.log") +
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())

## x-intercept count

xIntCombined.summary = xIntCombined %>% group_by(type) %>% summarise(sd = sd(signChanges, na.rm = TRUE), signChanges = mean(signChanges, na.rm = TRUE))

ggplot(xIntCombined, aes(x = type, y = signChanges)) +
  geom_point(aes(color = type), shape = 16, alpha = 0.3, size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5)) +
  geom_violin(aes(color = type), alpha = 0, linewidth = 0.8) +
  #geom_text_repel(aes(x=type, y=dArea, label = label), color = "red", max.overlaps = 40) +. ### cant use this unless also having label from control combined 
  theme_classic((base_size = 18)) +
  geom_pointrange(aes(ymin = signChanges - sd, ymax = signChanges + sd, color = type), size = 1, position = position_dodge(width = 0.5), data = xIntCombined.summary) +
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())

#ggsave(plot1, file = paste0(plotDirectory, 'spindleFraction.svg'), width = 6, height = 3)



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
                                                     convex_area,
                                                     x_intercepts,
                                                     concavity_ct,
                                                     circularity,
                                                     solidity,
                                                     #euler_number,
                                                     #form_factor,
                                                     equivalent_diameter_area,
                                                     area_bbox,
                                                     area_convex,
                                                     area_filled,
                                                     extent,
                                                     feret_diameter_max,
                                                     #minkowskiDimension,
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
                                                     convex_area,
                                                     x_intercepts,
                                                     concavity_ct,
                                                     circularity,
                                                     solidity,
                                                     #euler_number,
                                                     #form_factor,
                                                     equivalent_diameter_area,
                                                     area_bbox,
                                                     area_convex,
                                                     area_filled,
                                                     extent,
                                                     feret_diameter_max,
                                                     #minkowskiDimension,
                                                     compactness,
                                                     roundness,
                                                     convexity,
                                                     aspect_ratio,
                                                     concavity_ct.log,
                                                     perimeter_crofton) %>% mutate(ID = row_number())


controlCombined.metrics.noID = controlCombined.metrics %>% select(-ID)
controlCombined.metrics.noID = na.omit(controlCombined.metrics.noID)

controlCombined.metrics.pca = prcomp(controlCombined.metrics.noID, scale = TRUE)

controlCombined.metrics.pca$rotation = -1*controlCombined.metrics.pca$rotation

#controlCombined.pca.meta = as.data.frame(controlCombined.metrics.pca$x) %>% inner_join(as.data.frame(controlCombined.meta), by = character(0))

#write.csv(controlCombined.pca.meta, "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/output_metrics/UMAPs/20230918_age01_pca.csv")


controlCombined.pca.umap = umap(controlCombined.metrics.pca$x[,1:11])

controlCombined.umap.df = controlCombined.pca.umap$layout %>%
  as.data.frame() %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  #rename(UMAP1="X1", UMAP2="X2") %>%
  mutate(ID = row_number()) %>%
  inner_join(controlCombined.meta, by = "ID")

write.csv(controlCombined.umap.df, "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/output_metrics/UMAPs/20231221_umap_allAges.csv")


### making main UMAP plots

controlCombined.umap.df.plot = controlCombined.umap.df %>% filter(age == 24)
plot = ggplot(controlCombined.umap.df.plot, aes(x = UMAP1, y = UMAP2, color = as.factor(location))) +
  geom_point(size = 3) +
  #scale_color_continuous(name = "rms_avg") +  # Add a color legend with a name
  #theme_classic(base_size = 18) +
  theme_void() +
  theme(legend.position = "none") +
  #theme(legend.position = "right",  # Adjust the position of the legend
  #      strip.background = element_blank(),
  #      strip.text.x = element_blank()) +
  #geom_text_repel(aes(x=UMAP1, y=UMAP2, label = ID), color = "red", max.overlaps = 200) +
  scale_color_manual(values = c("Arch" = "#6b6ca3", "Desc" = "#c93f55"))
  #scale_color_viridis_d() + #for discrete colors
  #scale_color_viridis_c() + #for continuous colors
  #labs(x = "UMAP1", y = "UMAP2")
plot
ggsave(plot, file = paste0(plotDirectory, '52wk_dots.svg'), width = 6, height = 4)
ggsave(plot, file = paste0(plotDirectory, '52wk_dots.png'), width = 6, height = 4)

### density plots



controlCombined.umap.df.plot = controlCombined.umap.df %>% filter(age == 52)
plot = ggplot() +
  geom_point(data = controlCombined.umap.df.plot, aes(x = UMAP1, y = UMAP2), size = 1, color = "grey93") +
  geom_density_2d(data = controlCombined.umap.df.plot, aes(x = UMAP1, y = UMAP2, color = location), contour_var = "ndensity") +
  #scale_color_continuous(name = "rms_avg") +  # Add a color legend with a name
  #theme_classic(base_size = 18) +
  theme_void() +
  #labs(title = "12wk") +
  theme(legend.position = "none") +
  #theme(legend.position = "right",  # Adjust the position of the legend
  #      strip.background = element_blank(),
  #      strip.text.x = element_blank()) +
  #geom_text_repel(aes(x=UMAP1, y=UMAP2, label = ID), color = "red", max.overlaps = 200) +
  scale_color_manual(values = c("Arch" = "#6b6ca3", "Desc" = "#c93f55"))
  #scale_color_viridis_d() + #for discrete colors
  #scale_color_viridis_c() + #for continuous colors
  #labs(title = "52wk", x = "UMAP1", y = "UMAP2")
plot
ggsave(plot, file = paste0(plotDirectory, '52wk_density.svg'), width = 6, height = 4)
ggsave(plot, file = paste0(plotDirectory, '52wk_density.png'), width = 6, height = 4)



#### counting number of cells for each age

age_counts <- controlCombined.umap.df %>%
  group_by(age, location) %>%
  #group_by(location) %>%
  summarise(count = n())

print(age_counts)































ggplotly(plot) %>%
  layout(hovermode = "closest", hoverinfo = n)

## for minkowski dimension with legend for color!
ggplot(controlCombined.umap.df, aes(x = UMAP1, y = UMAP2, color = )) +
  geom_point() +
  scale_color_continuous(name = "area") +  # Add a color legend with a name
  theme_classic(base_size = 18) +
  theme(legend.position = "right",  # Adjust the position of the legend
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  #geom_text_repel(aes(x=UMAP1, y=UMAP2, label = ID), color = "red", max.overlaps = 200) +
  labs(x = "UMAP1", y = "UMAP2")

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

#### Plotting PCA plot colored by meta aspects of morphology

# Extract PCA scores
controlCombined.meta.naomit = na.omit(controlCombined.meta)
pca_scores <- predict(controlCombined.metrics.pca)

# Combine PCA scores with the "type" information from your original data
pca_df <- cbind.data.frame(pca_scores, age = controlCombined.meta$age)

# Plot using ggplot
ggplot(pca_df, aes(x = PC1, y = PC2, color = age)) +
  geom_point() +
  theme_classic(base_size = 18) +
  theme(legend.position = "none", strip.background = element_blank(), strip.text.x = element_blank()) +
  labs(x = "PC1", y = "PC2")



### elbow plot for pca above

explained_variance <- (controlCombined.metrics.pca$sdev)^2 / sum((controlCombined.metrics.pca$sdev)^2)

# Create a data frame with the principal component number and its explained variance
pca_data <- data.frame(PC = 1:length(explained_variance), Explained_Variance = explained_variance)

# Create the elbow plot
ggplot(pca_data, aes(x = PC, y = Explained_Variance)) +
  geom_line() +
  geom_point() +
  theme_classic(base_size = 18) +
  labs(x = "Principal Component", y = "Explained Variance")
#plot
#ggsave(plot, file = paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/plots/", '20230918_explainedVariance.png'), width = 6, height = 4)


## umap

umap_embeddings <- umap(controlCombined.metrics.noID)

plot_data <- data.frame(x = umap_embeddings$layout[, 1],
                        y = umap_embeddings$layout[, 2],
                        label = labels,
                        color = colors)

plot_data = umap_embeddings$layout %>%
  as.data.frame() %>%
  rename(UMAP1="V1", UMAP2="V2") %>%
  mutate(ID = row_number()) %>%
  inner_join(controlCombined.meta, by = "ID")

ggplot(plot_data, aes(x = x, y = y, label = ID, color = type)) +
  geom_point() +
  geom_text_repel() +
  theme_minimal()


plot(umap(controlCombined.metrics))


####### looking at PCAs to see how much variance is explained
# Compute cumulative variance
controlCombined.metrics.pca.variances <- controlCombined.metrics.pca$sdev^2
cumulative_variances <- cumsum(controlCombined.metrics.pca.variances) / sum(controlCombined.metrics.pca.variances)

plot(controlCombined.metrics.pca$sdev^2, type = "b", xlab = "Principal Component", ylab = "Eigenvalue")
plot(cumulative_variances, type = "b", xlab = "Number of Components", ylab = "Cumulative Proportion of Variance Explained")
screeplot(controlCombined.metrics.pca, type = "line", npcs = 10)

###OLD########################################################################################################################################################################################
###########################################################################################################################################################################################

#labeling for artificial nuclei:

#controlCombined.umap.df$rad <- factor(controlCombined.umap.df$rad)
#controlCombined.umap.df$n <- factor(controlCombined.umap.df$n)
#controlCombined.umap.df$edgy <- factor(controlCombined.umap.df$edgy)

#rad_colors <- c('0' = 'red', '0.5' = 'green', '1' = 'blue')
#n_colors <- c('4' = 'red', '9' = 'green', '15' = 'blue', '20' = 'purple')
#edgy_colors <- c('0' = 'red', '0.5' = 'green', '1' = 'blue')

#allCombined

allCombined.metrics = allCombined %>% select(eccentricity, perimeter, dArea, major_axis_length, minor_axis_length, rms_avg, x.intercepts, concavity_ct, circularity) %>% mutate(ID = row_number())

allCombined <- allCombined %>% mutate(ID=row_number())

allCombined.meta <- allCombined %>%
  select(ID, label, file, eccentricity, perimeter, dArea, major_axis_length, minor_axis_length, rms_avg, x.intercepts, concavity_ct, circularity, n, rad, edgy, id)

allCombined.metrics.noid = allCombined.metrics %>% select(-ID)

allCombined.metrics.pca = prcomp(allCombined.metrics.noid, scale = TRUE)

allCombined.metrics.pca$rotation = -1*allCombined.metrics.pca$rotation

allCombined.pca.umap = umap(allCombined.metrics.pca$x)

allCombined.umap.df <- allCombined.pca.umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID = row_number()) %>%
  inner_join(allCombined.meta, by = "ID")

ggplot(allCombined.umap.df, aes(x = UMAP1, y = UMAP2, color = age, shape = genotype)) +
  #geom_text_repel(aes(x=UMAP1, y=UMAP2, label = ID), color = "red", max.overlaps = 25) +
  geom_point(aes(color = age)) +
  theme_classic((base_size = 16)) +
  labs(x = "UMAP1",
       y = "UMAP2")


controlCombined.meta <- controlCombined %>% select(file, label, ID, n)

controlCombined.metrics = controlCombined %>% select(eccentricity,
                                                     perimeter,
                                                     dArea,
                                                     major_axis_length,
                                                     minor_axis_length,
                                                     rms_avg,
                                                     x.intercepts,
                                                     concavity_ct,
                                                     circularity,
                                                     solidity,
                                                     euler_number,
                                                     equivalent_diameter_area,
                                                     area_bbox,
                                                     area_convex,
                                                     area_filled,
                                                     extent,
                                                     feret_diameter_max,
                                                     minkowskiDimension,
                                                     perimeter_crofton) %>% mutate(ID = row_number())

controlCombined.metrics.noid = controlCombined.metrics %>% select(-ID)

controlCombined.metrics.pca = prcomp(controlCombined.metrics.noid, scale = TRUE)

controlCombined.umap <- controlCombined.metrics.pca %>%
  #select(where(is.numeric)) %>%
  #column_to_rownames("ID") %>%
  scale() %>%
  umap()

controlCombined.metrics.pca$rotation = -1*controlCombined.metrics.pca$rotation

controlCombined.pca.umap = umap(controlCombined.metrics.pca$x)

controlCombined.umap.df <- controlCombined.pca.umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID = row_number()) %>%
  inner_join(controlCombined.meta, by = "ID")

controlCombined.umap.df %>% head()

controlCombined.umap.df %>% ggplot(aes(x = UMAP1, y = UMAP2, color = n))+
  geom_point(aes(color = n)) +
  theme_classic((base_size = 18)) +
  labs(x = "UMAP1",
       y = "UMAP2")


###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################

ageList = allInducible %>% select(age) %>% unique() %>% mutate(ageReal = c("10-21", "10-30", "10-60","21-30","30-60","5-10","5-30", "5-60" ))

allInducibleAgeCorrected = inner_join(allInducible, ageList, by = "age") %>% select(-age) %>% rename(age = ageReal, distance = kNNType)

allInducibleAgeCorrected$age =  factor(allInducibleAgeCorrected$age, levels=c("5-10", "5-30", "5-60", "10-21", "10-30", "10-60", "21-30", "30-60"))
allInducibleAgeCorrected$aortaNumber =  factor(allInducibleAgeCorrected$aortaNumber, levels=c("aorta07", "aorta06", "aorta05", "aorta04", "aorta03", "aorta02", "aorta01"))

allInducibleAgeCorrected.summary = allInducibleAgeCorrected %>% group_by(age) %>% summarise(sd = sd(value), value = mean(value))

allInducible1NN.summary = allInducibleAgeCorrected %>% filter(distance == "1NN") %>% group_by(age) %>% summarise(sd = sd(value), value = mean(value))

plot1 = ggplot(allInducibleAgeCorrected %>% filter(distance == "1NN"), aes(x = age, y = value)) +
  geom_jitter(width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.5, color = "gray58") +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allInducible1NN.summary) +
  ylim(-2.5,2.5) +
  theme_classic((base_size = 18))+

ggsave(plot1, file = paste0(plotDirectory, 'inducibleRandomBaseline_knn1ALL.svg'), width = 6, height = 4)

plot1a = ggplot(allInducibleAgeCorrected %>% filter(distance == "1NN"), aes(x = age, y = value)) +
  geom_jitter(aes(color = aortaNumber), width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.5) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allInducible1NN.summary) +
  ylim(-2.5,2.5) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())

ggsave(plot1a, file = paste0(plotDirectory, 'inducibleRandomBaseline_knn1ALLColored.svg'), width = 6, height = 4)
ggsave(plot1a, file = paste0(plotDirectory, 'inducibleRandomBaseline_knn1ALLColoredLegend.svg'), width = 6, height = 4)


plot2 = ggplot(allInducibleAgeCorrected, aes(x = age, y = value)) +
  geom_jitter(color = "gray58", width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.4) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allInducibleAgeCorrected.summary) +
  ylim(-2.5,2.5) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())

ggsave(plot2, file = paste0(plotDirectory, 'inducibleRandomBaseline_knnAll.svg'), width = 6, height = 4)


plot2a = ggplot(allInducibleAgeCorrected, aes(x = age, y = value)) +
  geom_jitter(aes(color = distance), width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.4) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allInducibleAgeCorrected.summary) +
  ylim(-2.5,2.5) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())

ggsave(plot2a, file = paste0(plotDirectory, 'inducibleRandomBaseline_knnAllColored.svg'), width = 6, height = 4)
ggsave(plot2a, file = paste0(plotDirectory, 'inducibleRandomBaseline_knnAllColoredLegend.svg'), width = 6, height = 4)


###################
### Plotting Red and Orange channels on the same plot
allInducibleRed = as_tibble(read.table(paste0(dataDirectory, 'markedNormalized_kNN/', 'markedNormClusterskNN.csv'), stringsAsFactors=F, header = T, sep = ","))
allInducibleOrange = as_tibble(read.table(paste0(dataDirectory, 'markedNormalized_kNN/', 'markedNormClustersOrangekNN.csv'), stringsAsFactors=F, header = T, sep = ","))

ageListRed = allInducibleRed %>% select(age) %>% unique() %>% mutate(ageReal = c("10-21", "10-30", "10-60","21-30","30-60","5-10","5-30", "5-60" ))
ageListOrange = allInducibleOrange %>% select(age) %>% unique() %>% mutate(ageReal = c("10-21", "10-30", "5-10","5-30" ))

allInducibleAgeCorrectedRed = inner_join(allInducibleRed, ageListRed, by = "age") %>% select(-age) %>% rename(age = ageReal, distance = kNNType) %>% mutate(color = "red") %>% filter(age %in% c("10-21", "10-30", "5-10","5-30" ))
allInducibleAgeCorrectedOrange = inner_join(allInducibleOrange, ageListOrange, by = "age") %>% select(-age) %>% rename(age = ageReal, distance = kNNType) %>% mutate(color = "orange")

allInducibleAgeCorrectedOrangeRed = bind_rows(allInducibleAgeCorrectedRed,allInducibleAgeCorrectedOrange)
allInducibleAgeCorrectedOrangeRed$age =  factor(allInducibleAgeCorrectedOrangeRed$age, levels=c("5-10", "5-30", "10-21", "10-30"))

allInducibleAgeCorrectedOrangeRed.summary = allInducibleAgeCorrectedOrangeRed %>% group_by(age, color) %>% summarise(sd = sd(value), value = mean(value))
allInducibleAgeCorrectedOrangeRed1NN.summary = allInducibleAgeCorrectedOrangeRed %>% filter(distance == "1NN") %>% group_by(age,color) %>% summarise(sd = sd(value), value = mean(value))

plot3 = ggplot(allInducibleAgeCorrectedOrangeRed %>% filter(distance == "1NN"), aes(x = age, y = value)) +
  geom_jitter(aes(color = color), width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.5) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allInducibleAgeCorrectedOrangeRed1NN.summary) +
  ylim(-2.5,2.5) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot3, file = paste0(plotDirectory, 'inducibleRandomBaseline_knn1RedOrange.svg'), width = 3, height = 4)


plot4 = ggplot(allInducibleAgeCorrectedOrangeRed %>% filter(distance == "1NN"), aes(x = age, y = value)) +
  geom_point(aes(color = color), shape = 16, size = 0.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5)) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd, color = color), size = 1, position = position_dodge(width = 0.5), data = allInducibleAgeCorrectedOrangeRed1NN.summary) +
  ylim(-2.5,2.5) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot4, file = paste0(plotDirectory, 'inducibleRandomBaseline_knn1RedOrangeSideBySide.svg'), width = 5, height = 4)


###########################################################################################################
####################################### checking for division bias ########################################
###########################################################################################################

dataDirectory <- "/Volumes/GoogleDrive/My Drive/Northwestern/papers/PiBraunEtAl_growth/data/inducible rainbow/spatialPattern/kNN/"
plotDirectory <- '/Volumes/GoogleDrive/My Drive/Northwestern/papers/PiBraunEtAl_growth/plots/Supplementary/'

dataDirectory <- "/Users/yogeshgoyal/Library/CloudStorage/GoogleDrive-yogesh.goyal0308@gmail.com/My Drive/Northwestern/papers/PiBraunEtAl_growth/data/inducible rainbow/spatialPattern/kNN/"
plotDirectory <- '/Users/yogeshgoyal/Library/CloudStorage/GoogleDrive-yogesh.goyal0308@gmail.com/My Drive/Northwestern/papers/PiBraunEtAl_growth/plots/'


allDividingInducible = as_tibble(read.table(paste0(dataDirectory, 'dividingClusterskNN/', 'dividingClusterskNN_v3.csv'), stringsAsFactors=F, header = T, sep = ","))
allDividingInducible = allDividingInducible %>% mutate(type = "realData")

allDividingPos1 = as_tibble(read.table(paste0(dataDirectory, 'dividingClusterskNN/', 'dividedPositiveControlMidkNN1Run.csv'), stringsAsFactors=F, header = T, sep = ","))
allDividingPos2 = as_tibble(read.table(paste0(dataDirectory, 'dividingClusterskNN/', 'dividedPositiveControlOutkNN1Run.csv'), stringsAsFactors=F, header = T, sep = ","))
allDividingPos1= allDividingPos1 %>% mutate(type = "PositiveControl1")
allDividingPos2 = allDividingPos2 %>% mutate(type = "PositiveControl2")

allDividingNeg = as_tibble(read.table(paste0(dataDirectory, 'dividingClusterskNN/', 'dividedNegativeControlkNN1Run.csv'), stringsAsFactors=F, header = T, sep = ","))
allDividingNeg = allDividingNeg %>% mutate(type = "NegativeControl")


allDividingInducibleWithControls = bind_rows(allDividingInducible,allDividingPos1,allDividingPos2,allDividingNeg)

allDividingInducible.summary = allDividingInducible %>% group_by(age) %>% summarise(sd = sd(value), value = mean(value))
allInducible1NN.summary = allDividingInducible %>% filter(kNNType == "1NN") %>% group_by(age) %>% summarise(sd = sd(value), value = mean(value))

#allDividingInducibleWithControls.summary = allDividingInducibleWithControls %>% group_by(age,type) %>% summarise(sd = sd(value), value = mean(value))
allDividingInducibleWithControls.summary = allDividingInducibleWithControls %>% group_by(type) %>% summarise(sd = sd(value), value = mean(value))

#type                sd    value
#<fct>            <dbl>    <dbl>
#  1 NegativeControl  0.381  0.00273
#2 PositiveControl1 0.205 -0.446  
#3 PositiveControl2 0.290 -0.322  
# 4 realData         0.409  0.00791

allDividingInducibleWithControls$type = factor(allDividingInducibleWithControls$type, levels=c("realData", "NegativeControl","PositiveControl1", "PositiveControl2"))
allDividingInducibleWithControls.summary$type = factor(allDividingInducibleWithControls.summary$type, levels=c("realData", "NegativeControl","PositiveControl1", "PositiveControl2"))

allDividingInducible$age = factor(allDividingInducible$age, levels = c( "0-5", "0-10", "5-10","5-30", "5-60","10-21", "10-30", "10-60","21-30","30-60"))
allDividingInducible.summary$age = factor(allDividingInducible.summary$age, levels = c( "0-5", "0-10", "5-10","5-30", "5-60","10-21", "10-30", "10-60","21-30","30-60"))


plot5 = ggplot(allDividingInducible, aes(x = age, y = value)) +
  geom_jitter(aes(color = aorta), width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.4) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allDividingInducible.summary) +
  ylim(-4,4) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())

ggsave(plot5, file = paste0(plotDirectory, 'inducibleDividingRedAll.svg'), width = 8, height = 5)

set.seed(2059)

allDividingInducibleWithControls1 = sample_n(allDividingInducibleWithControls,25000) ## subsampling some rows otherwise the plots are pretty large size
allDividingInducibleWithControls.summary$type = factor(allDividingInducibleWithControls.summary$type, levels=c("realData", "NegativeControl","PositiveControl1", "PositiveControl2"))

real = 

wilcox.test(allDividingInducible$value, allDividingPos1$value ,
            alternative = "greater", paired = FALSE, conf.int = TRUE, conf.level = 0.95)

wilcox.test(allDividingInducible$value, allDividingPos2$value ,
            alternative = "greater", paired = FALSE, conf.int = TRUE, conf.level = 0.95)

wilcox.test(allDividingInducible$value, allDividingNeg$value ,
            alternative = "greater", paired = FALSE, conf.int = TRUE, conf.level = 0.95)


plot6 = ggplot(allDividingInducibleWithControls1, aes(x = type, y = value)) +
  geom_jitter(aes(color = type), width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.4) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allDividingInducibleWithControls.summary) +
  ylim(-4,4) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot6, file = paste0(plotDirectory, 'DividingWControls.svg'), width = 5, height = 4)

plot7 = ggplot(allDividingInducibleWithControls, aes(x = type, y = value)) +
  geom_jitter(aes(color = age), width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.4) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allDividingInducibleWithControls.summary) +
  ylim(-4,4) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot7, file = paste0(plotDirectory, 'DividingWControls_includingAgeColors.svg'), width = 5, height = 4)



# ####Do the next section only once to combine the files
# allInducible1NN = as_tibble(read.table(paste0(dataDirectory, 'markedNormalized_kNN/', 'markedNorm1NN.csv'), stringsAsFactors=F, header = T, sep = ",")) %>% add_column(distance = "1NN")
# allInducible2NN = as_tibble(read.table(paste0(dataDirectory, 'markedNormalized_kNN/', 'markedNorm2NN.csv'), stringsAsFactors=F, header = T, sep = ",")) %>% add_column(distance = "2NN")
# allInducible3NN = as_tibble(read.table(paste0(dataDirectory, 'markedNormalized_kNN/', 'markedNorm3NN.csv'), stringsAsFactors=F, header = T, sep = ",")) %>% add_column(distance = "3NN")
# allInducible4NN = as_tibble(read.table(paste0(dataDirectory, 'markedNormalized_kNN/', 'markedNorm4NN.csv'), stringsAsFactors=F, header = T, sep = ",")) %>% add_column(distance = "4NN")
# allInducible5NN = as_tibble(read.table(paste0(dataDirectory, 'markedNormalized_kNN/', 'markedNorm5NN.csv'), stringsAsFactors=F, header = T, sep = ",")) %>% add_column(distance = "5NN")
# 
# allInducible2130_3060 = as_tibble(read.table(paste0(dataDirectory, 'markedNormalized_kNN/', 'markedClusters_21-30_30-60_kNN.csv'), stringsAsFactors=F, header = T, sep = ",")) %>% rename(distance = kNNType) %>% select (-aortaNumber)
# 
# allInducibleALL= bind_rows(allInducible, allInducible2130_3060)
# allInducible = allInducibleALL
# 
# smallTable = allInducible %>% select(age,aorta) %>% unique()
# write.table(allInducible, file=paste0(dataDirectory, 'markedNormalized_kNN/', 'markedNormAllNNFINAL.csv'), col.names = TRUE, sep=',')
# write.table(smallTable, file=paste0(dataDirectory, 'markedNormalized_kNN/', 'smallTable.tsv'), col.names = TRUE, sep='/t')
# ####The previous section ends here and does not need to be repeated again
