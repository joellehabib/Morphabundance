###########################################################################
# morphospace.r
#
# Purpose: Build and visualise a morphospace of detritus objects
#          (data import -> feature selection -> transform -> PCA -> clustering)
#
# Author:  Joelle Habib
# Date:    2025-11-15 (header added)
#
# Requires: R packages: morphr, chroma, castr, patchwork, lubridate,
#           tidyverse (dplyr, tidyr, ggplot2), imager, FactoMineR, factoextra,
#           wkmeans, vegan, gridExtra
#
# Inputs:  - `/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/float_ecotaxa.tsv` (tab-separated table, you get from ecotaxa)
#          - `/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/uvp6_images/` (image folders)
#
# Outputs: - Several PNG/JPEG plots saved in working directory (e.g. `Eig.jpg`, `Dim1.jpg`)
#          - CSVs with kmeans centers and PCA/morphotype results in
#            `GIT/TRATLEQ/Data/Float_data/15_OCT2025/`

###########################################################################

library("morphr")
library("chroma")
library("castr")
library("patchwork")
library("lubridate")
library("tidyverse")
library(dplyr)
library(imager)
library("FactoMineR")
library("factoextra")


# my data (equator)
ecotaxa_df <- read.table("/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/float_ecotaxa.tsv", sep = "\t", header = TRUE)

# only keep equator data 
#equator_df = filter(ecotaxa_df, sample_cruise == "WMO6903095")

# only keep ascent profiles 
equator_df <- subset(ecotaxa_df, grepl("a", ecotaxa_df$sample_id))

# select only detritus rows 
equator_detritus <- equator_df %>% 
  filter(grepl("detritus", object_annotation_hierarchy) | grepl("feces", object_annotation_hierarchy))

# set a base directory 
base_dir <- "/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/uvp6_images/"

# output directory for generated images/plots
out_dir <- "/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/15_OCT2025/"
# ensure output directory exists
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# get path to images
imgs <- str_c(base_dir, equator_detritus$sample_id, "/1/", equator_detritus$object_id, ".png") # I had to add '/1/' to my files path (I don't know why I have this supplementary folder in my export)
img_ok <- file.exists(imgs) # check if files exist 
sum(!img_ok) # nb of files that missing 
# TODO check why one file is missing -> 914 files missing for me ...

equator_detritus$img_path <- imgs
equator_detritus <- equator_detritus[img_ok,]

# remove the "object_" part of column name 
names(equator_detritus) <- gsub("object_", "", names(equator_detritus), fixed = TRUE)

# select traits
traits <- equator_detritus |>
  # keept only traits
  select(area:skeleton_area) |>
  # remove useless features
  select(-angle) |>
  # remove badly behaved/difficult to interpret features
  select(-starts_with("nb")) |>
  # remove clearly correlated ones
  select(-area_exc, -convarea, -skelarea, -symetriehc, -symetrievc, -`.area`, -convperim_perim,
         -feretareaexc, -perimareaexc) |>
  # select only variable that was not normal after transformation and difficult to interpret
  select(-circex, -fcons, -height, 
         -starts_with("histcum"), -kurt_mean, -max, -mode, -min, -minor, -perimferet, 
         -perimmajor, -slope, -width, -x, -xm, -y, -ym) #-area, -median, -convperim, -convarea_area) # TO DO --> ask to Jean-Olivier in which measure we consider a variable "normal"

# check no trait has variance 0
which(sapply(traits, var, na.rm=TRUE) < 10^-5)
# skew_mean had a variance smaller than 10^-5 so I remove it 
traits <- select(traits, -skew_mean)

# transform for normality
traits |>
  pivot_longer(everything()) |>
  ggplot() + facet_wrap(~name, scale="free") +
  geom_histogram(aes(x=value), bins=50)
ggsave(file.path(out_dir, "traits_raw_hist.png"), width = 12, height = 8, dpi = 150)

traits_trimmed <- traits |>
  mutate(
    area = mask_extreme(area, c(0, 0.1)),
    `circ.` = mask_extreme(`circ.`, c(0, 0.05)),
    convarea_area = mask_extreme(convarea_area, c(0, 0.03)),
    convperim = mask_extreme(convperim, c(0, 0.02)),
    cv = mask_extreme(cv, c(0, 0.05)),
    elongation = mask_extreme(elongation, c(0, 0.02)),
    feret = mask_extreme(feret, c(0, 0.02)),
    fractal = mask_extreme(fractal, c(0, 0.01)),
    intden = mask_extreme(intden, c(0, 0.03)),
    kurt = mask_extreme(kurt, c(0, 0.03)),
    major = mask_extreme(major, c(0, 0.03)),
    mean = mask_extreme(mean, c(0.01)),
    meanpos = mask_extreme(meanpos, c(0.01, 0)),
    median = mask_extreme(median, c(0.01, 0)),
    median_mean = mask_extreme(median_mean, c(0, 0)),
    median_mean_range = mask_extreme(median_mean_range, c(0, 0)),
    `perim.`= mask_extreme(`perim.`, c(0, 0.06)),
    range = mask_extreme(range, c(0, 0)),
    skeleton_area = mask_extreme(skeleton_area, c(0, 0.05)),
    skew = mask_extreme(skew, c(0.01, 0)),
    sr = mask_extreme(sr, c(0, 0)),
    stddev = mask_extreme(stddev, c(0, 0.01)),
    symetrieh = mask_extreme(symetrieh, c(0, 0.02)),
    symetrieh_area = mask_extreme(symetrieh_area, c(0, 0.01)),
    symetriev = mask_extreme(symetriev, c(0, 0.02)),
    symetriev_area = mask_extreme(symetriev_area, c(0, 0.01)),
    thickr = mask_extreme(thickr, c(0, 0.015)),
    # with all the variable 
    #mode = mask_extreme(mode, c(0, 0.05)),
    #max = mask_extreme(max, c(0, 0.05)),
    #min = mask_extreme(min, c(0, 0.05)),
    #x = mask_extreme(x, c(0, 0.05)),
    #y = mask_extreme(y, c(0, 0.05)),
    #xm = mask_extreme(xm, c(0, 0.05)),
    #ym = mask_extreme(ym, c(0, 0.05)),
    #width = mask_extreme(width, c(0, 0.05)),
    #height = mask_extreme(height, c(0, 0.05)),
    #minor = mask_extreme(minor, c(0, 0.05)),
    #slope = mask_extreme(slope, c(0, 0.05)),
    #histcum1 = mask_extreme(histcum1, c(0, 0.05)),
    #histcum2 = mask_extreme(histcum2, c(0, 0.05)),
    #histcum3 = mask_extreme(histcum3, c(0, 0.05)),
    #fcons = mask_extreme(fcons, c(0, 0.05)),
    #perimferet = mask_extreme(perimferet, c(0, 0.05)),
    #perimmajor = mask_extreme(perimmajor, c(0, 0.05)),
    #circex = mask_extreme(circex, c(0, 0.05)),
    #kurt_mean = mask_extreme(kurt_mean, c(0, 0.05))
  )

traits_trimmed |>
  pivot_longer(everything()) |>
  ggplot() + facet_wrap(~name, scale="free") +
  geom_histogram(aes(x=value), bins=50)
ggsave(file.path(out_dir, "traits_trimmed_hist.png"), width = 12, height = 8, dpi = 150)

traits_transfo <- mutate_all(traits_trimmed, function(x) {yeo_johnson(x) |> as.numeric()})
# NB: wrap yeo_johnson() in a function to avoid the extra attributes which mess up some fo the rest (but they are useful if we want to reproject in the same space)

traits_transfo |>
  pivot_longer(everything()) |>
  ggplot() + facet_wrap(~name, scale="free") +
  geom_histogram(aes(x=value), bins=50)
# -> ~OK. will do for now
ggsave(file.path(out_dir, "traits_transfo_hist.png"), width = 12, height = 8, dpi = 150)

# Set the seed for reproducibility
set.seed(123)

# Perform a random subsample of 5000 observations
subsample <- traits_transfo[sample(nrow(traits_transfo), 5000), ]
# test normality 
for (col in names(subsample)) {
  sw_test <- shapiro.test(subsample[[col]])
  p_value <- sw_test$p.value
  if (p_value < 0.05) {
    cat(sprintf("Variable %s is not normally distributed (p-value = %f)\n", col, p_value))
  } else {
    cat(sprintf("Variable %s is normally distributed (p-value = %f)\n", col, p_value))
  }
}

subsample |>
  pivot_longer(everything()) |>
  ggplot() + facet_wrap(~name, scale="free") +
  geom_histogram(aes(x=value), bins=50)
ggsave(file.path(out_dir, "subsample_hist.png"), width = 12, height = 8, dpi = 150)

# Reshape the data from wide to long format
traits_transfo_long <- gather(traits_transfo)

# Create a normal probability plot for each variable using ggplot2
ggplot(traits_transfo_long, aes(sample = value)) +
  stat_qq() +
  facet_wrap(~ key, ncol = 5) +
  theme_bw()

ms <- morphospace(traits_transfo, weights=rep(1, nrow(traits_transfo)))

ms$eig |> head(10)
# -> relevant until 3 or 4 

objs <- ms$ind$coord |> as_tibble()
ggplot(objs) + geom_point(aes(Dim.1, Dim.2), shape=".") + coord_fixed()
ggsave(file.path(out_dir, "objs_scatter_dim1_dim2.png"), width = 8, height = 6, dpi = 150)

## Save FactoMineR PCA variable plot (base graphics)
jpeg(file.path(out_dir, "PCA_FactoMineR_var.jpg"), width = 8, height = 6, units = 'in', res = 300)
FactoMineR::plot.PCA(ms, choix="var")
dev.off()


objs <- ms$ind$coord |> as_tibble()
ggplot(objs) +
  geom_point(aes(Dim.1, Dim.2), shape=".", alpha=0.5) +
  coord_fixed() +
  scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0)

# Function to pre-process UVP images
preprocess <- function(x) {
  x |>
    # remove 31 pixels from the bottom (=the scale bar)
    img_chop(bottom=31) |>
    # change the gamma to see the light objects better
    img_adjust_gamma(gamma=0.7)
}

# visual representation of the morphospace
ggmorph_tile(ms, equator_detritus$img_path,  steps=18, n_imgs=3, fun=preprocess)


objs <- ms$ind$coord |> as_tibble()


# PCs to plot
PCs <- 1:2


# get coordinates of active variables/columns
var <- ms$var$coord[,PCs] |> as.data.frame() |> rownames_to_column()

ggmorph_tile(ms, equator_detritus$img_path, steps=18, n_imgs=3, fun=preprocess) +
  #geom_point(data=objs, aes(Dim.1, Dim.2), shape=".", alpha=0.5) +
  geom_segment(data=var, aes(x = 0, y = 0, xend = (Dim.1*5),
                                     yend = (Dim.2*5)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black") +
  geom_point(size = 3) +
  geom_label(aes(x=6, y=-3, label="Brightness", fontface = 2), size = 7) +
  geom_label(aes(x=6, y=3, label="Size", fontface = 2), size = 7) +
  geom_label(aes(x=-6, y=3, label="Homogeneity", fontface = 2), size = 7) +
  geom_label(aes(x=-6, y=-3, label="Circularity", fontface = 2), size = 7) +
  theme(axis.title=element_text(size=14,face="bold"))
ggtitle("Title")



pca.vars <- ms$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

circ <- circleFun(c(0,0),2,npoints = 500)

ggmorph_tile(ms, equator_detritus$img_path, steps=18, n_imgs=3, fun=preprocess) 
  #geom_path(data=circ, aes(x*6,y*6), lty = 2, color = "grey", alpha = 0.7) +
  #geom_hline(yintercept=0, lty=2, color="grey", alpha=0.9) +
  #geom_vline(xintercept=0, lty=2, color="grey", alpha=0.9) +
  geom_segment(data=pca.vars, aes(x = 0, xend = Dim.1*6, y = 0, yend = Dim.2*6),
               arrow=arrow(length = unit(0.025, "npc"), type = "open"),
               lwd=1, color="seagreen") 
  #geom_text(data=pca.vars, aes(x=Dim.1*6.15, y=Dim.2*6.15, label=vars),
            #color="seagreen", size=4, hjust=0) #+ # labels

## Save the annotated morphospace (with variable vectors/labels)
p_annot <- ggmorph_tile(ms, equator_detritus$img_path, steps=18, n_imgs=3, fun=preprocess) +
  geom_segment(data=var, aes(x = 0, y = 0, xend = (Dim.1*5),
                             yend = (Dim.2*5)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black") +
  geom_point(size = 3) +
  geom_label(aes(x=6, y=-3, label="Brightness", fontface = 2), size = 7) +
  geom_label(aes(x=6, y=3, label="Size", fontface = 2), size = 7) +
  geom_label(aes(x=-6, y=3, label="Homogeneity", fontface = 2), size = 7) +
  geom_label(aes(x=-6, y=-3, label="Circularity", fontface = 2), size = 7) +
  theme(axis.title=element_text(size=14,face="bold")) +
  ggtitle("Annotated morphospace")

ggsave(file.path(out_dir, "morphospace_annotated.png"), plot = p_annot, width = 10, height = 8, dpi = 150)

## Save the PCA circle plot with variable arrows
p_pca_circle <- ggmorph_tile(ms, equator_detritus$img_path, steps=18, n_imgs=3, fun=preprocess) +
  geom_segment(data=pca.vars, aes(x = 0, xend = Dim.1*6, y = 0, yend = Dim.2*6),
               arrow=arrow(length = unit(0.025, "npc"), type = "open"),
               lwd=1, color="seagreen")

ggsave(file.path(out_dir, "morphospace_pca_vectors.png"), plot = p_pca_circle, width = 10, height = 8, dpi = 150)



#Visualize the eigenvalues and the percent explained
jpeg(file.path(out_dir, "Eig.jpg"), width = 6, height = 6, units = 'in', res = 600)
fviz_eig(ms, addlabels = TRUE, ylim = c(0, 70))
dev.off()

# top 10 variables contributing to the principal components:
# Contributions of variables to PC1
jpeg(file.path(out_dir, "Dim1.jpg"), width = 4, height = 4, units = 'in', res = 600)
fviz_contrib(ms, choice = "var", axes = 1, top = 10)
dev.off()
# Contributions of variables to PC2
jpeg(file.path(out_dir, "Dim2.jpg"), width = 4, height = 4, units = 'in', res = 600)
fviz_contrib(ms, choice = "var", axes = 2, top = 10)
dev.off()
jpeg(file.path(out_dir, "Dim3.jpg"), width = 4, height = 4, units = 'in', res = 600)
fviz_contrib(ms, choice = "var", axes = 3, top = 10)
dev.off()
jpeg(file.path(out_dir, "Dim4.jpg"), width = 4, height = 4, units = 'in', res = 600)
fviz_contrib(ms, choice = "var", axes = 4, top = 10)
dev.off()

## Clustering ----

# clustering of the points into many morphs
# clust_k <- kmeans(ms$ind$coord[,1:4], centers=200, nstart=20, iter_max=30)
library("wkmeans")
#clust_k <- wkmeans(ms$ind$coord[,1:4], k=200, nstart=100, iter_max=30)

# tr if it is different if put k = 4 directly 
clust_k <- wkmeans(ms$ind$coord[,1:4], k=5, nstart=100, iter_max=30)

centers <- as_tibble(clust_k$centers)

write.csv(centers, "/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/15_OCT2025/kmeans_centers_morpho.csv", row.names = FALSE)
write.csv(objs, "/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/15_OCT2025/morphotypes_df.csv", row.names = FALSE)

objs <- read_csv("/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/15_OCT2025/morphotypes_df.csv")

# plot the morphs in the space
ggplot(objs) +
  geom_point(aes(Dim.1, Dim.2, colour=factor(clust_k$cluster)), size = 2) +
  #geom_point(data=centers, aes(Dim.1, Dim.2), alpha=0.5, shape=16) + #add to add cluster centers
  coord_fixed() + scale_colour_discrete(guide="none") +
  theme_minimal()+
  xlab("PC1 (41.9%)") +
  ylab("PC2 (19.4%)") +
  scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
  theme(axis.title=element_text(size=14,face="bold"))


ggsave(file.path(out_dir, "morphospace_200_clusters.png"), width=10, height=8, dpi=100)


# hierarchisation of the morphs
#clust_h <- hclust(dist(clust_k$centers), method="ward.D2")
#plot(clust_h, hang=-1)
#dev.print(device=pdf, file="cluster_tree.pdf", width=8, height=6)
# -> no real clear cutting point
#    let's cut at 6

#clust_reduced <- cutree(clust_h, k=5) |> factor()

# assign each original point to a reduced cluster
#equator_detritus <- mutate(equator_detritus, cluster_kmeans=factor(clust_k$cluster), cluster_reduced=NA)
#match_clusters <- tibble(
 # cluster_kmeans=factor(1:nrow(clust_k$centers)),
  #cluster_reduced=clust_reduced
#)
#equator_detritus <- left_join(select(equator_detritus, -cluster_reduced), match_clusters)

objs$cluster_kmeans <- clust_k$cluster
equator_detritus$cluster_kmeans <- clust_k$cluster
### adding columns with dimensions to your dataframe:

equator_detritus$Dim.1 <- objs$Dim.1
equator_detritus$Dim.2 <- objs$Dim.2
equator_detritus$Dim.3 <- objs$Dim.3
equator_detritus$Dim.4 <- objs$Dim.4


### adding cluster names to your dataframe:
### Note each time numbers may differ!
equator_detritus$Cluster <- NA
equator_detritus$Cluster[equator_detritus$cluster_kmeans == 5] <- "cluster 5"
equator_detritus$Cluster[equator_detritus$cluster_kmeans == 4] <- "cluster 4"
equator_detritus$Cluster[equator_detritus$cluster_kmeans == 2] <- "cluster 2"
equator_detritus$Cluster[equator_detritus$cluster_kmeans == 3] <- "cluster 3"
equator_detritus$Cluster[equator_detritus$cluster_kmeans == 1] <- "cluster 1"



################################################################################
# try to compute distance between indiv and cluster center #
################################################################################

# Try to compute distance of each individual with its associated cluster center

center_clust1 <- centers[1,] 
#names(center_clust1) <- gsub(".", "", names(center_clust1), fixed = TRUE)
center_clust2 <- centers[2,] 
#names(center_clust2) <- gsub(".", "", names(center_clust2), fixed = TRUE)
center_clust3 <- centers[3,] 
#names(center_clust3) <- gsub(".", "", names(center_clust3), fixed = TRUE)
center_clust4 <- centers[4,] 
#names(center_clust4) <- gsub(".", "", names(center_clust4), fixed = TRUE)
center_clust5 <- centers[5,] 


indiv_clust1 <- objs %>% 
  filter(cluster_kmeans == 1) %>%
  select(Dim.1, Dim.2, Dim.3, Dim.4)

indiv_clust2 <- objs %>% 
  filter(cluster_kmeans == 2) %>%
  select(Dim.1, Dim.2, Dim.3, Dim.4)

indiv_clust3 <- objs %>% 
  filter(cluster_kmeans ==3) %>%
  select(Dim.1, Dim.2, Dim.3, Dim.4)

indiv_clust4 <- objs %>% 
  filter(cluster_kmeans == 4) %>%
  select(Dim.1, Dim.2, Dim.3, Dim.4)

indiv_clust5 <- objs %>% 
  filter(cluster_kmeans == 5) %>%
  select(Dim.1, Dim.2, Dim.3, Dim.4)


# compute euclidean distance between each individual coordinates and cluster center coordinates 
library(vegan)

# cluster 1 
distances1 <- dist(rbind(center_clust1, indiv_clust1), method = "euclidean")[2:(nrow(indiv_clust1) + 1)]

q25 <- quantile(distances1, 0.25)
#q95 <- quantile(distances1, 0.95)

indiv_clust1 <- indiv_clust1 %>% mutate(distance_to_center = distances1,
                                        mean_distance_to_center = mean(distances1), 
                                        q25 = q25)



# cluster 2
distances2 <- dist(rbind(center_clust2, indiv_clust2), method = "euclidean")[2:(nrow(indiv_clust2) + 1)]
q25_2 <- quantile(distances2, 0.25)
indiv_clust2 <- indiv_clust2 %>% mutate(distance_to_center = distances2,
                                        mean_distance_to_center = mean(distances2), 
                                        q25 = q25_2)

# cluster 3
distances3 <- dist(rbind(center_clust3, indiv_clust3), method = "euclidean")[2:(nrow(indiv_clust3) + 1)]
q25_3 <- quantile(distances3, 0.25)
indiv_clust3 <- indiv_clust3 %>% mutate(distance_to_center = distances3,
                                        mean_distance_to_center = mean(distances3),
                                        q25 = q25_3)

# cluster 4
distances4 <- dist(rbind(center_clust4, indiv_clust4), method = "euclidean")[2:(nrow(indiv_clust4) + 1)]
q25_4 <- quantile(distances4, 0.25)
indiv_clust4 <- indiv_clust4 %>% mutate(distance_to_center = distances4,
                                        mean_distance_to_center = mean(distances4), 
                                        q25 = q25_4)

# cluster 5
distances5 <- dist(rbind(center_clust5, indiv_clust5), method = "euclidean")[2:(nrow(indiv_clust5) + 1)]
q25_5 <- quantile(distances5, 0.25)
indiv_clust5 <- indiv_clust5 %>% mutate(distance_to_center = distances5,
                                        mean_distance_to_center = mean(distances5), 
                                        q25 = q25_5)

# compile all data set 
indiv_clust <- rbind(indiv_clust1, indiv_clust2) %>%
  rbind(., indiv_clust3)%>%
  rbind(., indiv_clust4)%>%
  rbind(., indiv_clust5)%>%
  select(Dim.1, distance_to_center, mean_distance_to_center, q25)

# add the distance_to_center column to the objs dataset 
objs <- merge(objs, indiv_clust, by = 'Dim.1')

objs <- objs %>% mutate(keep = distance_to_center <= q25)

# add the image path column to objs 
objs <- objs %>% select(Dim.1, distance_to_center, q25, keep)
equator_detritus <- merge(equator_detritus, objs, by = 'Dim.1')


objs_keep <- equator_detritus %>% filter(keep == TRUE)
objs_delete <- equator_detritus %>% filter(keep == FALSE)

# plot the station in the space
ggplot() +
  geom_point(aes(x = Dim.1, y = Dim.2, colour=factor(cluster_kmeans)), data  = objs_keep, size = 2)+
  geom_point(aes(x = Dim.1, y = Dim.2, colour=factor(cluster_kmeans)), data  = objs_delete, alpha = 0.05, size = 2) +
  #geom_label(aes(x = Dim.1, y = Dim.2, label = CYCLE_NUMBER), data = objs_delete, color = as.factor(objs_delete$cluster), alpha = 0.1)+
  geom_point(data=centers, aes(Dim.1, Dim.2), alpha=0.5, shape=16, size = 4, color = 'yellow') + #add to add cluster centers
  coord_fixed() + scale_colour_discrete(guide="none") +
  theme_minimal()+
  xlab("PC1 (41.9%)") +
  ylab("PC2 (19.4%)") +
  scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
  theme(axis.title=element_text(size=14,face="bold"))

ggplot() +
  geom_point(aes(x = Dim.1, y = Dim.2, colour=factor(cluster_kmeans)), data  = objs_keep, size = 2)+
  geom_point(aes(x = Dim.1, y = Dim.2, colour=factor(cluster_kmeans)), data  = objs_delete, alpha = 0.2, size = 2) +
  #geom_label(aes(x = Dim.1, y = Dim.2, label = CYCLE_NUMBER), data = objs_delete, color = as.factor(objs_delete$cluster), alpha = 0.1)+
  geom_point(data=centers, aes(Dim.1, Dim.2), alpha=1, shape=16, size = 4, color = 'yellow') + #add to add cluster centers
  coord_fixed() + scale_colour_manual('',values=c("#00BFC4", "#F8766D", "#CD8BFF", "#7AAF00","#CD8BFF")) +
  theme_minimal()+
  xlab("PC1 (41.9%)") +
  ylab("PC2 (19.4%)") +
  scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
  theme(axis.title=element_text(size=14,face="bold"), 
        legend.position = "none")

ggplot() +
  geom_point(aes(x = Dim.1, y = Dim.3, colour=factor(cluster_kmeans)), data  = objs_keep, size = 2)+
  geom_point(aes(x = Dim.1, y = Dim.3, colour=factor(cluster_kmeans)), data  = objs_delete, alpha = 0.2, size = 2) +
  #geom_label(aes(x = Dim.1, y = Dim.2, label = CYCLE_NUMBER), data = objs_delete, color = as.factor(objs_delete$cluster), alpha = 0.1)+
  geom_point(data=centers, aes(Dim.1, Dim.3), alpha=1, shape=16, size = 4, color = 'yellow') + #add to add cluster centers
  coord_fixed() + scale_colour_manual('',values=c("#00BFC4", "#F8766D", "#CD8BFF", "#7AAF00","#CD8BFF")) +
  theme_minimal()+
  xlab("PC1 (41.9%)") +
  ylab("PC3 (17.5%)") +
  scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
  theme(axis.title=element_text(size=14,face="bold"), 
        legend.position = "none")

ggplot() +
  geom_point(aes(x = Dim.2, y = Dim.3, colour=factor(cluster_kmeans)), data  = objs_keep, size = 2)+
  geom_point(aes(x = Dim.2, y = Dim.3, colour=factor(cluster_kmeans)), data  = objs_delete, alpha = 0.2, size = 2) +
  #geom_label(aes(x = Dim.1, y = Dim.2, label = CYCLE_NUMBER), data = objs_delete, color = as.factor(objs_delete$cluster), alpha = 0.1)+
  geom_point(data=centers, aes(Dim.2, Dim.3), alpha=1, shape=16, size = 4, color = 'yellow') + #add to add cluster centers
  coord_fixed() + scale_colour_manual('',values=c("#00BFC4", "#F8766D", "#CD8BFF", "#7AAF00","#CD8BFF")) +
  theme_minimal()+
  xlab("PC2 (19.4%)") +
  ylab("PC3 (17.5%)") +
  scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
  theme(axis.title=element_text(size=14,face="bold"), 
        legend.position = "none")


################################################################################
#
################################################################################


ggplot(equator_detritus) +
  geom_point(aes(Dim.1, Dim.2, colour=as.factor(cluster_kmeans)), size = 2) +
  coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
  scale_colour_discrete("Cluster", guide=guide_legend(override.aes=list(shape=16)))
ggsave(file.path(out_dir, "morphospace_clusters_reduced.png"), width=10, height=8, dpi=100)




## Display and analysis of clusters -----

# get some example images from each reduced cluster
clust_imgs <- equator_detritus |>
  group_by(cluster_kmeans) |>
  sample_n(size=min(n(), 100)) |>
  group_map(.f=function(x, y) {
    ggimg_grid(x$img_path, fun=preprocess, scale=0.002) +
      labs(title=str_c("Cluster ", y$cluster_kmeans)) +
      theme(plot.title=element_text(hjust=0.5))
  })
wrap_plots(clust_imgs, ncol=2)
ggsave(file.path(out_dir, "cluster_images.png"), width=10, height=6, dpi=200)

# compare with only selected individuals (distance to center =< mean(distance to center)) 
clust_imgs <- objs_keep |>
  group_by(cluster_kmeans) |>
  sample_n(size=min(n(), 100)) |>
  group_map(.f=function(x, y) {
    ggimg_grid(x$img_path, fun=preprocess, scale=0.002) +
      labs(title=str_c("Cluster ", y$cluster_kmeans)) +
      theme(plot.title=element_text(hjust=0.5))
  })
wrap_plots(clust_imgs, ncol=2)




# Box plot of the size, shape, brightness and structure of each cluster 

# size ----

size <- equator_detritus %>% select(perim., Cluster)
size$perim. <- log10(size$perim.)

cluster_colors <- c("cluster 1" = "#DE64B1",
                    "cluster 2" = "#57A7B3",
                    "cluster 3" = "#F5F28E",
                    "cluster 4" = "grey33",
                    "cluster 5" = "#A6DAE6")

size$Cluster <- factor(size$Cluster, levels = c("cluster 1", "cluster 2", "cluster 3", "cluster 4", "cluster 5"))

p1 <- ggplot(size, aes(x = Cluster, y = perim., fill = Cluster))+
  geom_boxplot() +
  scale_fill_manual(values = cluster_colors) +
  ggtitle("Size") +
  xlab("Cluster Name") +
  ylab("log(Perimeter)") +
  theme_bw()+
  theme(legend.position = "none",
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15))
p1

# shape ----


shape <- equator_detritus %>% select(circ., Cluster)
shape$circ. <- shape$circ.

shape$Cluster <- factor(shape$Cluster, levels = c("cluster 1", "cluster 2", "cluster 3", "cluster 4", "cluster 5"))

p2 <- ggplot(shape, aes(x = Cluster, y = circ., fill = Cluster))+
  geom_boxplot() +
  scale_fill_manual(values = cluster_colors) +
  ggtitle("Shape") +
  xlab("Cluster Name") +
  ylab("Circularity") +
  theme_bw()+
  theme(legend.position = "none",
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15))

p2

# Brightness ----

brightness <- equator_detritus %>% select(mean, Cluster)

brightness$Cluster <- factor(brightness$Cluster, levels = c("cluster 1", "cluster 2", "cluster 3", "cluster 4", "cluster 5"))

p3 <- ggplot(brightness, aes(x = Cluster, y = mean, fill = Cluster))+
  geom_boxplot() +
  scale_fill_manual(values = cluster_colors) +
  ggtitle("Brightness") +
  xlab("Cluster Name") +
  ylab("mean grey level") +
  theme_bw()+
  theme(legend.position = "none",
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15))
p3

# Structure ----

structure <- equator_detritus %>% select(kurt, Cluster)

structure$Cluster <- factor(structure$Cluster, levels = c("cluster 1", "cluster 2", "cluster 3", "cluster 4", "cluster 5"))

p4 <- ggplot(structure, aes(x = Cluster, y = kurt, fill = Cluster))+
  geom_boxplot() +
  scale_fill_manual(values = cluster_colors) +
  ggtitle("Structure") +
  xlab("Cluster Name") +
  ylab("kurtosis") +
  theme_bw()+
  theme(legend.position = "none",
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15))

p4


# put the four plot into one figure ----
library(gridExtra)
## save the single size boxplot and a 4-panel grid to the output folder
ggsave(file.path(out_dir, "size_boxplot.png"), plot = p1, width = 5, height = 6, dpi = 200)

# create a 4-panel grob and save it
grid_plot <- arrangeGrob(p1, p2, p3, p4, ncol = 4)
ggsave(file.path(out_dir, "boxplots_grid.png"), plot = grid_plot, width = 14, height = 6, dpi = 200, limitsize = FALSE)

# time series based on PCA results ----

# select only columns that I'm interested in 

equator_detritus$datetime <- paste(equator_detritus$date, equator_detritus$time, sep = "")
equator_detritus <- equator_detritus %>% select(-date, -time)
equator_detritus$datetime <- as.POSIXct(equator_detritus$datetime, format = "%Y%m%d%H%M%S")

ts_df <-  equator_detritus %>% select(id, datetime, depth_min, Dim.1, Dim.2, Dim.3, Dim.4, Cluster, keep)

# save it 

write.csv(ts_df, '/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/15_OCT2025/k-means_pca_results_all.csv', row.names = FALSE)
write.csv(ts_df, '/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/15_OCT2025/k-means_pca_results_q25.csv', row.names = FALSE)
