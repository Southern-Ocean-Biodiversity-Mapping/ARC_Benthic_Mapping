################################################################################
#                                 TODO
################################################################################
# - 

################################################################################
#                                 INIT
################################################################################
# Packages
library(CCAMLRGIS)
library(raster)
library(SOmap)
#library(sp)
#library(proj4)
library(ggplot2)
library(dplyr)
library(sp)
library(rgdal)
library(geosphere)

library(dismo)
library(tripack)
library(rgeos)

# Set working directory
setwd("C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/vme_index")

# Source functions
source("../../../Useful_Functions_Tools/SOmap_functions_JJ.R")
source("utils.R")

# Path towards data
path_raster <- "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/biodata_step4.csv"
path_metdata <- "C:/Users/cgros/code/IMAS/ARC_Data/annotation/Circumpolar_Annotation_Env_Data.RData"
################################################################################
#                                PROCESS
################################################################################
# Read Data
df_raster <- read.csv(path_raster)
list_cellID <- df_raster$cellID
load(path_metdata)
df <- cell_metadata_env[cell_metadata_env$cellID %in% list_cellID, ]

# Init map
JJ_SOmap()
#SOplot(x = df$longitude, y = df$latitude, pch = 19, col = 2)

# Load ASDs
ASDs <- load_ASDs()

# Create area
df$section_col <- "black"
df$section_name <- "out_of_scope"
# ASD 58.4.1
df_5841 <- get_points_in_asd(df_=df,
                            asd_object=ASDs,
                            asd_name="5841")
df[(df$cellID %in% df_5841$cellID) & (df$lon < 130), "section_col"] = "#00188f"
df[(df$cellID %in% df_5841$cellID) & (df$lon < 130), "section_name"] = "Mawson sea"
df[(df$cellID %in% df_5841$cellID) & (df$lon > 130), "section_col"] = "#bad80a"
df[(df$cellID %in% df_5841$cellID) & (df$lon > 130), "section_name"] = "D'Urville sea"

# ASD 881
df_881 <- get_points_in_asd(df_=df,
                            asd_object=ASDs,
                            asd_name="881")
df[df$cellID %in% df_881$cellID, "section_col"] = "#00b294"
df[df$cellID %in% df_881$cellID, "section_name"] = "Ross sea"
# ASD 882
df_882 <- get_points_in_asd(df_=df,
                            asd_object=ASDs,
                            asd_name="882")
df[df$cellID %in% df_882$cellID, "section_col"] = "#00b294" #"#bad80a"
df[df$cellID %in% df_882$cellID, "section_name"] = "Ross sea"
# ASD 481
df_481 <- get_points_in_asd(df_=df,
                            asd_object=ASDs,
                            asd_name="481")
df[df$cellID %in% df_481$cellID, "section_col"] = "#ec008c"
df[df$cellID %in% df_481$cellID, "section_name"] = "Bellinghausen sea"
# ASD 482
df_482 <- get_points_in_asd(df_=df,
                            asd_object=ASDs,
                            asd_name="482")
df[df$cellID %in% df_482$cellID, "section_col"] = "#ff8c00"
df[df$cellID %in% df_482$cellID, "section_name"] = "Scotia sea"
# ASD 485
df_485 <- get_points_in_asd(df_=df,
                            asd_object=ASDs,
                            asd_name="485")
df[df$cellID %in% df_485$cellID, "section_col"] = "#68217a"
df[df$cellID %in% df_485$cellID, "section_name"] = "Weddell sea"
# ASD 486
df_486 <- get_points_in_asd(df_=df,
                            asd_object=ASDs,
                            asd_name="486")
df[df$cellID %in% df_486$cellID, "section_col"] = "#e81123"
df[df$cellID %in% df_486$cellID, "section_name"] = "Lazarev sea"

df[(df$section_name == "out_of_scope") & (df$proj_coord_x > 0), "section_col"] = "#bad80a"
df[(df$section_name == "out_of_scope") & (df$proj_coord_x > 0), "section_name"] = "D'Urville sea"

df[(df$section_name == "out_of_scope") & (df$proj_coord_x < -2000000), "section_col"] = "#ec008c"
df[(df$section_name == "out_of_scope") & (df$proj_coord_x < -2000000), "section_name"] = "Bellinghausen sea"

df[(df$section_name == "out_of_scope") & (df$proj_coord_x > -700000) & (df$proj_coord_x < -400000), "section_col"] = "#e81123"
df[(df$section_name == "out_of_scope") & (df$proj_coord_x > -700000) & (df$proj_coord_x < -400000), "section_name"] = "Lazarev sea"

df[(df$section_name == "out_of_scope"), "section_col"] = "#68217a"
df[(df$section_name == "out_of_scope"), "section_name"] = "Weddell sea"

df[df$section_name == "out_of_scope", ]

df[(df$section_name == "Weddell sea") & (df$proj_coord_x < -2000000), "section_col"] = "#ec008c"
df[(df$section_name == "Weddell sea") & (df$proj_coord_x < -2000000), "section_name"] = "Bellinghausen sea"

# Plot
SOplot(x = df$lon, y = df$lat,  pch = 19, col = df$section_col, cex = 0.9)
JJ_SOleg(col = unique(df$section_col),
      type = "discrete",
      ladj = -0.1,
      tlabs = unique(df$section_name))



# Get sampling effort
df_ <- df[c("section_name", "counts_area", "section_col")]
dff <- as.data.frame(df_ 
                     %>% group_by(section_name) 
                     %>% summarise(across("counts_area",
                                          function(x,...)
                                            {if (!all(is.na(x))){sum(na.omit(x))} 
                                            else{NA}})))
dff$section_col <- df_$section_col[match(dff$section_name, df_$section_name)]
dff

# Plot Sampling effort
p <- dff %>%
  ggplot( aes(x=section_name, y=counts_area, fill=section_col)) +
  geom_bar(stat="identity", fill=dff$section_col, width=.4) +
  #geom_segment( aes(xend=section_name, yend=0)) +
  #geom_point( size=5, color=dff$section_col) +
  xlab("Subareas") +
  theme_bw() +
  ylab("Sampling effort (m2)") +
  theme(legend.position = "none")
p

df_out <- df[c("cellID", "lon", "lat", "proj_coord_x", "proj_coord_y", "section_col", "section_name")]

# convert data to a SpatialPointsDataFrame object
xy <- SpatialPointsDataFrame(
  df_out[, c("lon", "lat")], data.frame(cellID=df_out[, c("cellID")]),
  proj4string=CRS("+proj=longlat +datum=WGS84"))


# use the distm function to generate a geodesic distance matrix in meters
mdist <- distm(xy)

# cluster all points using a hierarchical clustering approach
hc <- hclust(as.dist(mdist), method="complete")

# define the distance threshold, in this case 40 m
d=50000

# define clusters based on a tree "height" cutoff "d" and add them to the SpDataFrame
xy$clust <- cutree(hc, h=d)

df_out$cluster <- xy$clust
for (i in 1:max(df_out$cluster)) {
  print(i)
  a <- df_out[df_out$cluster == i, "section_name"]
  print(unique(a))
  print(length(a))
  print(" ")
}


write.csv(df_out, "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/biodata_sea_area.csv", row.names = FALSE)










