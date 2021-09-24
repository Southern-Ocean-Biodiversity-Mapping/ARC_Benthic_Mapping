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

# Set working directory
setwd("C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/vme_index")

source("../../../Useful_Functions_Tools/SOmap_functions_JJ.R")

# Path towards data
path_raster <- "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20210920_raster_cover/raster_data.csv"
#source("../data_preparation/csv_to_raster.R")

################################################################################
#                                XX
################################################################################
df <- read.csv(path_raster)
JJ_SOmap()
#SOplot(x = df$longitude, y = df$latitude, pch = 19, col = 2)

ASDs <- load_ASDs()

# Create area
df$section_col <- "black"
df$section_name <- "out_of_scope"
# ASD 58.4.1
df_5841 <- get_points_in_asd(df_=df,
                            asd_object=ASDs,
                            asd_name="5841")
df[df$cellID %in% df_5841$cellID, "section_col"] = "#00188f"
df[df$cellID %in% df_5841$cellID, "section_name"] = "58.4.1"
# ASD 881
df_881 <- get_points_in_asd(df_=df,
                            asd_object=ASDs,
                            asd_name="881")
df[df$cellID %in% df_881$cellID, "section_col"] = "#00b294"
df[df$cellID %in% df_881$cellID, "section_name"] = "88.1"
# ASD 882
df_882 <- get_points_in_asd(df_=df,
                            asd_object=ASDs,
                            asd_name="882")
df[df$cellID %in% df_882$cellID, "section_col"] = "#bad80a"
df[df$cellID %in% df_882$cellID, "section_name"] = "88.2"
# ASD 481
df_481 <- get_points_in_asd(df_=df,
                            asd_object=ASDs,
                            asd_name="481")
df[df$cellID %in% df_481$cellID, "section_col"] = "#ec008c"
df[df$cellID %in% df_481$cellID, "section_name"] = "48.1"
# ASD 482
df_482 <- get_points_in_asd(df_=df,
                            asd_object=ASDs,
                            asd_name="482")
df[df$cellID %in% df_482$cellID, "section_col"] = "#ff8c00"
df[df$cellID %in% df_482$cellID, "section_name"] = "48.2"
# ASD 485
df_485 <- get_points_in_asd(df_=df,
                            asd_object=ASDs,
                            asd_name="485")
df[df$cellID %in% df_485$cellID, "section_col"] = "#68217a"
df[df$cellID %in% df_485$cellID, "section_name"] = "48.5"
# ASD 486
df_486 <- get_points_in_asd(df_=df,
                            asd_object=ASDs,
                            asd_name="486")
df[df$cellID %in% df_486$cellID, "section_col"] = "#e81123"
df[df$cellID %in% df_486$cellID, "section_name"] = "48.6"

SOplot(x = df$longitude, y = df$latitude,  pch = 19, col = df$section_col, cex = 0.9)

JJ_SOleg(col = unique(df$section_col),
      type = "discrete",
      ladj = -0.5,
      tlabs = unique(df$section_name))


df_ <- df[c("section_name", "area", "section_col")]
dff <- as.data.frame(df_ 
                     %>% group_by(section_name) 
                     %>% summarise(across("area",
                                          function(x,...)
                                            {if (!all(is.na(x))){sum(na.omit(x))} 
                                            else{NA}})))
dff$section_col <- df_$section_col[match(dff$section_name, df_$section_name)]

p <- dff %>%
  ggplot( aes(x=section_name, y=area, fill=section_col)) +
  geom_bar(stat="identity", fill=dff$section_col, width=.4) +
  #geom_segment( aes(xend=section_name, yend=0)) +
  #geom_point( size=5, color=dff$section_col) +
  xlab("Subareas") +
  theme_bw() +
  ylab("Sampling effort (m2)") +
  theme(legend.position = "none")
p
