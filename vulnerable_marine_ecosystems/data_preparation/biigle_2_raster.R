################################################################################
#                                 TODO
################################################################################
# Take empty images into account

################################################################################
#                                 INIT
################################################################################
# Packages
library(stringr)

# Set working directory
setwd("C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/data_preparation")

# Path towards BIIGLE report
path_biigle <- "C:/Users/cgros/code/biigle_scripts/20210908_biigle_vme.csv"


################################################################################
#                                 XX
################################################################################
# Read BIIGLE data
df_biigle <- read.csv(path_biigle)
df_biigle <- df_biigle[, c("label_hierarchy",
                        "filename",
                        "image_longitude",
                        "image_latitude",
                        "attributes",
                        "points",
                        "shape_name")]
# Clean label names
df_biigle$label <- gsub(" > ", "_", df_biigle$label_hierarchy)
# Add survey column
df_biigle$survey <- str_split_fixed(df_biigle$filename, "_", 2)[, 1]
# Get image width, height, area
df_biigle$width <- as.integer(str_split_fixed(str_split_fixed(df_biigle$attributes,
                                                              '"width":', 2)[, 2],
                                              ',', 2)[, 1])
df_biigle$height <- as.integer(str_split_fixed(str_split_fixed(df_biigle$attributes,
                                                              '"height":', 2)[, 2],
                                              ',', 2)[, 1])
df_biigle$area <- as.double(str_split_fixed(str_split_fixed(df_biigle$attributes,
                                                               '"area":', 2)[, 2],
                                               ',', 2)[, 1])
# Fill missing values of width and height




# Remove unnecessary columns
df_biigle <- df_biigle[,!names(df_biigle) %in% c("label_hierarchy", "filename")]

colnames(dataframe)[which(names(dataframe) == "columnName")] <- "newColumnName"


