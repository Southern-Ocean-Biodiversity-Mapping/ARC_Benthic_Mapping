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

# Source functions
source("utils.R")

# Path towards data
path_raster <- "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20210927_raster_cover/raster_data_ccamlr.csv"
#path_raster <- "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20210927_raster_cover/raster_data.csv"
#path_scores <- "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/vme_index/morpho_taxa_scores.csv"
path_scores <- "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/vme_index/taxa_scores.csv"
################################################################################
#                                PROCESS
################################################################################
# Read Data
df <- read.csv(path_raster)
df_taxa <- read.csv(path_scores)
if ("morpho_taxon" %in% colnames(df_taxa)) {
  list_taxa_all <- gsub("-", "..", df_taxa$morpho_taxon)
} else {
  list_taxa_all <- NA
  #list_taxa_all <- tolower(df_taxa$Taxon)
}

# Load ASDs
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

# Select columns
if (is.na(list_taxa_all)) {
  list_taxa_selected <- colnames(df)[!(colnames(df) %in% c("n_annotation", "area_pix", "area", "cellID", "survey", "longitude", "latitude", "filename", "proj_coord_x", "proj_coord_y", "section_col", "section_name"))]
} else {
  list_taxa_selected <- intersect(list_taxa_all, colnames(df))
}
df_ <- df[, c(list_taxa_selected, "section_name", "section_col")]





# Get statistics
df_mean <- as.data.frame(df_ 
                     %>% group_by(section_name) 
                     %>% summarise(across(list_taxa_selected,
                                          function(x,...)
                                          {if (!all(is.na(x))){mean(na.omit(x))} 
                                            else{NA}})))
df_sd <- as.data.frame(df_ 
                           %>% group_by(section_name) 
                           %>% summarise(across(list_taxa_selected,
                                                function(x,...)
                                                {if (!all(is.na(x))){sd(na.omit(x))} 
                                                  else{NA}})))

df_mean <- as.data.frame(df_ 
                         %>% summarise(across(list_taxa_selected,
                                              function(x,...)
                                              {if (!all(is.na(x))){mean(na.omit(x))} 
                                                else{NA}})))

df_sd <- as.data.frame(df_ 
                       %>% summarise(across(list_taxa_selected,
                                            function(x,...)
                                            {if (!all(is.na(x))){sd(na.omit(x))} 
                                              else{NA}})))

df_violin <- as.data.frame(df_ %>% tidyr::pivot_longer(
  cols = all_of(list_taxa_selected),
  names_to = "taxa", 
  values_to = "percent_cover"))
#df_violin$morpho_taxa <- gsub("\\.", "-", df_violin$morpho_taxa)
#df_violin <- separate(data = df_violin, col = morpho_taxa, into = c("morpho", "taxa"), sep = "--")

# Plot Sampling effort
p <- df_violin %>%
  #left_join(sample_size) %>%
  #mutate(myaxis = paste0(name, "\n", "n=", num)) %>%
  ggplot( aes(x=taxa, y=percent_cover, fill=taxa)) +
  geom_violin(width=1.4) +
  #geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  #scale_fill_viridis(discrete = TRUE) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A Violin wrapping a boxplot") +
  ylim(c(0, 5)) +
  xlab("")
p



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
