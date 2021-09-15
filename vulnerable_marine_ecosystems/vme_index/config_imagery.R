## DATA
path_bio_data = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20210914_raster_cover/"

resolution_raster = 500

## METHOD PARAMETERS
# Reference system
#proj_def = "+proj=longlat +datum=WGS84"
# Indicator taxa vulnerability scores
path_taxa_scores = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/vme_index/morpho_taxa_scores.csv"
# Aggregation method for vulnerability scores across taxa
agg_vulnerability_score = "quadratic_mean"
# Number of categories for the abundance scores
n_abundance_categories = 3
# Aggregation method for VME indexes across taxa
vme_index_agg = "max"
# Number of categories for the VME indexes
n_index_categories = 3

## PACKAGES
library(raster)
library(RColorBrewer)
library(SOmap)
library(dplyr)
library(tidyr)
library(CCAMLRGIS)
library(BAMMtools)
