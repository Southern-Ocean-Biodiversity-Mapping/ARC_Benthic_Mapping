## DATA
# CCAMLR VME Registry
path_ccamlr_registry = "C:/Users/cgros/data/20210806_ccamlr_records/CCAMLR_VME_Registry_09032021.xlsx"
#path.raster = "C:/Users/cgros/data/20210806_ccamlr_records/Circumpolar_EnvData_bathy500m_shelf_gebco2020_depth.gri"
#resolution.raster = 3000

## METHOD PARAMETERS
# Reference system
proj_def = "+proj=longlat +datum=WGS84"
# Indicator taxa vulnerability scores
path_taxa_scores = "C:/Users/cgros/data/20210806_ccamlr_records/taxa_scores.csv"
# Aggregation method for vulnerability scores across taxa
agg_vulnerability_score = "mean"
# Number of categories for the abundance scores
n_abundance_categories = 5
# Aggregation method for VME indexes across taxa
vme_index_agg = "max"
# Number of categories for the VME indexes
n_index_categories = 3