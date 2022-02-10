################################################################################
#                                 TODO
################################################################################
# - Account for BIIGLE 254
# - Add other covariates: date, image quality, acqusition method

################################################################################
#                                 INIT
################################################################################
# Packages
library(CCAMLRGIS)
library(raster)
library(SOmap)
library(sp)
library(proj4)
library(dplyr)

# Set working directory
setwd("C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems")

# Path to BIO data
path_bio_data <- "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/biodata_step4.csv"

# Path to REF Raster
path_ref_raster <- "C:/Users/cgros/data/SO_env_layers/derived/Circumpolar_EnvData_500m_shelf_bathy_gebco_depth"

# Path to Coastline
path_coastline <- "C:/Users/cgros/code/IMAS/ARC_Data/prep_environment/Circumpolar_Coastline.Rdata"
load(path_coastline)

################################################################################
#                                 READ INPUT DATA
################################################################################

# Read Bio data
cat("\nReading data:", path_bio_data, "...\n")
df_abd <- read.csv(path_bio_data)

# Create PA data
lst_taxa <- colnames(df_abd)
lst_taxa <- lst_taxa[!lst_taxa == "cellID"]
df_pa <- data.frame(df_abd)
for (t in lst_taxa) {
  df_pa[ , t] <- as.integer(as.logical(df_pa[ , t]))
}

# Load REF raster
r_ref <- raster(path_ref_raster)

################################################################################
#                                RASTERIZE
################################################################################

r_bio = setValues(r_ref, NA)

r_bio[df_abd$cellID] = df_abd[[lst_taxa[1]]]
names(r_bio) <- lst_taxa[1]

for (idx_taxa in 2:4) {
  print(lst_taxa[idx_taxa])
  r_tmp = setValues(r_ref, NA)
  r_tmp[df_abd$cellID] = df_abd[[lst_taxa[idx_taxa]]]
  names(r_tmp) <- lst_taxa[idx_taxa]
  r_bio <- stack(r_bio, r_tmp)
  remove(r_tmp)
}

buffer_size <- 10000
rnd_cell <- sample(1:nrow(df_abd), 1)
coords_zoom_center <- xyFromCell(r_ref, df_abd$cellID[[rnd_cell]])
df_zoom <- data.frame(x = coords_zoom_center[[1]], y = coords_zoom_center[[2]])
coordinates(df_zoom) <- c("x", "y")
sf_zoom <- st_as_sf(df_zoom, crs = crs(r_ref))
buffer <- st_buffer(sf_zoom, buffer_size)
r_crop <- crop(r_bio, extent(buffer))
plot(r_crop)


r[is.na(r[])] <- 0 

library(mapview)
library(leaflet)
library(leafem)

leaflet() %>% 
  addRasterImage(r_ref, layerId = "values") %>% 
  addMouseCoordinates() %>%
  addImageQuery(r_ref, type="mousemove", layerId = "values")
