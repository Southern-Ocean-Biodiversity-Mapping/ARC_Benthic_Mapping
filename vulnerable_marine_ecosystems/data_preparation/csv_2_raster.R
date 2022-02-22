################################################################################
#                                 TODO
################################################################################
# - Account for BIIGLE 254
# - Add other covariates: date, image quality, acqusition method

################################################################################
#                                 INIT
################################################################################
# Packages
#library(CCAMLRGIS)
#library(SOmap)
library(sp)
#library(proj4)
library(dplyr)

library(mapview)
library(raster)

# Set working directory
setwd("C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems")

# Path to BIO data
path_bio_data <- "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20220218161751_002/df_vme_idx.csv"

# Path to REF Raster
path_ref_raster <- "C:/Users/cgros/data/SO_env_layers/derived/Circumpolar_EnvData_500m_shelf_bathy_gebco_depth"

# Path to Coastline
#path_coastline <- "C:/Users/cgros/code/IMAS/ARC_Data/prep_environment/Circumpolar_Coastline.Rdata"
#load(path_coastline)

################################################################################
#                                 READ INPUT DATA
################################################################################

# Read Bio data
cat("\nReading data:", path_bio_data, "...\n")
df_idx <- read.csv(path_bio_data)

# Load REF raster
r_ref <- raster(path_ref_raster)

df_xy <- xyFromCell(r_ref, df_idx$cellID)
df_idx$x <- df_xy[, 1]
df_idx$y <- df_xy[, 2]

col_names <- colnames(df_idx)
col_not_coords <- col_names[col_names != "x" & col_names != "y"]

area_1 = list("xmin"=-3e6, "xmax"=-2e6,
              "ymin"=0, "ymax"=3e6)
area_2 = list("xmin"=-1e6, "xmax"=1e6,
              "ymin"=-2.5e6, "ymax"=-1e6)
area_3 = list("xmin"=1e6, "xmax"=3e6,
              "ymin"=-3e6, "ymax"=-1e6)
area_4 = list("xmin"=-1.5e6, "xmax"=1e6,
              "ymin"=0.5e6, "ymax"=2.5e6)
area = area_1
df_idx_area = df_idx[(df_idx$x > area$xmin) &
                       (df_idx$x < area$xmax) &
                       (df_idx$y > area$ymin) &
                       (df_idx$y < area$ymax), ]
min(df_idx_area$x)
max(df_idx_area$x)
min(df_idx_area$y)
max(df_idx_area$y)
# TODO: work on how to get all raster layers
r_idx_area <- rasterFromXYZ(df_idx_area[c("x", "y", "VME.index")], res=res(r_ref), crs=crs(r_ref))

# TODO: work on how to 
mapview(r_idx_area, maxpixels=9225970, na.color="transparent")













################################################################################
#                                RASTERIZE
################################################################################

r_bio = setValues(r_ref, NA)

r_bio[df_idx$cellID] = df_idx$VME.index
names(r_bio) <- "VME index"

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
