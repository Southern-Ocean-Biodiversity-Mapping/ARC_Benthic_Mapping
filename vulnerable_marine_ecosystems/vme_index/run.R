# Set working directory
setwd("C:/Users/cgros/code/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/vme_index")

# Load config file
readRenviron("config.R")
path.ccamlr.registry = Sys.getenv("path.ccamlr.registry")
path.taxa.scores = Sys.getenv("path.taxa.scores")
resolution.raster = as.numeric(Sys.getenv("resolution.raster"))
score.agg = Sys.getenv("score.agg")
n_abundance_categories = as.numeric(Sys.getenv("n_abundance_categories"))
vme_index_agg = Sys.getenv("vme_index_agg")
n_index_categories = as.numeric(Sys.getenv("n_index_categories"))

# Load functions
source("utils.R")

# Packages
library(raster)
library(RColorBrewer)
library(SOmap)
library(dplyr)
library(tidyr)
library(CCAMLRGIS)
library(BAMMtools)

# Read CCAMLR registry
ccamlr.registry <- read_excel_allsheets(path.ccamlr.registry)
# Curate data
ccamlr.data <- curate_ccamlr_registry(ccamlr.registry)
# Read VME Taxa Scores
scores.data <- read.csv(path.taxa.scores)

# Compute VME indicator score
if (score.agg == "mean") {
  scores.data$score <- rowMeans(subset(scores.data, select = 5:11), na.rm = TRUE)
}
# TODO: other method eg quadratic mean
# Order for visualisation
scores.data <- scores.data[order(scores.data$score, decreasing = TRUE),]
print(scores.data[ , c("Taxon", "score")])

# Get VME Risk Areas
df.risk.areas <- ccamlr.data[["vme.risk.areas"]]
df.risk.areas.taxa <- ccamlr.data[["vme.risk.areas.taxa"]]
df.risk.areas.taxa <- df.risk.areas.taxa %>% separate(VMEIndicatorTaxon, c("TaxonName", "TaxonCode"), "\\(")
df.risk.areas.taxa <- df.risk.areas.taxa %>% separate(TaxonCode, c("TaxonCode"), "\\)")
# Remove records of taxa for which score is not available
taxa_2_rm <- unique(df.risk.areas.taxa$TaxonCode[!(df.risk.areas.taxa$TaxonCode %in% scores.data$Taxon_Code)])
df.risk.areas.taxa <- df.risk.areas.taxa[!(df.risk.areas.taxa$TaxonCode %in% taxa_2_rm), ]
list_taxa_code <- unique(df.risk.areas.taxa$TaxonCode)

# List of VME risk areas
vme_risk_areas_code <- intersect(df.risk.areas.taxa$VMECode, df.risk.areas$VMECode)
# Build dataset
abundance_df <- df.risk.areas[df.risk.areas$VMECode %in% vme_risk_areas_code, ]
abundance_df[ ,c("Number VME taxa", "VME-indicator units", "DepthMid")] <- list(NULL)

# Get Abundance data
for (row in 1:nrow(abundance_df)) {
  dat_area <- df.risk.areas.taxa[df.risk.areas.taxa$VMECode == abundance_df[row, "VMECode"], c("TaxonCode", "VMESpecimenWeight")]
  taxa_cur <- unique(dat_area$TaxonCode)
  for (taxon in taxa_cur) {
    abundance_df[row, taxon] <- sum(as.numeric(dat_area[dat_area$TaxonCode == taxon, "VMESpecimenWeight"]))
  }
}

# Get Jenks breaks for abundance data
abundance_vals <- as.numeric(df.risk.areas.taxa$VMESpecimenWeight)
abundance_breaks <- getJenksBreaks(abundance_vals, n_abundance_categories+1)

# Get Abundance scores
df_abundance_score <- data.frame(abundance_df)
for (row in 1:nrow(abundance_df)) {
  for (taxon in list_taxa_code) {
    abundance_taxon <- df_abundance_score[row, taxon]
    if (! is.na(abundance_taxon)) {
      for (abundance_break_idx in 1:(length(abundance_breaks)-1)) {
        abundance_break <- abundance_breaks[abundance_break_idx]
        if (abundance_taxon >= abundance_break) {
          df_abundance_score[row, taxon] <- abundance_break_idx
        } else {
          break
        }
      }
    }
  }
}

# Modulate abundance score by indicator score
df_indicator_abundance_score <- data.frame(df_abundance_score)
for (taxon in list_taxa_code) {
  taxon_score <- scores.data[scores.data$Taxon_Code == taxon, "score"]
  df_indicator_abundance_score[ , taxon] <- df_abundance_score[ , taxon] * taxon_score
}

# Aggregate scores to compute VME index
if (vme_index_agg == "max") {
  df_indicator_abundance_score$vme_index <- do.call(pmax, c(df_indicator_abundance_score[list_taxa_code], list(na.rm=TRUE)))
}
# TODO: add other methods
print(head(df_indicator_abundance_score[c("VMECode", "vme_index")]))

# Get Jenks breaks for vme indexes
vme_index_vals <- as.numeric(df_indicator_abundance_score$vme_index)
vme_index_breaks <- getJenksBreaks(vme_index_vals, n_index_categories+1)

# Get final VME index
df_vme_index <- data.frame(df_indicator_abundance_score)
df_vme_index <- df_vme_index[ , !(colnames(df_vme_index) %in% list_taxa_code)]
for (row in 1:nrow(df_vme_index)) {
  vme_index_area <- df_vme_index[row, "vme_index"]
  for (vme_index_break_idx in 1:(length(vme_index_breaks)-1)) {
    vme_index_break <- vme_index_breaks[vme_index_break_idx]
    if (vme_index_area >= vme_index_break) {
      df_vme_index[row, "vme_index"] <- vme_index_break_idx
    } else {
      break
    }
  }
}

# Get CCAMLR coordinates
coords <- data.frame(x = as.numeric(df_vme_index$LongitudeMid), 
                     y= as.numeric(df_vme_index$LatitudeMid))
coordinates(coords) <- c("x","y")
projection(coords) <- "+proj=longlat +datum=WGS84"
coords <- SOproj(coords)
SOmap_auto(coords, input_points = FALSE, input_lines = FALSE)
plot(coords, add = TRUE, pch = 19, col = 3, cex=1)

# Create mask around one of the records, with some buffer
mask <- buffer(coords[20, ], 50000)
#mask <- buffer(coords[1, ], 50000)
#mask <- buffer(coords[2, ], 50000)
#mask <- buffer(coords[75, ], 50000)
#mask <- buffer(coords[48, ], 50000)
SOmap_auto(coords, input_points = FALSE, input_lines = FALSE)
plot(coords, add = TRUE, pch = 19, col = 3, cex=1)
plot(mask, add=TRUE)

# Crop raster around this mask
raster.crop <- crop(SmallBathy, mask)
# Aesthetic
blue.col <- colorRampPalette(c("darkblue", "lightblue"))
blue.br <- seq(from=-11000, to=0, by=1000)
# Plot
plot(raster.crop, col=blue.col(22), breaks=blue.br, legend=FALSE)
plot(coords, add = TRUE, pch = 19, col = 3, cex=1)

# Create Risk areas circle shapes
df_risk_areas <- data.frame(long = as.numeric(df_vme_index$LongitudeMid),
                            lat = as.numeric(df_vme_index$LatitudeMid))
sf_risk_areas <- st_as_sf(df_risk_areas, coords = c("long", "lat"), crs="+proj=longlat +datum=WGS84")
# Buffer circles by 1 nmile radius
sf_circles <- st_buffer(sf_risk_areas, dist = 1852)
sf_circles <- SOproj(sf_circles)
sf_circles$vme <- as.numeric(df_vme_index$vme_index)
yellow.col <- colorRampPalette(c("yellow", "orange"))
yellow.br <- seq(from=1, to=max(n_index_categories), by=1)
plot(sf_circles, add=TRUE, col=yellow.col(n_index_categories), legend=TRUE)
# Final resolution
if (resolution.raster < res(raster.crop)[1]) {
  resample_fact <- res(raster.crop)[1] / resolution.raster
  raster.crop.res <- disaggregate(raster.crop, fact=resample_fact, method='bilinear')
} else if (resolution.raster > res(raster.crop)[1]) {
  resample_fact <- resolution.raster / res(raster.crop)[1]
  raster.crop.res <- aggregate(raster.crop, fact=resample_fact, method='bilinear')
} else {
  raster.crop.res <- raster.crop
}
crs(raster.crop.res) <- "+proj=longlat +datum=WGS84"
# Rasterize
raster_circles <- rasterize(sf_circles, raster.crop.res, field="vme", fun="max") #, fun=mean)
unique(values(raster_circles))
unique(sf_circles$vme)
crs(raster_circles) <- "+proj=longlat +datum=WGS84"
# TODO: MOve between sf_circles and raster_circles
yellow.col <- colorRampPalette(c("yellow", "orange"))
yellow.br <- seq(from=1, to=max(n_index_categories)+1, by=1)
par(mar=c(2, 2, 2, 2))
plot(raster.crop.res, col=blue.col(22), breaks=blue.br, legend=FALSE)
plot(rasterToPolygons(raster_circles), add=TRUE, col=yellow.col(n_index_categories), breaks=yellow.br, legend=FALSE, border="black")
legend("topright", legend=unique(raster_circles), col=yellow.col(n_index_categories), pch=15, border="black", title = "VME index")
plot(sf_circles, add=TRUE, col=yellow.col(6), legend=TRUE)

par(mar=c(2, 2, 2, 2))
plot(raster.crop.res, col=blue.col(22), breaks=blue.br, legend=FALSE)
plot(rasterToPolygons(raster_circles), add=TRUE, col=yellow.col(n_index_categories), breaks=yellow.br, legend=FALSE, border="black")



