################################################################################
#                                 INIT
################################################################################
# Set working directory
setwd("C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/vme_index")

# Source method variables and load packages
source("config.R")
# Source utils functions
source("utils.R")

################################################################################
#                                 LOAD DATA
################################################################################
# Read CCAMLR Registry
ccamlr_registry <- read_excel_allsheets(path_ccamlr_registry)
# Curate data
data_ccamlr <- curate_ccamlr_registry(ccamlr_registry)
# Read VME Taxa Scores
vulnerability_scores_taxa <- read.csv(path_taxa_scores)
head(vulnerability_scores_taxa)

################################################################################
#                   COMPUTE VULNERABILITY SCORE OF EACH TAXA
################################################################################
cat("Computing vulnerability score of each taxa using",
    toupper(agg_vulnerability_score),
    "aggregation method...")
if (agg_vulnerability_score == "mean") {
  result <- apply(vulnerability_scores_taxa[, 5:11], 1,
                  function(x) mean(x, na.rm=TRUE))
} else if (agg_vulnerability_score == "quadratic_mean") {
  quadratic_mean <- function(x) {
    x_ <- x[!is.na(x)]
    return(sqrt(sum(x_^2)/length(x_)))
  }
  result <- apply(vulnerability_scores_taxa[, 5:11], 1,
                  function(x) quadratic_mean(x))
} else {
  cat("\tERROR: Unknown aggregation method:", agg_vulnerability_score,
      "\n\t\tPlease choose among:", c("mean", "quadratic_mean"))
  exit()
}
vulnerability_scores_taxa$score <- result
vulnerability_scores_taxa <- vulnerability_scores_taxa[order(vulnerability_scores_taxa$score,
                                                             decreasing = TRUE),]
print(vulnerability_scores_taxa[ , c("Taxon", "score")])

################################################################################
#                           CURRATNG RISK AREAS DATA
################################################################################
cat("Curating VME risk areas data...")
# Get VME Risk Areas data
df_risk_areas <- data_ccamlr[["vme.risk.areas"]]
df_risk_areas_taxa <- data_ccamlr[["vme.risk.areas.taxa"]]
df_risk_areas_taxa <- df_risk_areas_taxa %>% separate(VMEIndicatorTaxon, c("TaxonName", "TaxonCode"), "\\(")
df_risk_areas_taxa <- df_risk_areas_taxa %>% separate(TaxonCode, c("TaxonCode"), "\\)")
# Remove records of taxa for which score is not available
taxa_2_rm <- unique(df_risk_areas_taxa$TaxonCode[!(df_risk_areas_taxa$TaxonCode %in% vulnerability_scores_taxa$Taxon_Code)])
cat("\tWARNING: Removing the following taxa records:", taxa_2_rm,
    " ... because vulnerability score is not available for these taxa.")
df_risk_areas_taxa <- df_risk_areas_taxa[!(df_risk_areas_taxa$TaxonCode %in% taxa_2_rm), ]
# Remove areas where abundance data is not available
df_risk_areas_taxa$VMESpecimenWeight <- as.numeric(df_risk_areas_taxa$VMESpecimenWeight)
cat("\tWARNING: Removing", sum(is.na(df_risk_areas_taxa$VMESpecimenWeight)),
    "taxa records because VMESpecimenWeight is not available for them.")
df_risk_areas_taxa <- df_risk_areas_taxa[!is.na(df_risk_areas_taxa$VMESpecimenWeight), ]
# List taxa included in the analysis
list_taxa_code <- unique(df_risk_areas_taxa$TaxonCode)
cat("\tINFO: Available data include", length(list_taxa_code), "VME indicator taxa.")
# List of VME risk areas
list_vme_risk_areas_code <- intersect(df_risk_areas_taxa$VMECode, df_risk_areas$VMECode)
cat("\tINFO: Available data include", length(list_vme_risk_areas_code), "VME Risk Areas.")
df_abundance <- df_risk_areas[df_risk_areas$VMECode %in% list_vme_risk_areas_code, ]
df_abundance[ ,c("Number VME taxa", "VME-indicator units", "DepthMid")] <- list(NULL)

################################################################################
#                           CURRATNG RISK AREAS DATA
################################################################################
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



