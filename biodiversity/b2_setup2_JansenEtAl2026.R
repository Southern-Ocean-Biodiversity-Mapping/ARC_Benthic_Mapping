############################################################
# BASIC CHECKS
############################################################
if (!dir.exists(pred.model.dir)){stop("Predictive map directory not found: ", pred.model.dir)}
if (!dir.exists(cont.dir))      {stop("Continuous output directory not found: ", cont.dir)}
if (!dir.exists(pct.dir))       {stop("Percentile output directory not found: ", pct.dir)}
if (!dir.exists(hot.dir))       {stop("Hotspot output directory not found: ", hot.dir)}
if (!dir.exists(domain.dir))    {stop("Domain output directory not found: ", domain.dir)}

############################################################
# LOAD SPECIES-LEVEL PREDICTIVE MAPS
############################################################
# Load all PA raster files for the selected model
pa.files <- list.files(  pred.model.dir, pattern = "^PA_.*\\.tif$", full.names = TRUE)
if (length(pa.files) == 0) {stop("No PA_*.tif files found in: ", pred.model.dir)}

# Load only the layers needed in this document:
# 1 = mean, 2 = median, 3 = standard error
hmsc.maps.pa.env.mean   <- rast(pa.files, lyrs = 1)
hmsc.maps.pa.env.median <- rast(pa.files, lyrs = 2)
hmsc.maps.pa.env.se     <- rast(pa.files, lyrs = 3)

############################################################
# LOAD TOTAL ABUNDANCE
############################################################
# Load current continuous abundance outputs
abundance.mean       <- rast(file.path(cont.dir, "abundance_mean.tif"))
abundance.median     <- rast(file.path(cont.dir, "abundance_median.tif"))
abundance.median.100 <- rast(file.path(cont.dir, "abundance_median_100.tif"))
abundance.se         <- rast(file.path(cont.dir, "abundance_se.tif"))

# Rebuild abundance object used in the figure code
hmsc.maps.totalabundance.env <- c(
  abundance.median,
  abundance.median.100,
  abundance.se
)
names(hmsc.maps.totalabundance.env) <- c("median", "median.100", "se")

############################################################
# LOAD DERIVED PREDICTIVE MAPS
############################################################

# Richness outputs
richness.median <- rast(file.path(cont.dir, "richness_median.tif"))
richness.boot.median <- rast(file.path(cont.dir, "richness_bootstrap_median.tif"))
richness.boot.sd <- rast(file.path(cont.dir, "richness_bootstrap_sd.tif"))

# # Richness / abundance percentiles
# perc.richness.median  <- rast(file.path(pct.dir, "richness_percentiles.tif"))
# perc.abundance.median <- rast(file.path(pct.dir, "abundance_percentiles.tif"))
# perc.rich_per_abund   <- rast(file.path(pct.dir, "richness_per_abundance_percentiles.tif"))
# perc.abund_per_rich   <- rast(file.path(pct.dir, "abundance_per_richness_percentiles.tif"))

biodiversity.hot <- rast(file.path(hot.dir, "biodiversity_hotspots.tif"))

top.percentiles <- rast(file.path(hot.dir, "all_hotspots_file.tif"))


# Domain-level percentiles
perc.richness.domains  <- rast(file.path(domain.dir, "domain_richness_percentiles.tif"))
perc.abundance.domains <- rast(file.path(domain.dir, "domain_abundance_percentiles.tif"))

# Domain-level hotspot classes reconstructed from current outputs
domain.rich.hot.class <- rast(file.path(domain.dir, "domain_richness_hotspot_class.tif"))
domain.abund.hot.class <- rast(file.path(domain.dir, "domain_abundance_hotspot_class.tif"))
domain.biodiv.hot.class <- rast(file.path(domain.dir, "domain_biodiversity_hotspot_class.tif"))

# Domain-level hotspot thresholds used for polygon reconstruction
domain.rich.hot.top10 <- rast(file.path(domain.dir, "domain_richness_hotspot_top10.tif"))
domain.abund.hot.top10 <- rast(file.path(domain.dir, "domain_abundance_hotspot_top10.tif"))
domain.biodiv.hot.top10 <- rast(file.path(domain.dir, "domain_biodiversity_hotspot_top10.tif"))

top.percentiles.domains <- c(
  domain.rich.hot.class,
  domain.abund.hot.class,
  domain.biodiv.hot.class
)
names(top.percentiles.domains) <- c(
  paste0("richness domain ", c(1, 3, 4, 7, 8, 9)),
  paste0("abundance domain ", c(1, 3, 4, 7, 8, 9)),
  paste0("biodiversity domain ", c(1, 3, 4, 7, 8, 9))
)

# Beta-diversity products are not currently part of the updated pipeline.
# Create placeholder rasters so downstream code can still run unless those
# sections are explicitly used.
median.betadiv <- rast(richness.median)
median.betadiv[] <- NA
names(median.betadiv) <- "median.betadiv"

perc.betadiv.median <- rast(richness.median)
perc.betadiv.median[] <- NA
names(perc.betadiv.median) <- "perc.betadiv.median"


############################################################
# MATCH NA STRUCTURE ACROSS DERIVED RASTERS
############################################################
ra.na.sel <- which(is.na(richness.median[]))

median.betadiv[ra.na.sel]    <- NA
richness.boot.median[ra.na.sel] <- NA
richness.boot.sd[ra.na.sel]     <- NA
perc.betadiv.median[ra.na.sel] <- NA

############################################################
# CREATE HOTSPOT POLYGONS
############################################################

# Circumpolar polygons are written directly by the hotspot pipeline
top.richness.polygons <- vect(file.path(poly.dir, "richness_top10_hotspots.gpkg"))
top.abundance.polygons <- vect(file.path(poly.dir, "abundance_top10_hotspots.gpkg"))
top.biodiv.polygons <- vect(file.path(poly.dir, "biodiversity_top10_hotspots.gpkg"))

# Domain polygons are reconstructed from the domain top-10 hotspot rasters
top.dom.rich.polygons <- as.polygons(domain.rich.hot.top10[[1]], dissolve = TRUE)
for (i in 2:nlyr(domain.rich.hot.top10)) {
  top.dom.rich.polygons <- c(
    top.dom.rich.polygons,
    as.polygons(domain.rich.hot.top10[[i]], dissolve = TRUE)
  )
}

top.dom.abund.polygons <- as.polygons(domain.abund.hot.top10[[1]], dissolve = TRUE)
for (i in 2:nlyr(domain.abund.hot.top10)) {
  top.dom.abund.polygons <- c(
    top.dom.abund.polygons,
    as.polygons(domain.abund.hot.top10[[i]], dissolve = TRUE)
  )
}

top.dom.biodiv.polygons <- as.polygons(domain.biodiv.hot.top10[[1]], dissolve = TRUE)
for (i in 2:nlyr(domain.biodiv.hot.top10)) {
  top.dom.biodiv.polygons <- c(
    top.dom.biodiv.polygons,
    as.polygons(domain.biodiv.hot.top10[[i]], dissolve = TRUE)
  )
}

############################################################
# DATA FOR VIOLIN PLOTS
############################################################
## first match the NAs for cells between the environmental and biodiversity layers
r.stack.subset <- mask(r.stack.subset, richness.median)

# #### values benthic biodiveristy hotspots vs all other areas
# ## extract environmental conditions
# r.stack.insideAES    <- mask(r.stack.subset,        top.biodiv.polygons)
# r.stack.outsideAES.1 <- mask(r.stack.subset[[1:2]], top.biodiv.polygons, inverse=TRUE)
# r.stack.outsideAES.2 <- mask(r.stack.subset[[3:4]], top.biodiv.polygons, inverse=TRUE)
# r.stack.outsideAES.3 <- mask(r.stack.subset[[5:6]], top.biodiv.polygons, inverse=TRUE)
# r.stack.outsideAES.4 <- mask(r.stack.subset[[7:8]], top.biodiv.polygons, inverse=TRUE)
# r.stack.outsideAES.5 <- mask(r.stack.subset[[9:11]],top.biodiv.polygons, inverse=TRUE)
# r.stack.outsideAES <- c(r.stack.outsideAES.1,r.stack.outsideAES.2,r.stack.outsideAES.3,r.stack.outsideAES.4,r.stack.outsideAES.5)
# rm(r.stack.outsideAES.1,r.stack.outsideAES.2,r.stack.outsideAES.3,r.stack.outsideAES.4,r.stack.outsideAES.5)
# # ## extract values
# vals.in <- list()
# vals.out <- list()
# for(i in 1:nlyr(r.stack.insideAES)){
#   vals.in[[i]] <- values(r.stack.insideAES[[i]], na.rm=TRUE)
#   vals.out[[i]] <- values(r.stack.outsideAES[[i]], na.rm=TRUE)
# }
# names(vals.in) <- names(vals.out) <- names(r.stack.subset)
# ## save output
# save(vals.in, vals.out,file=file.path(pred.dir, "BiodiversityHotspots_EnvironmentalSetting_CafeNPPandFAM.Rdata"))
load(file.path(pred.dir, "BiodiversityHotspots_EnvironmentalSetting_CafeNPPandFAM.Rdata"))

## translate flux to normalised values, displayed on log scale
all_vals <- c(exp(vals.out$log.flux.mean.cafe), exp(vals.in$log.flux.mean.cafe))
global_median <- median(all_vals, na.rm = TRUE)

vals.out$rel_median_flux.mean.cafe <- exp(vals.out$log.flux.mean.cafe) / global_median
vals.in$rel_median_flux.mean.cafe  <- exp(vals.in$log.flux.mean.cafe)  / global_median

vals.out$log_rel_median_flux.mean.cafe <- log10(vals.out$rel_median_flux.mean.cafe)
vals.in$log_rel_median_flux.mean.cafe  <- log10(vals.in$rel_median_flux.mean.cafe)

## translate sed to normalised values
all_vals_sed <- c(vals.out$sed.mean.cafe, vals.in$sed.mean.cafe)
global_median_sed <- median(all_vals_sed, na.rm = TRUE)

vals.out$rel_median_sed.mean.cafe <- vals.out$sed.mean.cafe / global_median_sed
vals.in$rel_median_sed.mean.cafe  <- vals.in$sed.mean.cafe  / global_median_sed


#############
#### values for predator hotspots and overlap with benthic biodiveristy hotspots
# r.stack.predatoroverlap <- mask(r.stack.subset, top.combined.only.polygons)
# ## which values are inside the polygons
# vals.o.in <- list()
# for(i in 1:nlyr(r.stack.predatoroverlap)){
#   vals.o.in[[i]] <- values(r.stack.predatoroverlap[[i]], na.rm=TRUE)
# }
# r.stack.predators <- mask(r.stack.subset, predator.shelf.top.polygons)
# ## extract values that are inside the polygons
# vals.p.in <- list()
# for(i in 1:nlyr(r.stack.predators)){
#   vals.p.in[[i]] <- values(r.stack.predators[[i]], na.rm=TRUE)
# }
# save(vals.o.in, vals.p.in, file=file.path(pred.dir, "BiodiversityHotspots_EnvironmentalSetting_CafeNPPandFAM_PredatorsAndOverlap.Rdata"))
load(file.path(pred.dir, "BiodiversityHotspots_EnvironmentalSetting_CafeNPPandFAM_PredatorsAndOverlap.Rdata"))

############
#### values for hotspots in each planning domain
# # Select variables
# domain_ids <- c(1, 3, 4, 7, 8, 9)
# domain_indices <- 1:6 # Corresponds to perc.richness.domains and top.dom.biodiv.polygons
# 
# # sel <- c(1, 6, 9, 10)
# # var_names <- c("depth", "curr", "temp", "flux")
# sel <- 1:11
# var_names <- c("depth","slope","tpi","distance2canyons","log.flux","sed.mean","npp_mean","curr_mean","curr_residual","salinity","temp")
# 
# # Create masked subsets
# r.stack.subset.dom <- lapply(domain_indices, function(i) {
#   mask(r.stack.subset[[sel]], perc.richness.domains[[i]])
# })
# names(r.stack.subset.dom) <- paste0("dom", domain_ids)
# 
# # Domain-specific hotspots and outside areas
# r.stack.domtop <- list()
# r.stack.domtopout <- list()
# for (i in seq_along(domain_ids)) {
#   gc()
#   print(i)
#   dom <- paste0("dom", domain_ids[i])
#   subset_dom <- r.stack.subset.dom[[i]]
#   poly <- top.dom.biodiv.polygons[[i]]
#   ## inside
#   r.stack.domtop[[dom]] <- mask(subset_dom, poly)
#   ## outside
#   out_a <- mask(subset_dom[[1:2]], poly, inverse = TRUE)
#   out_b <- mask(subset_dom[[3:4]], poly, inverse = TRUE)
#   out_c <- mask(subset_dom[[5:6]], poly, inverse = TRUE)
#   out_d <- mask(subset_dom[[7:8]], poly, inverse = TRUE)
#   out_e <- mask(subset_dom[[9:11]], poly, inverse = TRUE)
#   r.stack.domtopout[[dom]] <- c(out_a, out_b,out_c, out_d, out_e)
# }
# 
# # Circumpolar hotspots and outside areas
# r.stack.circtop <- list()
# r.stack.circtopout <- list()
# for (i in seq_along(domain_ids)) {
#   gc()
#   print(i)
#   dom <- paste0("dom", domain_ids[i])
#   subset_dom <- r.stack.subset.dom[[i]]
#   ## inside
#   r.stack.circtop[[dom]] <- mask(subset_dom, top.biodiv.polygons)
#   ## outside
#   out_a <- mask(subset_dom[[1:2]], top.biodiv.polygons, inverse = TRUE)
#   out_b <- mask(subset_dom[[3:4]], top.biodiv.polygons, inverse = TRUE)
#   out_c <- mask(subset_dom[[5:6]], top.biodiv.polygons, inverse = TRUE)
#   out_d <- mask(subset_dom[[7:8]], top.biodiv.polygons, inverse = TRUE)
#   out_e <- mask(subset_dom[[9:11]],top.biodiv.polygons, inverse = TRUE)
#   r.stack.circtopout[[dom]] <- c(out_a, out_b,out_c, out_d, out_e)
# }
# 
# ## Setup functions to extract values
# extract_values <- function(stack) {
#   lapply(1:11, function(i) values(stack[[i]], na.rm = TRUE))
# }
# 
# ## Extract values
# vals.domtop <- lapply(names(r.stack.domtop), function(dom) {extract_values(r.stack.domtop[[dom]])})
# names(vals.domtop) <- names(r.stack.domtop)
# vals.domtopout <- lapply(names(r.stack.domtopout), function(dom) {extract_values(r.stack.domtopout[[dom]])})
# names(vals.domtopout) <- names(r.stack.domtopout)
# 
# vals.circtop <- lapply(names(r.stack.circtop), function(dom) {extract_values(r.stack.circtop[[dom]])})
# names(vals.circtop) <- names(r.stack.circtop)
# vals.circtopout <- lapply(names(r.stack.circtopout), function(dom) {extract_values(r.stack.circtopout[[dom]])})
# names(vals.circtopout) <- names(r.stack.circtopout)
# save(vals.domtop, vals.circtop, vals.domtopout, vals.circtopout, file=file.path(pred.dir, "BiodiversityHotspots_EnvironmentalSetting_CafeNPPandFAM_Domains.Rdata"))
load(file.path(pred.dir, "BiodiversityHotspots_EnvironmentalSetting_CafeNPPandFAM_Domains.Rdata"))






















###########
# # Continuous surfaces
# rich.per.abund <- rast(file.path(cont.dir, "richness_per_abundance.tif"))
# abund.per.rich <- rast(file.path(cont.dir, "abundance_per_richness.tif"))
# 
# # Hotspot rasters
# richness.hot <- rast(file.path(hot.dir, "richness_hotspots.tif"))
# abundance.hot <- rast(file.path(hot.dir, "abundance_hotspots.tif"))
# richness.only.hot <- rast(file.path(hot.dir, "richness_only_hotspots.tif"))
# abundance.only.hot <- rast(file.path(hot.dir, "abundance_only_hotspots.tif"))
# rich.per.abund.hot <- rast(file.path(hot.dir, "richness_per_abundance_hotspots.tif"))
# abund.per.rich.hot <- rast(file.path(hot.dir, "abundance_per_richness_hotspots.tif"))
# 
# top.percentiles <- c(
#   richness.hot[[5]],
#   abundance.hot[[5]],
#   biodiversity.hot[[5]],
#   richness.only.hot,
#   abundance.only.hot,
#   rich.per.abund.hot[[5]],
#   abund.per.rich.hot[[5]]
# )
# names(top.percentiles) <- c(
#   "richness",
#   "abundance",
#   "biodiversity",
#   "richness_only",
#   "abundance_only",
#   "richness_per_abundance",
#   "abundance_per_richness"
# )
# ##### calculate hotspot threshold values
# ## biodiversity
# biodiv.rich.thresholds <- round(c(min(richness.median[top.percentiles$biodiversity[]==4]),
#                                   min(richness.median[top.percentiles$biodiversity[]==3]),
#                                   min(richness.median[top.percentiles$biodiversity[]==2]),
#                                   min(richness.median[top.percentiles$biodiversity[]==1]),
#                                   min(richness.median[top.percentiles$biodiversity[]==1])),0)
# biodiv.abund.thresholds <- round(c(min(abundance.median.100[top.percentiles$biodiversity[]==4]),
#                                    min(abundance.median.100[top.percentiles$biodiversity[]==3]),
#                                    min(abundance.median.100[top.percentiles$biodiversity[]==2]),
#                                    min(abundance.median.100[top.percentiles$biodiversity[]==1]),
#                                    min(abundance.median.100[top.percentiles$biodiversity[]==1])),0)
# levels(top.percentiles$biodiversity) <-
#   rev(c(paste0("> ", biodiv.abund.thresholds[1:4], " %-cover & > ",biodiv.rich.thresholds[1:4], " morphospecies"),
#         paste0("< ", biodiv.abund.thresholds[5],   " %-cover & < ",biodiv.rich.thresholds[5],   " morphospecies")))
# ## richness
# rich.rich.thresholds <- round(c(min(richness.median[top.percentiles$richness[]==4]),
#                                 min(richness.median[top.percentiles$richness[]==3]),
#                                 min(richness.median[top.percentiles$richness[]==2]),
#                                 min(richness.median[top.percentiles$richness[]==1]),
#                                 min(richness.median[top.percentiles$richness[]==1])),0)
# levels(top.percentiles$richness) <-
#   rev(c(paste0("> ",rich.rich.thresholds[1:4], " morphospecies"), paste0("< ",rich.rich.thresholds[5], " morphospecies")))
# ## abundance
# abund.abund.thresholds <- round(c(min(abundance.median.100[top.percentiles$abundance[]==4]),
#                                   min(abundance.median.100[top.percentiles$abundance[]==3]),
#                                   min(abundance.median.100[top.percentiles$abundance[]==2]),
#                                   min(abundance.median.100[top.percentiles$abundance[]==1]),
#                                   min(abundance.median.100[top.percentiles$abundance[]==1])),0)
# levels(top.percentiles$abundance) <-
#   rev(c(paste0("> ", abund.abund.thresholds[1:4], " %-cover"), paste0("< ", abund.abund.thresholds[5], " %-cover")))
# ## richness only
# rich.only.rich.thresholds <- round(c(min(richness.median[top.percentiles$richness_only[]==3]),
#                                      min(richness.median[top.percentiles$richness_only[]==2]),
#                                      min(richness.median[top.percentiles$richness_only[]==1]),
#                                      min(richness.median[top.percentiles$richness_only[]==1])),0)
# rich.only.abund.thresholds <- round(c(min(abundance.median.100[top.percentiles$richness_only[]==3]),
#                                       min(abundance.median.100[top.percentiles$richness_only[]==2]),
#                                       min(abundance.median.100[top.percentiles$richness_only[]==1]),
#                                       min(abundance.median.100[top.percentiles$richness_only[]==1])),0)
# levels(top.percentiles$richness_only) <-
#   rev(c(paste0("> ", rich.only.abund.thresholds[1:3], " %-cover & > ",rich.only.rich.thresholds[1:3], " morphospecies"),
#         paste0("< ", rich.only.abund.thresholds[4],   " %-cover & < ",rich.only.rich.thresholds[4],   " morphospecies")))
# ## abundance only
# abund.only.rich.thresholds <- round(c(min(richness.median[top.percentiles$abundance_only[]==2]),
#                                       min(richness.median[top.percentiles$abundance_only[]==1]),
#                                       min(richness.median[top.percentiles$abundance_only[]==1])),0)
# abund.only.abund.thresholds <- round(c(min(abundance.median.100[top.percentiles$abundance_only[]==2]),
#                                        min(abundance.median.100[top.percentiles$abundance_only[]==1]),
#                                        min(abundance.median.100[top.percentiles$abundance_only[]==1])),0)
# levels(top.percentiles$abundance_only) <-
#   rev(c(paste0("> ", abund.only.abund.thresholds[1:2], " %-cover & > ",abund.only.rich.thresholds[1:2], " morphospecies"),
#         paste0("< ", abund.only.abund.thresholds[3],   " %-cover & < ",abund.only.rich.thresholds[3],   " morphospecies")))
# ## richness per abundance
# rich.per.abund.thresholds <- round(c(min(rich.per.abund[top.percentiles$richness_per_abundance[]==4]),
#                                      min(rich.per.abund[top.percentiles$richness_per_abundance[]==3]),
#                                      min(rich.per.abund[top.percentiles$richness_per_abundance[]==2]),
#                                      min(rich.per.abund[top.percentiles$richness_per_abundance[]==1]),
#                                      min(rich.per.abund[top.percentiles$richness_per_abundance[]==1])),1)
# levels(top.percentiles$richness_per_abundance) <-
#   rev(c(paste0("> ",rich.per.abund.thresholds[1:4], " morphospecies for each 1% cover"),
#         paste0("< ",rich.per.abund.thresholds[5],   " morphospecies for each 1% cover")))
# ## abundance per richness
# abund.per.rich.thresholds <- round(c(min(abund.per.rich[top.percentiles$abundance_per_richness[]==4]),
#                                      min(abund.per.rich[top.percentiles$abundance_per_richness[]==3]),
#                                      min(abund.per.rich[top.percentiles$abundance_per_richness[]==2]),
#                                      min(abund.per.rich[top.percentiles$abundance_per_richness[]==1]),
#                                      min(abund.per.rich[top.percentiles$abundance_per_richness[]==1])),1)
# levels(top.percentiles$abundance_per_richness) <-
#   rev(c(paste0("< ",abund.per.rich.thresholds[1:4], " morphospecies for each 1% cover"),
#         paste0("> ",abund.per.rich.thresholds[5],   " morphospecies for each 1% cover")))
# ## output
