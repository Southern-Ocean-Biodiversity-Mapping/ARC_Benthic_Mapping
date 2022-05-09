# Load config file
readRenviron("config.R")
path.ccamlr.registry = Sys.getenv("path.ccamlr.registry")
proj.ref = Sys.getenv("proj_def")
path.raster = Sys.getenv("path.raster")
resolution.raster = as.numeric(Sys.getenv("resolution.raster"))

# Load functions
source("utils.R")

# Packages
library(raster)
library(fasterize)
library(RColorBrewer)
library(SOmap)

# Import raster
raster.500 <- raster(path.raster)
# Resample it to resolution.raster
raster.res <- raster(crs=crs(raster.500), ext=extent(raster.500), res=resolution.raster)
# Raster Ross Sea
raster.ross.sea <- raster(crs=crs(raster.500), xmn=-13500, xmx=533500, ymn=-2037000, ymx=-1537000, res=res(raster.res))
SOmap_auto(raster.ross.sea)

# Read CCAMLR registry
ccamlr.registry <- read_excel_allsheets(path.ccamlr.registry)
# Curate data
ccamlr.data <- curate_ccamlr_registry(ccamlr.registry)

# Get VME data
df.vme <- ccamlr.data[["vme"]]
## TODO: Include Longitude Start / End etc.
# Get CCAMLR coordinates
coords <- data.frame(x = as.numeric(df.vme$LongitudeMid), 
                     y= as.numeric(df.vme$LatitudeMid))
coordinates(coords) <- c("x","y")
projection(coords) <- proj.ref
coords <- SOproj(coords)
# Get CCAMLR data
data <- data.frame(vme = as.numeric(df.vme$NumberVMETaxa))
# Assign CCAMLR values
cells <- cellFromXY(raster.res, coords)
raster.res[cells] <- data$vme
# Plot
SOmap_auto(raster.res, input_points = FALSE, input_lines = FALSE)
plot(raster.res, add=TRUE, col=brewer.pal(n = 3, name = "YlOrRd"))
## TODO: Some Zoom-in with finer raster
## To plot without auto scaling
#plot(base_map)
## To plot without raster
#plot(spatial.dat, add = TRUE, pch = 19, col = 3, cex=1)
SOmap_auto(coords, input_points = FALSE, input_lines = FALSE)
plot(coords, add = TRUE, pch = 19, col = 3, cex=1)



# Get VME Risk Areas
df.risk.areas <- ccamlr.data[["vme.risk.areas"]]
df.risk.areas.taxa <- ccamlr.data[["vme.risk.areas.taxa"]]
## TODO: Include Longitude Start / End etc.
## TODO: Combine Number VME taxa and VME-indicator units
# Get CCAMLR coordinates
coords <- data.frame(x = as.numeric(df.risk.areas$LongitudeMid), 
                     y= as.numeric(df.risk.areas$LatitudeMid))
coordinates(coords) <- c("x","y")
projection(coords) <- proj.ref
coords <- SOproj(coords)
# Get CCAMLR data
data <- data.frame(vme = as.numeric(df.risk.areas$`Number VME taxa`))
# Assign CCAMLR values
cells <- cellFromXY(raster.res, coords)
raster.res[cells] <- data$vme
# Plot
SOmap_auto(raster.res, input_points = FALSE, input_lines = FALSE)
plot(raster.res, add=TRUE, col=brewer.pal(n = 3, name = "YlOrRd"))
# Get CCAMLR data
data <- data.frame(vme = as.numeric(df.risk.areas$`VME-indicator units`))
# Assign CCAMLR values
cells <- cellFromXY(raster.res, coords)
raster.res[cells] <- data$vme
# Plot
SOmap_auto(raster.res, input_points = FALSE, input_lines = FALSE)
plot(raster.res, add=TRUE, col=brewer.pal(n = 3, name = "YlOrRd"))
## TODO: Some Zoom-in with finer raster
## To plot without auto scaling
#plot(base_map)
## To plot without raster
#plot(spatial.dat, add = TRUE, pch = 19, col = 3, cex=1)
SOmap_auto(coords, input_points = FALSE, input_lines = FALSE)
plot(coords, add = TRUE, pch = 19, col = 3, cex=1)




dat_st <- data.frame(long = as.numeric(df.risk.areas$LongitudeMid), 
                     lat = as.numeric(df.risk.areas$LatitudeMid))
#vme = as.numeric(df.risk.areas$`Number VME taxa`))
dat_sf <- st_as_sf(dat_st, coords = c("long", "lat"), crs = proj.ref)
# Buffer circles by 1 nmile radius
dat_circles <- st_buffer(dat_sf, dist = 1852)
dat_circles_proj <- SOproj(dat_circles)
dat_circles_proj$vme <- as.numeric(df.risk.areas$`Number VME taxa`)
#SOmap_auto(raster.ross.sea)
#plot(dat_circles_proj, add=TRUE)
#raster_circles <- rasterize(dat_circles_proj, raster.ross.sea, field="vme", background = NA_real_)
#SOmap_auto(raster.ross.sea)
#plot(raster_circles, add=TRUE)














library(CCAMLRGIS)
#Load ASDs
ASDs=load_ASDs()
#Subsample ASDs to only keep Subarea 88.1
S881=ASDs[ASDs$GAR_Short_Label=='881',]
#Crop bathymetry to match the extent of S881
B881=raster::crop(SmallBathy,S881)
#Optional: get the maximum depth in that area to constrain the color scale
minD=raster::minValue(S881)
#Set the figure margins as c(bottom, left, top, right)
par(mai=c(0.2,0.4,0.2,0.55))
#Plot the bathymetry
plot(S881,breaks=Depth_cuts,col=Depth_cols,legend=F,axes=F,box=F)
#Add color scale
add_Cscale(height=80,fontsize=0.7,offset=300,width=15,lwd=0.5,minVal=minD,maxVal=-1)
#Add coastline (for Subarea 48.6 only)
plot(Coast[Coast$ID=='88.1',],col='grey',lwd=0.01,add=T)
#Add reference grid
add_RefGrid(bb=bbox(B881),ResLat=5,ResLon=10,fontsize=0.75,lwd=0.75,offset = 100000)
#Add Subarea 88.1 boundaries
plot(S881,add=T,lwd=1,border='red')
#Add a -2000m contour
raster::contour(B881,levels=-2000,add=T,lwd=0.5,labcex=0.3)
#Add single label at the centre of the polygon (see ?Labels)
text(Labels$x[Labels$t=='88.1'],Labels$y[Labels$t=='88.1'],labels='88.1',col='red',cex=1.5)
data_vme <- data.frame(Lon = as.numeric(df.risk.areas$LongitudeMid), 
                       Lat = as.numeric(df.risk.areas$LatitudeMid),
                       VME_taxa = as.numeric(df.risk.areas$`Number VME taxa`))
#coordinates(data_vme) <- c("Lon","Lat")
#data_vme <- SOproj(data_vme)
MyGrid=create_PolyGrids(data_vme,dlon=1,dlat=2)
plot(MyGrid,col=MyGrid$VME_taxa,main='Example 1',cex.main=0.75,lwd=0.1)
box()


raster.s881 <- raster::crop(raster.500, S881)
raster.s881.50km <- raster::aggregate(raster.s881, fact=100)
#raster.s881.50km <- raster(crs=crs(raster.s881.50km), ext=extent(raster.s881.50km), res=res(raster.s881.50km))
#dat_sf$vme <- as.numeric(df.risk.areas$`Number VME taxa`)
#a <- rasterize(dat_sf, raster.s881.50km, field="vme", fun=sum)
#plot(a)

# Get VME FSR
df.fsr <- ccamlr.data[["vme.fsr"]]
## TODO: get coordinates right
## Plot




#projection(spatial.dat) <- SOproj(proj.ref)

#SOmap_auto(spatial.dat, bathy = "space")
#SOplot(spatial.dat, col = "green", pch=20)

#legend("topright",
#       legend = unique(S$Survey),
#       col = 2:cmpt,
#       pch=3,
#       cex = 0.8)

#so.map <- SOgg(SOmap_auto(spatial.dat))
#so.map$bathy <- NULL
#myggplot <- plot(so.map)
#Y.na <- SXY
#Y.na[Y.na == 0] <- NA
#for (g in c("Sponges", "Hydroids")) {
#  print(myggplot + geom_point(data=as.data.frame(Y.na), aes_string("Longitude", "Latitude", color=g)) + scale_colour_distiller(palette = "YlOrRd", direction=1, limits = c(1e-7, 10)) + labs(title = g, subtitle = "Percentage cover"))
#}
