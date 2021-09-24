library(readxl)    
read_excel_allsheets <- function(filename) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

curate_ccamlr_registry <- function(data) {
  # VME Summary
  vme.summary <- data[[2]]
  colnames(vme.summary) <- vme.summary[2, ]
  vme.summary <- vme.summary[-c(1, 2), ]
  vme.summary <- vme.summary[, -c(1, 2, 5)]
  # VME risk areas
  vme.risk.areas <- data[[4]]
  colnames(vme.risk.areas) <- vme.risk.areas[2, ]
  vme.risk.areas <- vme.risk.areas[-c(1, 2), ]
  vme.risk.areas <- vme.risk.areas[, -c(1, 2, 5, 11)]
  # VME risk areas taxa
  vme.risk.areas.taxa <- data[[5]]
  colnames(vme.risk.areas.taxa) <- vme.risk.areas.taxa[2, ]
  vme.risk.areas.taxa <- vme.risk.areas.taxa[-c(1, 2), ]
  # VME Fine Scale rectangles
  vme.fsr <- data[[6]]
  colnames(vme.fsr) <- vme.fsr[2, ]
  vme.fsr <- vme.fsr[-c(1, 2), ]
  # Store in list
  list(vme = vme.summary,
       vme.risk.areas = vme.risk.areas,
       vme.risk.areas.taxa = vme.risk.areas.taxa,
       vme.fsr = vme.fsr)
}

exit <- function() { invokeRestart("abort") }

get_jenks_breaks <- function(vals, n_categories) {
  vals <- vals[!is.na(vals)]
  return(getJenksBreaks(vals, n_categories+1))
}

apply_jenks_breaks <- function(vals, list_jenks_breaks, n_category) {
  categories <- c()
  for (val in vals) {
    if (is.na(val)) {
      category <- NA
    } else {
      category <- 1
      if (val == tail(list_jenks_breaks, n=1)) {
        category <- n_category
      } else {
        for (break_idx in 1:(length(list_jenks_breaks)-1)) {
          break_ <- list_jenks_breaks[break_idx]
          if (val > break_) {
            category <- break_idx
          } else {
            break
          }
        }
      }
    }
    categories <- c(categories, category)
  }
  return(categories)
}

get_points_in_asd <- function(df_, asd_object, asd_name) {
  asd_interest <- asd_object[asd_object$GAR_Short_Label==asd_name, ]
  asd_interest_spPoly = as(asd_interest, "SpatialPolygons")
  
  pts <- df_
  coordinates(pts) <- c("proj_coord_x", "proj_coord_y")
  projection(pts) <- crs(asd_interest)
  
  pts_in = pts[!is.na(over(pts, asd_interest_spPoly)), ]
  df_interest = as.data.frame(pts_in)
  
  cat("\nSelecting points in ASD #", asd_name, " ...\n")
  cat("Number of points: ", nrow(df_interest), "\n")
  return(df_interest)
}