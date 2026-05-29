# ----------------------------------------------------------
# Function: topX.calc
# ----------------------------------------------------------
# Identifies the top (or bottom) X proportion of raster values
#
# Arguments:
#   ra  : input raster
#   x   : proportion threshold (e.g. 0.1 = top 10%)
#   low : if TRUE selects lowest X%, otherwise highest
#
# Returns:
#   Boolean raster (1 = hotspot, NA = masked cells)
# ----------------------------------------------------------
topX.calc <- function(ra, x = 0.1, low = FALSE) {
  # total number of values
  n_cells <- sum(!is.na(values(ra)))
  # Select exactly top X% cells
  topX <- selectHighest(ra, ceiling(x * n_cells), low = low)
  # Convert NA (outside selection) to 0
  topX[is.na(topX)] <- 0
  # Restore NA where original raster had NA
  topX[is.na(ra)] <- NA
  
  return(topX)
}

# ----------------------------------------------------------
# Function: calc_percentile_raster
# ----------------------------------------------------------
# Converts raster values into percentiles (0–1)
#
# Each cell value represents its rank relative to all cells
#
# Returns:
#   Raster with values between 0 and 1
# ----------------------------------------------------------
calc_percentile_raster <- function(ra) {
  vals <- values(ra)
  # Rank values (ignoring NA)
  ranks <- rank(vals, na.last = "keep")
  # Convert ranks to percentiles
  percentiles <- ranks / length(na.omit(vals))
  # Assign back to raster
  out <- setValues(ra, percentiles)
  
  return(out)
}
