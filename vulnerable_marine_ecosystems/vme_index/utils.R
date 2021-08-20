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
