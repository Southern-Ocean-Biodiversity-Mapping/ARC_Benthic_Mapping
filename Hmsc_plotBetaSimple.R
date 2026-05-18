# Requires the fields package for image.plot
# Install if needed: install.packages("fields")
# Load the package: library(fields)

plotBetaSimple <- function(
    hM,                    # HMSC model object
    post,                  # Posterior distribution object
    param = "Support",     # Type of parameter to plot: "Mean", "Support", or "Sign"
    supportLevel = 0.9,    # Threshold for posterior support
    speciesNames = TRUE,   # Whether to show species names
    covariateNames = TRUE, # Whether to show covariate names
    colors = colorRampPalette(c("blue", "white", "red")), # Color scheme
    colorLevels = NULL,    # Number of color levels
    cex = c(0.7, 0.7),    # Text size for covariate and species names
    mar = c(8, 15, 2, 2), # Margins
    spSort = "none",      # Sorting option: "none", "alphabetical", or "reverse"
    main = NULL,          # Title for the plot
    ynames = TRUE         # labels on y-axis
) {
  # Check if fields package is loaded
  if (!requireNamespace("fields", quietly = TRUE)) {
    stop("Please install and load the 'fields' package to use this function")
  }
  
  # Adjust color levels based on parameter
  if (is.null(colorLevels)) {
    colorLevels <- switch(param,
                          "Sign" = 3,
                          "Support" = 3,
                          "Mean" = 200
    )
  }
  
  # Adjust colors based on parameter
  if (param == "Support") {
    colors <- colorRampPalette(c("blue", "white", "red"))(3)
  } else {
    colors <- colors(colorLevels)
  }
  
  # Get species and covariate names
  spNames <- if(speciesNames) hM$spNames else paste0("S", 1:hM$ns)
  covNames <- if(covariateNames) hM$covNames else paste0("C", 1:hM$nc)
  
  # Create order based on sorting option
  spOrder <- switch(spSort,
                    "alphabetical" = order(spNames),
                    "reverse" = order(spNames, decreasing = TRUE),
                    "none" = 1:length(spNames)
  )
  
  spNames <- spNames[spOrder]
  
  mbeta <- post$mean
  betaP <- post$support
  
  toPlot <- switch(param,
                   "Sign" = sign(mbeta),
                   "Mean" = mbeta,
                   "Support" = 2 * betaP - 1
  )
  
  toPlot <- toPlot * ((betaP > supportLevel) + (betaP < (1 - supportLevel)) > 0)
  betaMat <- matrix(toPlot, nrow = hM$nc, ncol = ncol(hM$Y))
  
  # Apply the sorting to the matrix columns
  betaMat <- betaMat[, spOrder]
  
  # Ensure betaMat is treated as a matrix even if it has only one column
  if (is.null(dim(betaMat))) {
    betaMat <- matrix(betaMat, nrow = hM$nc, ncol = 1)
  }
  
  rownames(betaMat) <- covNames
  colnames(betaMat) <- spNames
  
  par(mar = mar)
  
  X <- t(betaMat)
  zlim <- if(all(is.na(X)) || sum(abs(X)) == 0) {
    c(-1, 1)
  } else {
    c(-max(abs(range(X))), max(abs(range(X))))
  }
  
  # Create sequence for positioning
  x_pos <- seq(0.1, 0.9, length.out = ncol(X))
  y_pos <- seq(0.1, 0.9, length.out = nrow(X))
  
  # Customize legend labels based on parameter
  legend.args <- if(param == "Support") {
    list(
      labels = c("Negative", "None", "Positive"),
      at = c(-0.67, 0, 0.67),  # Centered positions for the three categories
      cex.axis = cex[1]
    )
  } else {
    list(cex.axis = cex[1])
  }
  
  # Plot without axes
  fields::image.plot(
    x = x_pos,
    y = y_pos,
    z = t(X),
    xlab = "",
    ylab = "",
    col = colors,
    zlim = zlim,
    axes = FALSE,
    main = main,
    axis.args = legend.args
  )
  
  # Add custom axes with only our labels
  axis(1, at = x_pos, labels = covNames, las = 2, cex.axis = cex[1])
  if(ynames) axis(2, at = y_pos, labels = spNames, las = 2, cex.axis = cex[2])
  
  # Add box around the plot
  box()
}
