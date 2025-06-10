#-----------------------------------------------------------------------------80
# scatter_w_histogram
#-----------------------------------------------------------------------------80
#' Two-dimensional scatter plot with histograms along the x and y axes.
#'
#' This function visualizes two-dimensional data using a scatter plot, 
#'   with histograms displayed along the x and y axes.
#'
#' @param dataframe2D A data frame with two numeric columns.
#' @param logscale Logical: If TRUE, histogram frequencies are displayed on a
#'   log1p scale.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(x = rnorm(100), y = 0.5 * rnorm(100) + rnorm(100))
#' scatter_w_histogram(dataframe2D = df, logscale = TRUE)
#'
scatter_w_histogram <- function(dataframe2D, logscale = FALSE){
  #-----------------------------------------------
  # 0. Preparation
  #-----------------------------------------------
  x <- dataframe2D[[1]]
  y <- dataframe2D[[2]]
  layout(matrix(c(1, 0, 2, 3), nrow = 2, byrow = TRUE),
         widths = c(4, 1.5), heights = c(1.5, 4))
  #-----------------------------------------------
  # 1. Histogram for x-axis
  #-----------------------------------------------
  par(mar = c(0, 4, 2, 0.5))
  h_x <- hist(x, plot = FALSE, breaks = 20)
  if(logscale){
    barplot(log1p(h_x$counts), axes = FALSE, space = 0, col = "gray80",
            border = "white")
  }else{
    barplot(h_x$counts, axes = FALSE, space = 0, col = "gray80",
            border = "white")
  }
  axis(2, las = 1, cex.axis = 0.7)
  #-----------------------------------------------
  # 2. Scatter plot
  #-----------------------------------------------
  par(mar = c(4, 4, 0.5, 0.5))
  plot(x, y, pch = 19, col = "steelblue", xlab = "x", ylab = "y")
  #-----------------------------------------------
  # 3. Histogram for y-axis
  #-----------------------------------------------
  par(mar = c(4, 0.5, 0.5, 2))
  h_y <- hist(y, plot = FALSE, breaks = 20)
  if(logscale){
    barplot(log1p(h_y$counts), horiz = TRUE, axes = FALSE, space = 0,
            col = "gray80", border = "white")
  }else{
    barplot(h_y$counts, horiz = TRUE, axes = FALSE, space = 0,
            col = "gray80", border = "white")
  }
  axis(1, cex.axis = 0.7)
}
#-----------------------------------------------------------------------------80
# scatter_w_gradient
#-----------------------------------------------------------------------------80
#' Show a 2D scatter plot showing gene expression levels with a color gradient.
#'
#' This function visualizes the expression of a specified gene across cells
#'   on a two-dimensional plot. Each point is colored based on the normalized
#'   expression level of the gene, using a smooth gradient to highlight
#'   differences.
#'
#' @param dataframe2D A data frame containing two numeric columns representing
#'   coordinates (e.g., PCA embeddings) for each cell.
#' @param mat A gene-by-cell expression matrix (e.g., raw counts or
#'   normalized values), where columns correspond to the rows of dataframe2D.
#' @param gene The name of the gene to visualize, which must be present in the
#'   rownames of mat.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(x = rnorm(50), y = 0.5 * rnorm(50) + rnorm(50))
#' mat <- matrix(rpois(10 * 50, lambda = 5), nrow = 10)
#' rownames(mat) <- paste0("GENE_", 1:10)
#' colnames(mat) <- paste0("Cell_", 1:50)
#' scatter_w_gradient(dataframe2D = df, mat = mat, gene = "GENE_5")
#'
scatter_w_gradient <- function(dataframe2D, mat, gene){
  #-----------------------------------------------
  # Extract gene expression vector.
  #-----------------------------------------------
  expr <- mat[gene, ]

  #-----------------------------------------------
  # Normalize expression values to [0,1] for color mapping.
  #-----------------------------------------------
  expr_norm <- (expr - min(expr)) / (max(expr) - min(expr))

  #-----------------------------------------------
  # Define color palette (approximate viridis "plasma").
  #-----------------------------------------------
  colors <- c("#0D0887", "#6A00A8", "#B12A90", "#E16462", "#FCA636")
  color_palette <- colorRampPalette(colors)(100)
  #-----------------------------------------------
  # Map normalized expression to colors.
  #-----------------------------------------------
  point_colors <- color_palette[as.numeric(cut(expr_norm, breaks = 100))]

  #-----------------------------------------------
  # Expand x-axis limits more to reserve space for colorbar and labels.
  #-----------------------------------------------
  xlim_expanded <- range(dataframe2D[,1])
  x_span <- diff(xlim_expanded)
  # Increased margin on right.
  xlim_expanded[2] <- xlim_expanded[2] + x_span * 0.22  
  #-----------------------------------------------
  # Plot scatter points with colors.
  #-----------------------------------------------
  plot(dataframe2D, col = point_colors, pch = 19, xlim = xlim_expanded,
       main = paste("Expression of", gene))
  #-----------------------------------------------
  # Enable overlay drawing for colorbar.
  #-----------------------------------------------
  par(new = TRUE)

  #-----------------------------------------------
  # Get current plot coordinates.
  #-----------------------------------------------
  usr <- par("usr")  # c(xmin, xmax, ymin, ymax)

  #-----------------------------------------------
  # Define colorbar rectangle coordinates on right margin.
  #-----------------------------------------------
  xleft <- usr[2] - x_span * 0.18
  xright <- usr[2] - x_span * 0.13
  #-----------------------------------------------
  # Add small vertical margin so colorbar won't get clipped.
  #-----------------------------------------------
  ybottom <- usr[3] + (usr[4] - usr[3]) * 0.02  
  ytop <- usr[4] - (usr[4] - usr[3]) * 0.02
  #-----------------------------------------------
  # Use original color order: dark = low, bright = high.
  #-----------------------------------------------
  colorbar_colors <- color_palette
  n_colors <- length(colorbar_colors)
  rect_height <- (ytop - ybottom) / n_colors
  #-----------------------------------------------
  # Draw colorbar as stacked rectangles.
  #-----------------------------------------------
  for(i in seq_len(n_colors)){
    rect(xleft, ybottom + (i - 1) * rect_height, xright,
         ybottom + i * rect_height, col = colorbar_colors[i], border = NA)
  }
  #-----------------------------------------------
  # Add legend labels next to colorbar with adjusted positions to avoid cutoff.
  #-----------------------------------------------
  text_x_pos <- xright + x_span * 0.03  # shift further right
  text(y = c(ybottom, (ybottom + ytop)/2, ytop), x = text_x_pos,
       labels = c("Low", "Medium", "High"), adj = 0, cex = 0.8)
  #-----------------------------------------------
  # Reset overlay drawing flag.
  #-----------------------------------------------
  par(new = FALSE)
}

