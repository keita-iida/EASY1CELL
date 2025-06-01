#-----------------------------------------------------------------------------80
# scatter_w_histogram
#-----------------------------------------------------------------------------80
#' Two-dimensional scatter plot with histograms along the x and y axes.
#'
#' This function visualizes two-dimensional data using a scatter plot, 
#' with histograms displayed along the x and y axes.
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

