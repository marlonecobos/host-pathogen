if (!require("scales")) {
  install.packages("scales")
  Sys.sleep(1)
  library(scales)
}

if (!require("ks")) {
  install.packages("ks")
  Sys.sleep(1)
  library(ks)
}


#' Kernel density plots to explore data used in niche_signal tests
#' 
#' @param niche_signal_list list of results from niche_signal.
#' @param variables (character) name of variable(s) to used in plots. 
#' The default, NULL, uses the variable used in univariate tests or the first 
#' two variables used in PERMANOVA tests.
#' @param col_all a color lines representing all values. Default = "black".
#' @param col_positive a color lines representing positive values. Default = "red".
#' @param alpha new alpha level, representing transparency. Default = 0.5.
#' This value is adjusted (x 0.75) for all data and (x 1.25) for positives.
#' @param nlevels number of levels to be plotted in 2D kernel contour lines.
#' @param lwd_all (numeric) line width for all values. See options in 
#' \code{\link[graphics]{par}}. Default = 1.5.
#' @param lwd_positive (numeric) line width for positive values, default = 0.7.  
#' @param ... other arguments for the general plotting area.
#' 
#' @rdname plot_niche_signal_compsarison
#' @usage 
#' plot_data_kernel(niche_signal_list, variables = NULL, 
#'                  col_all = "black", col_positive = "red", 
#'                  alpha = 0.5, nlevels = 20, lwd_all = 1.5, 
#'                  lwd_positive = 0.7, ...)


plot_data_kernel <- function(niche_signal_list, variables = NULL, 
                             col_all = "black", col_positive = "red", 
                             alpha = 0.5, nlevels = 20, lwd_all = 1.5, 
                             lwd_positive = 0.7, ...) {
  
  # method using to produce niche_signal_list
  method <- niche_signal_list$summary$method
  
  if (!is.null(method)) {
    cond <- niche_signal_list$summary$condition 
    
    # plotting results from general analysis
    if (method == "univariate") {
      # variables
      if (is.null(variables)) {
        vars <- niche_signal_list$summary$variables[1] 
      } else {
        vars <- variables[1]
      }
      
      condi <- niche_signal_list$data[, cond] == 1
      data <- list(all = niche_signal_list$data[, vars], 
                   positive = niche_signal_list$data[condi, vars])
      clas <- "uni"
    } else {
      data <- niche_signal_list$modified_data
      clas <- "mult"
    }
  } else {
    # plotting results from univariate analysis
    if (!is.null(niche_signal_list$analysis_results)) {
      # variables
      if (is.null(variables)) {
        vars <- niche_signal_list$summary$variable 
      } else {
        vars <- variables[1]
      }
      
      condi <- niche_signal_list$data[, cond] == 1
      data <- list(all = niche_signal_list$data[, vars], 
                   positive = niche_signal_list$data[condi, vars])
      clas <- "uni"
    } else {
      # plotting results from permanova results
      if (!is.null(niche_signal_list$permanova_results)) {
        data <- niche_signal_list$modified_data
        clas <- "mult"
      } 
    }
  }
  
  # plots
  # colors
  alppos <- ifelse((alpha * 1.25) > 1, alpha, alpha * 1.25)
  call <- scales::alpha(col_all, alpha = alpha * 0.75)
  cpos <- scales::alpha(col_positive, alpha = alppos)
  
  
  if (clas == "uni") {
    # density
    dall <- density(data$all)
    dpos <- density(data$positive, )
    
    # limits
    lims <- apply(cbind(c(dall$x, dpos$x), c(dall$y, dpos$y)), 2, range)
    colnames(lims) <- c(vars, "Density")
    
    # plot
    plot(lims, type = "n", ...)
    lines(dall, col = call, lwd = lwd_all)
    lines(dpos, col = cpos, lwd = lwd_positive)
    
  } else {
    # variables
    if (is.null(variables)) {
      vars <- niche_signal_list$summary$variables[1:2] 
    } else {
      vars <- variables[1:2]
    }
    
    # kernels
    dall <- ks::kde(data[data[, cond] == 0, vars])
    dpos <- ks::kde(data[data[, cond] == 1, vars])
    
    # limits
    lims <- apply(do.call(cbind, dall$eval.points), 2, range)
    colnames(lims) <- vars
    
    # plot
    plot(lims, type = "n", ...)
    contour(dall$eval.points[[1]], dall$eval.points[[2]], z = dall$estimate, 
            col = call, drawlabels = FALSE, nlevels = nlevels, 
            lwd = lwd_all, add = TRUE)
    contour(dpos$eval.points[[1]], dpos$eval.points[[2]], z = dpos$estimate, 
            col = cpos, drawlabels = FALSE, nlevels = nlevels, 
            lwd = lwd_positive, add = TRUE)
  }
  
}
