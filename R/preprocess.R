#‘ Smooth raw intensity data by loess for one gel
#'
#' Initial smoothing using the obtained raw intensities for all sample lanes in one gel
#'
#' @param dat_raw Raw data obtained from \code{\link{read_IPdata}}. First Column is from 0 to 1, evenly
#' split grids; Other columns record intensity values obtained from all the tested sample lanes.
#' @param span The parameter tha controls the degree of smoothing. See \code{\link[stats]{loess}}
#'
#' @return A list comprising two elements, \code{res} for smoothed intensities, one column
#' per sample lane, and \code{sigma} for the standard deviation
#' @export
smooth_gel <- function(dat_raw,span){
  curr_dat0 <- dat_raw
  curr_dat  <- dat_raw[,-1] # remove the first linear increasing reference.
  res   <- dat_raw
  sigma <- dat_raw
  for (lane in 1:ncol(curr_dat)){
    fit <- msir::loess.sd(curr_dat0[,1],curr_dat[,lane],
                          nsigma = 1,degree=0, span=span)
    x   <- curr_dat0[,1]
    y   <- fit$y
    res[,lane+1]   <- y
    sigma[,lane+1] <- fit$sd
  }
  list(res=res,sigma=sigma)
}


#’ Binning smoothed autoradio intensity data
#'
#' Force multiple gels to share identical bins
#'
#' @param dat_smoothed Smoothed autoradiointensity data obtained from \code{\link{smooth_gel}}.
#' First column contains a set of evenly split grid points from 0 to 1; other columns
#' contain smoothed (e.g., by loess) intensity values.
#' @param break_points The break points defining the bins (shorter than the no. of rows )
#'
#' @return A data frame with binned data, calculated from all intensity data
#' belonging to the same bin.
#'
#' @export
bin_gel <- function(dat_smoothed,break_points){
  group <- cut(dat_smoothed$rf,break_points,include.lowest = TRUE) # <-- binning according to break points
  form  <- stats::as.formula(paste0("cbind(rf,",
                             paste(colnames(dat_smoothed)[-1],collapse=", "),
                             ")~group"))
  stats::aggregate(form,dat_smoothed,FUN = mean)
}


#' remove the right end of the gel image
#'
#' Generally performed before pattern recognition, because light molecules cannot
#' be distinguished and hence clustered at the righ end in the experiment.
#'
#' @param dat_binned Binned autoradiointensity data
#' @param percent Default is 4.4. The percentage of bins at the right end to be cropped.
#' @param bp Break points that defined the binning. See \code{\link{bin_gel}}.
#'
#' @return A list of two elements:
#' \itemize{
#' \item \code{dat_binned_cropped} Binned data after cropping
#' \item \code{break_points_cropped} Break points after cropping
#' }
#' @export
crop_end <- function(dat_binned,percent=4.4,bp){
  n_g <- ncol(dat_binned)-2
  n_bin <- nrow(dat_binned)
  n_bin_to_rm <- floor(n_bin*percent/100)
  list(dat_binned_cropped=dat_binned[1:(n_bin-n_bin_to_rm),],
       break_points_cropped=bp[-1][1:(n_bin-n_bin_to_rm)])
}



