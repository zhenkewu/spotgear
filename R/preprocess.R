#‘ Smooth raw intensity data by loess for one gel
#'
#' Initial smoothing using the obtained raw intensities for all sample lanes in one gel
#'
#' @param dat_raw Raw data obtained from \code{\link{read_IPdata}}. First column represents
#' a equal-spaced grid over 0 and 1; Other columns record intensity values obtained from all
#' the tested sample lanes.
#' @param span The parameter that controls the degree of smoothing. See \code{\link[stats]{loess}}
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
#' belonging to the same bin. The first two columns are bin labels and equal-spaced grids, respectively.
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
#' @param bp Break points that defined the binning in \code{\link{bin_gel}}.
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


#' Create data for Bayesian 2D dewarping
#'
#' Construct data frame that store all detected bands (NB: the result can be trimmed to exclude specific rows.)
#'
#' @param score_mat_list A list of matrices. Each matrix has B rows and N_g columns. NA means non-peaks, a number,
#' usually 3, represents peaks.
#' @param analysis_id choose the gels in this analysis; user's choice.
#' @param rf_grid A vector of length B (=670) equal-spaced grid from 0 to 1.
#' @param N_g_vec A named numeric vector showing the number of lanes (incl. Lane 1s)
#'   for each gel with the name being the original ID before being subset by analysis_id.
#'
#' @return A data frame comprised of 5 elements:
#' \itemize{
#' \item \code{gel_ID}
#' \item \code{lane_ID}
#' \item \code{band_ID}
#' \item \code{Y} Location on Gel for the peaks
#' \item \code{peak_bin_index} The bin index for the peaks
#' \item \code{lane_ID_stacked} cumulative counts.
#' }
#'
#'
#' @export
build_data_for_dewarping <- function(score_mat_list,analysis_id,rf_grid,N_g_vec){
  peak_ind    <- apply(do.call(cbind,score_mat_list[analysis_id]),2,
                       function(v) which(!is.na(v))) # <- find the bin number for peaks.
  Y          <- rf_grid[unlist(peak_ind)]   # <-- normalized location for a peak.

  J_vec    <- unlist(lapply(peak_ind,length)) # no. of peaks per lane. Length is sum(n_g_vec). w/ lane 1.
  n_g_vec  <- N_g_vec[analysis_id] # no. of lanes per gel, w/ lane 1

  lane_ID_stacked   <- rep(1:length(J_vec),J_vec) # <----- temporary ID for serum lane, repeated for peaks in the same lane.
  gel_ID     <- as.factor(rep(as.numeric(names(n_g_vec)),times=n_g_vec)[lane_ID_stacked]) # <-- ID for gel.
  band_ID    <- unlist(sapply(1:length(J_vec),function(s) {
    if (J_vec[s]==0){
      print(paste0("==[spotgear] Found a lane without detected bands."))
      return(NULL)}
    1:(J_vec[s])
    })) # <--- ID for autoantibody bands.
  lane_ID    <- lane_ID_stacked-c(0,cumsum(n_g_vec)[-length(analysis_id)])[as.numeric(gel_ID)]

  data.frame(gel_ID    = gel_ID,
             lane_ID   = lane_ID,
             band_ID   = band_ID,
             Y         = Y,
             peak_bin_index = unlist(peak_ind),
             lane_ID_stacked = lane_ID_stacked) # for animation, not necesary for dewarping.
}
















