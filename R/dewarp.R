if(getRversion() >= "2.15.1") utils::globalVariables(c("gel_ID","lane_ID"))

#' Piecewise linear warping function for one index in the reference
#'
#' Do piecewise stretch and compression by matching knots_query to knots_ref, using
#' the latter as anchors. The entire stretch and compression is done by \code{\link{correct_batch_effect_by_ref}}.
#'
#' @param i The i-th bin in the reference
#' @param knots_ref The knots for the reference
#' @param knots_query The knots for the lane that is to be aligned to the reference
#'
#' @return The index in the query lane that is to be aligned to the i-th bin in the reference lane
#' @export
pwl <- function(i,knots_ref,knots_query){ # must include endpoints; knots are in query.
  if(length(knots_ref)!=length(knots_query)){stop("==unequal matching points==")}
  k <- as.numeric(cut(i,breaks = knots_ref,include.lowest =TRUE))
  w <- (i-knots_ref[k])/(knots_ref[k+1]-knots_ref[k])
  w*knots_query[k+1]+(1-w)*knots_query[k]
}


#' Piecewise-linear warping to correct batch effects
#'
#' Correct the error that the same sample with known molecules can show up at
#' slightly different locations when autoradiographed on multiple gels.
#'
#' @param dat_binned_query The binned data to be compressed or stretched towards the reference
#' @param dat_binned_ref The binned data to be aligned to. It is usually arbitrarily chosen
#' or the one least visually warped (could be very subjective). We sugguest in one particular
#' analysis, fix this data to be one specific gel.
#' @param peak_ind_mat_query A matrix of values, with the first-column entries that equal \code{score_thres} considered
#' as peaks in the query lane
#' @param peak_ind_mat_ref Same as \code{peak_ind_mat_query} but for the reference gel
#' to be aligned to.
#' @param bp The break points that binned the smoothed data into bins.
#' @param score_thres The score threshold above or equal to which we call peak candidate bins.
#' @return A list of two elements, the first is the vector of break points at reference
#' peak locations (large to small; NB:check this), the second is a vector of the same length
#' as the legnth of reference lane, with each element \code{warped_ind[i]} being the index in the query lane
#' that corresponds for reference index \code{i}.
#'
#' @export
#'
correct_batch_effect_by_ref <- function(dat_binned_query,dat_binned_ref,
                                        peak_ind_mat_query,peak_ind_mat_ref,
                                        bp,score_thres){
  ref_index    <- which(peak_ind_mat_ref[,1]==score_thres) # <- first lane has to be reference.
  knots_ref   <- c(1,ref_index,nrow(dat_binned_ref))
  ref_loc      <- rev(sort(bp[ref_index]))

  knots_query <- c(1,which(peak_ind_mat_query[,1]==score_thres),nrow(dat_binned_query))
  warped_ind  <- sapply(1:nrow(dat_binned_query),pwl,knots_ref,knots_query)
  make_list(ref_loc,warped_ind)
}



#' Create the binned data after piecewise linear warping
#'
#' @param dat_binned_list A list of data each being the binned data
#' @param score_mat_peak_list A list of score matrix with peak locations filled
#'  with numeric values
#' @param query_id The gel number that is to be aligned to the reference gel
#' @param ref_id The number of the reference gel
#' @param bp The break points that binned the smoothed data into bins.
#' @param score_thres The score threshold above or equal to which we call peak candidate bins.
#' @return We can name the results \code{dat_binned_batch_aligned_list}
#' @export
create_batch_aligned_data <- function(dat_binned_list,
                                      score_mat_peak_list,
                                      query_id,
                                      ref_id,
                                      bp,
                                      score_thres){
  if (query_id == ref_id){
    res <- dat_binned_list[[ref_id]]
  } else{
    out_batch_align <- correct_batch_effect_by_ref(dat_binned_list[[query_id]],
                                                   dat_binned_list[[ref_id]],
                                                   score_mat_peak_list[[query_id]],
                                                   score_mat_peak_list[[ref_id]],
                                                   bp, score_thres)
    res         <- dat_binned_list[[query_id]][out_batch_align$warped_ind,]
    res$rf      <- dat_binned_list[[ref_id]]$rf
  }
  res
}

#' Piecewise linear warping after aligning the landmarks
#'
#' Because at the right side of the gel, we have few peaks, data contain few
#' information about the dewarping function. We downsample the gel, so that we only
#' need to consider the dewarping results on the J=50 landmarks' dewarping results.
#' NB: assumed that all gels have 19 non-reference lanes. Need to change otherwise.
#'
#' @param g_id Gel id
#' @param l_id Lane id
#' @param normalized_rf The normalized_rf that is between 0 and 1. Note this subtle change compared
#' after cropping the right end of the gel.
#' @param dat_without_lane1 The data frame without peaks obtained from the reference lanes. Because
#' when doing warping, we know these reference molecules do not exactly correspond to the
#' molecules that would be observed in the serum samples.
#' @param data_lcm_all A binary matrix of the number of rows being the number of serum samples
#' and the number of columns being J+2 (J for intermediate landmarks and 2 end points)
#' @param zero_to_one_v_landmark The landmark location on the zero to one scale.
#'
#' @return A list of two elements:
#' \itemize{
#' \item \code{ref_index} Among \code{normalized_rf}, the index of the one that is
#' closest to each of the peak-matched landmarks. Because we want to use these index as anchors to do
#' piecewise linear warping.
#' \item \code{warped_ind} Warping mapping: bin index for post-dewarping data |-> bin index for the pre-dewarping data.
#' }
#' @export
pwl_after_dewarping    <- function(g_id,l_id,normalized_rf,dat_without_lane1,data_lcm_all,
                                   zero_to_one_v_landmark){
  n_total_bins  <- length(normalized_rf)
  ref           <- data_lcm_all[19*(g_id-1)+l_id-1,]
  ref_index     <- sapply(zero_to_one_v_landmark[ref > 0],function(v) which.min(abs(normalized_rf-v))) #<-- find among normalized_rf, which ones are close to each of landmarks.
  knots_base    <- c(1,ref_index,n_total_bins)

  knots_query   <- c(1,subset(dat_without_lane1,gel_ID==g_id & lane_ID==l_id)$peak_bin_index,n_total_bins)
  warped_ind    <- sapply(1:n_total_bins,pwl,knots_base,knots_query)
  make_list(ref_index,warped_ind)
}
