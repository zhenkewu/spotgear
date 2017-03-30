#' Get the two-sided end points
#'
#' Centered at the i-th bin, find the bin hald_width away to the left and to the right;
#' Truncate at the 1st or the biggest bin when hitting the boundary. Essential
#' for evaluating if a bin is likely to be the location of a peak (aka band)
#'
#' @param i i-th bin
#' @param half_width The number of bins to the left or right of the current bin (specified by
#' parameter \code{i})
#' @param n The total number of bins
#'
#' @return A vector of length two, representing the left and right interval end points
#' @export
#'
edges<- function(i,half_width,n){ # n is the entire length of bins.
  res <- c(NA,NA)
  if (i < half_width+1){
    res <- c(1,i+half_width)
  } else if (i+half_width>n){
    res <- c(i-half_width,n)
  } else{
    res <- c(i-half_width, i+half_width)
  }
  res
}

#' Calculate relative scores for peak calling
#'
#' Sign calculation for each bin in the binned intensity data. It is used
#' to define contiguous regions in which we will find a peak.
#'
#' @param y The binned intensity data
#' @param half_width The half window size used in \code{\link{edges}}
#' @param elevation The minimum intensity difference a peak needs to exceed the minimum
#' intensity within the window defined by \code{\link{edges}}
#'
#' @return A vector of scores, one for each intensity value in a bin
#' @export
#'
sign_score <-
  function(y,half_width,elevation){
    n   <- length(y)
    res <- rep(NA,n)
    for (i in 1:n){
      tmp <- edges(i,half_width,n)
      left <- tmp[1]
      right <- tmp[2]
      neib_vec <- y[left:right]
      res[i]   <- sign(y[i]-y[left])+sign(y[i]-y[right])+
        sign(y[i]-min(neib_vec)-elevation)        # <--- significant differences.
    }
    res
  }



#' Compute local relative scores:
#'
#' Compute using sign_score for each lane to find peak candidate regions
#'
#' @param dat_binned Binned autoradiointensity data
#' @param half_width See \code{\link{edges}}
#' @param lanes The index of lanes to be analyzed
#' @param elev_vec  See the argument \code{elevation} in \code{\link{sign_score}}.Default 0.014.
#'
#' @return a matrix of relative scores (nrow=number of bins, ncol=number of lanes analyzed)
#' @export
compute_local_relative_score <- function(dat_binned,half_width,lanes,
                                         elev_vec=rep(0.014,length(lanes))){
  n_g <- length(lanes)
  score_mat <- matrix(NA,nrow=nrow(dat_binned),ncol=n_g)
  for (l in lanes){
    y <- dat_binned[,l+2] # 2: 1st column for interval group, 2nd column for rf.
    score_mat[,l] <- sign_score(y,half_width,elev_vec[l])
  }
  score_mat
}


#' Get the bins for bands:
#'
#' Given each contiguous peak candidate regions, find the bin that we call peak (aka band)
#'
#'
#' @param score_mat The matrix of scores calcualted by \code{\link{compute_local_relative_score}}
#' @param score_thres The threshold greater than or equal to which we call the peak-candidate regions (
#' must be contiguous)
#' @param dat_binned Binned autoradiointensity data
#' @param gap_considered_identical_peak Default is 5 bins. To merge nearby contiguous peak-candidate
#' regions if they are less than or equal to \code{gap_considered_indentical_peak} away.
#'
#' @return A matrix with \code{NA} and \code{score_thres}. Those entries equaling \code{score_thres} are bins for peaks.
#' @export

compute_score_mat_peak <- function(score_mat,score_thres,dat_binned,
                                   gap_considered_identical_peak=5){
  
  # score_mat <- score_mat_list[[s]]
  # dat_binned <- dat_binned_list[[s]]
  # gap_considered_identical_peak <- 5
  
  n_g <- ncol(score_mat)
  score_mat_peak <- matrix(NA,nrow=nrow(dat_binned),ncol=n_g)
  for (lane in 1:n_g){
    curr_score <- score_mat[,lane]
    res        <- rle(curr_score)
    humps      <- which(res$values >= score_thres & res$lengths>=1)

    runs_lengths_cumsum <- cumsum(res$lengths)
    ends0 <- runs_lengths_cumsum[humps]

    newindex <- ifelse(humps>1, humps-1, 0)
    starts0 <- runs_lengths_cumsum[newindex] + 1
    if (0 %in% newindex) starts0 = c(1,starts0)

    #starts0 <- c(2,5,6,7,8,10,21)
    #ends0  <-  c(3,5,6,7,9,14,30)
    # put very nearby peaks together: (when computing peaks from the scores, we join very close contiguous bin segments).
    if (length(starts0)>1){
      small_gap_index <- which(starts0[-1]-ends0[-length(ends0)] <=gap_considered_identical_peak)
      if (length(small_gap_index)==0){
        starts <- starts0
        ends   <- ends0
      } else{
        ends <- ends0
        starts <- starts0
        for (j in rev(seq_along(small_gap_index))){
          ends[small_gap_index[j]] <- ends[small_gap_index[j]+1]
        }
        ends <- ends[-(small_gap_index+1)]
        starts <- starts[-(small_gap_index+1)]
      }
    } else{
      starts <- starts0
      ends   <- ends0
    }

    # print(starts)
    # print(ends)
    duration <- ends-starts+1
    starts <- starts[duration>=3]
    ends <- ends[duration>=3]
    for (r in seq_along(starts)){
      local_max_ind <- c(starts[r]:ends[r])[which.max(dat_binned[starts[r]:ends[r],lane+2])]
      score_mat_peak[local_max_ind,lane] <- score_thres
    }
  }
  score_mat_peak
}
