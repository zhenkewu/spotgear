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
#' @param bugs_gel_id Starts from 1.
#' @param bugs_in_gel_id Starts from 1 to say 19.
#' @param N_per_gel Analyzed number of lanes per gel.
#' @param Z_map a vector of length equal to total number of peaks detected; each element 
#' represents the landmark a peak is aligned to
#' @param n_vec A vector of length equal to the total lane number; each element is the number of peaks within a lane.
#'
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
                                   zero_to_one_v_landmark,bugs_gel_id,bugs_in_gel_id,N_per_gel,
                                   Z_map,n_vec){
  g_id <- g
  l_id <- l
  bugs_gel_id <- s
  bugs_in_gel_id <- l-1
  dat_without_lane1 <- curr_dat
  
  n_total_bins  <- length(normalized_rf)
  G             <- length(N_per_gel)
  ref           <- data_lcm_all[c(0,cumsum(N_per_gel)[-G])[bugs_gel_id]+bugs_in_gel_id,]
  ref_index     <- sapply(zero_to_one_v_landmark[ref > 0],function(v) which.min(abs(normalized_rf-v))) #<-- find among normalized_rf, which ones are close to each of landmarks.
  knots_base    <- c(ifelse(min(ref_index)==1,NA,1),ref_index,
                     ifelse(max(ref_index)==n_total_bins,NA,n_total_bins))
  knots_base    <- knots_base[!is.na(knots_base)]
  
  knots_query   <- c(ifelse(min(ref_index)==1,NA,1),
                     subset(dat_without_lane1,gel_ID==g_id & lane_ID==l_id)$peak_bin_index,
                     ifelse(max(ref_index)==n_total_bins,NA,n_total_bins))
  knots_query    <- knots_query[!is.na(knots_query)]
  
  if (length(knots_base) < length(knots_query)){ # two observed peaks being aligned to identical landmark.
    bugs_stacked_lane_id <- c(0,cumsum(N_per_gel)[-G])[bugs_gel_id]+bugs_in_gel_id
    peak_to_landmark_vec_with_overlap <- Z_map[lookup[bugs_stacked_lane_id,1:n_vec[bugs_stacked_lane_id]]] # <- need n_vec, and Z_map.
    for (lmk in unique(peak_to_landmark_vec_with_overlap)){
       tmp_query     <- subset(dat_without_lane1,gel_ID==g_id & lane_ID==l_id)$peak_bin_index
       ind_to_unique <- which(peak_to_landmark_vec_with_overlap==lmk)
       tmp_query[ind_to_unique] <- floor(mean(tmp_query[ind_to_unique]))
       knots_query <- unique(tmp_query)    
    }
    
    knots_query   <- c(ifelse(min(ref_index)==1,NA,1),
                       knots_query,
                       ifelse(max(ref_index)==n_total_bins,NA,n_total_bins))
    knots_query    <- knots_query[!is.na(knots_query)]
    
  }
  warped_ind    <- sapply(1:n_total_bins,pwl,knots_base,knots_query)
  make_list(ref_index,warped_ind)
}


#' Fit Bayesian 2D image dewarping models
#'
#' @details This function prepares data, specifies hyperparameters in priors, 
#' initializes the posterior sampling chain, writes the model file (for JAGS syntax), 
#' and fits the model. Features:
#' 
#' If running JAGS on windows, please go to control panel to add the directory to
#' jags into ENVIRONMENTAL VARIABLE!
#'
#' @param data_dewarp Data frame with columns: gel_ID (true id), lane_ID (1 to say 20, 1 will be removed),
#' band_ID (counts within each lane), Y (t-scale location of peaks), peak_bin_index (index of the bin corresponding
#' to the peaks; requires common binning across images), lane_ID_stacked (cummulative lane number upon stacking
#' gels one-by-one; mainly for plotting).
#' @param mcmc_options A list of Markov chain Monte Carlo (MCMC) options.
#'
#' \itemize{
#' \item \code{debugstatus} Logical - whether to pause WinBUGS after it finishes
#' model fitting;
#' \item \code{n.chains} Number of MCMC chains;
#' \item \code{n.burnin} Number of burn-in samples;
#' \item \code{n.thin} To keep every other \code{n.thin} samples after burn-in period;
#' \item \code{result.folder} Path to folder storing the results;
#' \item \code{bugsmodel.dir} Path to WinBUGS model files;
#' }
#' @return BUGS fit results.
#' 
#' @import stats
#' 
#' @family model fitting functions 
#' 
#' @export
dewarp2d <- 
  function(data_dewarp,mcmc_options){
    gel_with_lane1 <- unique(subset(data_dewarp[,c("gel_ID","lane_ID")],lane_ID==1))[,1]
    if (sum(data_dewarp$lane_ID==1)>0){
      warning(paste0("==[spotgear] 'data_dewarp' contains Lane 1s for Gel ",
                     paste(gel_with_lane1,collapse=", "),". They are now removed for gel 2d dewarping. =="))
      data_dewarp <- subset(data_dewarp,lane_ID!=1)
    }
    
    # Record the settings of current analysis:
    cat("==[spotgear] Results stored in: ==","\n",mcmc_options$result.folder,"\n")
    ##model_options:
    #dput(model_options,file.path(mcmc_options$result.folder,"model_options.txt"))
    #mcmc_options:
    dput(mcmc_options,file.path(mcmc_options$result.folder,"mcmc_options.txt"))
    
    analysis_id <- as.numeric(levels(unique(data_dewarp$gel_ID))[as.numeric(unique(data_dewarp$gel_ID))])
    
    
    #
    # Prepare Data:
    #
    curr_dat      <- data_dewarp # <--------- curr_dat is got rid of "Lane 1"s with known molecules.
    Y             <- curr_dat$Y
    curr_gel_lane <- unique(curr_dat[,c("gel_ID","lane_ID")])
    N             <- nrow(curr_gel_lane)
    n_vec         <- rep(NA,N)
    for (gl in seq_along(n_vec)){
      n_vec[gl] <- nrow(subset(curr_dat,gel_ID==curr_gel_lane[gl,1] & lane_ID==curr_gel_lane[gl,2]))
    }
    
    if (sum(n_vec!=0)!=N){stop("==[gel] some individuals have no peaks detected!==")}
    
    curr_dat$index <- 1:nrow(curr_dat)
    lookup <- as.matrix(reshape2::dcast(curr_dat,gel_ID+lane_ID~band_ID,value.var = "index"))[,-c(1,2)]
    class(lookup) <- "numeric"
    
    
    N_gel     <- length(analysis_id)
    G <- N_gel # sorry, repeated :).
    # get the number of lanes per gel (excluding Lane 1s):
    N_per_gel <- sapply(1:N_gel,function(i)  sum(unique(curr_dat[,c("gel_ID","lane_ID")])$gel_ID==analysis_id[i])) # <--- minus one for reference lane.
    
    #
    # for each gel construct bases for 2-dimensional warping functions:
    #
    
    # vertical:
    T1 <- 6  # <-- no. of bases in "lane" direction.
    
    ZBu <- array(0,c(max(N_per_gel),T1,G))
    u_obs0_list <- u_obs_list <- list()
    for (g in 1:G){
      N_homo   <- N_per_gel[g] # <-- assume equal lane numbers.
      u_obs0   <- (1:N_homo)/(N_homo+1)
      u_obs    <- (u_obs0-mean(u_obs0))/sd(u_obs0) # <-- vertical standardization.
      
      knots_u  <- quantile(u_obs,seq(0,1,length = (T1 - 2))[-c(1,(T1 - 2))])
      ZBu[1:N_homo,,g]      <- matrix(splines::bs(u_obs,knots= knots_u,
                                                  intercept=TRUE),nrow=length(u_obs))
      u_obs0_list[[g]] <- u_obs0
      u_obs_list[[g]] <- u_obs
    }
    
    # horizontal:
    T2 <- 10  # <-- no. of bases in "gel" direction.
    L  <- 100  # <-- no. of interior landmarks.
    Y_std       <- (Y-mean(curr_dat$Y))/sd(curr_dat$Y) # <-- horizontal standardization.
    leftend_std <- (0-mean(curr_dat$Y))/sd(curr_dat$Y)
    rightend_std <- (1-mean(curr_dat$Y))/sd(curr_dat$Y)
    # landmark grids (plus two boundary points):
    v_landmark <- seq(leftend_std,rightend_std,length.out=L+2)
    knots_v    <- quantile(Y_std,seq(0,1,length = (T2 - 2))[-c(1,(T2 - 2))])
    ZBv        <- matrix(splines::bs(v_landmark,knots= knots_v,
                                     intercept=TRUE),nrow=length(v_landmark))
    
    # get coefficients from identical transformation:
    ident_fit  <- lm(v_landmark~-1+ZBv)
    beta_ident <- ident_fit$coef
    
    h          <- diff(v_landmark)[1] # <- difference in distance between neighboring grid points.
    
    if (L==100){
      grid_lb    <- pmax(v_landmark[1],v_landmark-1.5*h*2)
      grid_ub    <- pmin(v_landmark[L+2],v_landmark+1.5*h*2)
    } 
    if (L==50){
      grid_lb    <- pmax(v_landmark[1],v_landmark-1.5*h)
      grid_ub    <- pmin(v_landmark[L+2],v_landmark+1.5*h)
    }
    
    prec_beta1 <- rep(1/(min(diff(beta_ident))/3)^2,G)
    
    N_all_gel <- sum(N_per_gel)
    gel_id    <- as.numeric(curr_gel_lane[,1])
    in_gel_id <- c(unlist(sapply(N_per_gel,function(n) 1:n)))
    
    # in_data:
    in_data <- c("Y_std","N_all_gel","G","gel_id","in_gel_id","N_per_gel",
                 "n_vec","lookup","L","v_landmark","h",
                 "beta_ident","grid_lb","grid_ub","prec_beta1",
                 "ZBu","T1","ZBv","T2")
    
    # out_parameter:
    out_parameter <- c("beta","Z","sigma",
                       "taubeta",
                       "probmat","Lambda","normalized_probmat",
                       "flexible_select",
                       "p_flexible"
    )
    
    # in_init:
    beta_init <- t(replicate(beta_ident,n=T1))
    beta_init[,1] <- beta_init[,T2] <- NA # because the beta's on both ends are fixed.
    beta_init_array <- array(NA,c(nrow(beta_init),ncol(beta_init),G))
    for (g in 1:G){beta_init_array[,,g] <- beta_init}
    
    #probmat_init <- matrix(0,nrow=C,ncol=J+2)
    #probmat_init[,unique(sapply(1:length(Y_std),function(i)
    #  which.min(abs(v_landmark-Y_std[i]))))] <- 1
    in_init <- function(){
      list(beta  = beta_init_array,
           Z0    = sapply(1:length(Y_std),function(i)
             which.min(abs(v_landmark-Y_std[i])))#,
           #subset = subset_init
           #probmat = probmat_init
      )
    }
    
    # write model:
    model_func <- "
    model{ #BEGIN MODEL
	for (i in 1:N_all_gel){
    n_vec[i] ~ dpois(Lambda)
    
    # alignment:
    for (p in 1:n_vec[i]){Z0[lookup[i,p]] ~ dcat(probmat[])} 
    Z[(lookup[i,1]):(lookup[i,n_vec[i]])] <- sort(Z0[(lookup[i,1]):(lookup[i,n_vec[i]])])
    
    Y_std[lookup[i,1]] ~ dnorm(mu[i,1],inv_sigmasq)T(lb[lookup[i,1]],ub[lookup[i,1]])
    lb[lookup[i,1]] <- grid_lb[Z[lookup[i,1]]]
    ub[lookup[i,1]] <- grid_ub[Z[lookup[i,1]]]
    mu[i,1] <- S_grid[in_gel_id[i],Z[lookup[i,1]],gel_id[i]]
    
    for (p in 2:n_vec[i]){
    Y_std[lookup[i,p]] ~ dnorm(mu[i,p],inv_sigmasq)T(lb[lookup[i,p]],ub[lookup[i,p]])
    lb[lookup[i,p]] <- max(Y_std[lookup[i,p-1]],grid_lb[Z[lookup[i,p]]])
    ub[lookup[i,p]] <- grid_ub[Z[lookup[i,p]]] 
    mu[i,p] <- S_grid[in_gel_id[i],Z[lookup[i,p]],gel_id[i]]
    }
  }
    
    # prior for response rates in each subset (shared across gels):
    for (j in 1:(L+2)){
    probmat[j] ~ dnorm(0,inv_scale_mu0[j])T(0,)  # <--- not really sensitivity, but to be scaled to reflect lambda(t)/int lambda(u) du.
    inv_scale_mu0[j] ~ dgamma(10E-4,10E-4)
    }
    Lambda ~ dgamma(0.1,0.1)
    normalized_probmat <- probmat/sum(probmat)*Lambda
    
    for (g in 1:G){
    # warping function:
    #S_grid[1:19,1:(J+2),g] <- ZBu%*%beta[,,g]%*%t(ZBv)
    S_grid[1:N_per_gel[g],1:(L+2),g] <- ZBu[1:N_per_gel[g],,g]%*%beta[,,g]%*%t(ZBv)
    
    #
    # priors and parameters:
    #
    # prior for B-spline coefficients: first-order penaltY_std matrix:
    for (s in 1:T1){ 
    beta[s,1,g]  <- v_landmark[1]
    beta[s,T2,g] <- v_landmark[L+2]
    }
    
    ####### not used; just to make sure the first index is filled with number.
    taubeta[1,g] <- 10000
    flexible_select[1,g] <- 0
    ind_flex_select[1,g] <- 2
    taubeta_inv[1,g] <- 100
    ########
    for (t in 2:(T2-1)){
    beta[1,t,g] ~ dnorm(beta[1,t-1,g]-beta_ident[t-1]+beta_ident[t],prec_beta1[g])T(beta[1,t-1,g],v_landmark[L+2])
    for (s in 2:T1){
    beta[s,t,g] ~ dnorm(beta[s-1,t,g],taubeta[t,g])T(beta[s,t-1,g],v_landmark[L+2]) # constraints.
    }
    # select flexible semiparametric regression:
    taubeta0[t,1,g]    ~ dgamma(3,2)              # <-------- flexible fit.
    taubeta_inv[t,g]   ~ dpar(1.5,1/400)          # <--------constant fit.
    taubeta0[t,2,g]      <- pow(taubeta_inv[t,g],-1)
    flexible_select[t,g] ~ dbern(p_flexible[g])
    ind_flex_select[t,g] <- 2-flexible_select[t,g]
    taubeta[t,g]         <- taubeta0[t,ind_flex_select[t,g],g] # controls the smoothness of the warping function at each v basis.
    }
    #prec_beta1 <- 1600  # controls the warping function prior for the first lane. depends on the number of knots.
    
    p_flexible[g] ~ dbeta(1,1) # fraction of flexible curves (out of total number of gel direction knots).
    }
    inv_sigmasq <- pow(3/h,2) # Measurement error beyond warping.
    
    #inv_sigmasq ~ dgamma(?,?)
    sigma <- pow(inv_sigmasq,-0.5)
    } # END OF MODEL.
    " 
  
    model_bugfile_name <- "model_gel_dewarp.bug"
    
    filename <- file.path(mcmc_options$bugsmodel.dir,model_bugfile_name)
    writeLines(model_func, filename)
    
    here <- environment()
    in_data.list <- lapply(as.list(in_data),get, envir=here)
    names(in_data.list) <- in_data
    dump(names(in_data.list), append = FALSE, envir = here,
         file = file.path(mcmc_options$result.folder,"jagsdata.txt"))
    
    load.module("glm")
    gs <- R2jags::jags2(data   = in_data,                                           # <------------ Bayesian image dewarping.
                        inits  = in_init,
                        parameters.to.save = out_parameter,
                        model.file = filename,
                        working.directory = mcmc_options$result.folder,
                        n.iter         = as.integer(mcmc_options$n.itermcmc),
                        n.burnin       = as.integer(mcmc_options$n.burnin),
                        n.thin         = as.integer(mcmc_options$n.thin),
                        n.chains       = as.integer(mcmc_options$n.chains),
                        DIC            = FALSE,
                        clearWD        = FALSE#,              #<--- special to JAGS.
    )
  }








