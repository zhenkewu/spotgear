#' Read data and other model information from a folder that stores model results.
#'
#' @param DIR_DEWARP File path to the folder containing posterior samples
#'
#' @return A list with data, options and posterior samples.
#' \itemize{
#' \item \code{bugs.dat}
#' \item \code{model_options}
#' \item \code{clean_otions}
#' \item \code{Nd}; \code{Nu}; \code{Y}; \code{Mobs}; 
#' \item \code{res_dewarp}.
#' }
#' @importFrom R2jags jags2
#'
#' @export

read_dewarp_folder <- function(DIR_DEWARP){
  
  #DIR_DEWARP <- "/Users/zhenkewu/Downloads/gel_bayes_spotgear/dewarp_all/2000\ result\ jags"
  
  #
  # Read data from DIR_DEWARP:
  #
  new_env <- new.env()
  source(file.path(DIR_DEWARP,"jagsdata.txt"),local=new_env)
  bugs.dat <- as.list(new_env)
  rm(new_env)
  res_dewarp <- coda::read.coda(file.path(DIR_DEWARP,"CODAchain1.txt"),
                                file.path(DIR_DEWARP,"CODAindex.txt"),
                                quiet=TRUE)
  for (bugs.variable.name in names(bugs.dat)) {
    assign(bugs.variable.name, bugs.dat[[bugs.variable.name]])
  }
  
  res <- list(bugs.dat = bugs.dat,
              res_dewarp = res_dewarp)
}

