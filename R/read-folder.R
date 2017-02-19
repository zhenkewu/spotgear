#' Read data and other model information from a folder that stores model results.
#'
#' @param DIR_DEWARP File path to the folder containing posterior samples
#'
#' @return A list with data, options and posterior samples.
#' \itemize{
#' \item \code{bugs.dat}
#' \item \code{res_dewarp}.
#' }
#' @importFrom R2jags jags2
#'
#' @export

read_dewarp_folder <- function(DIR_DEWARP){
  #
  # Read data from DIR_DEWARP:
  #
  new_env   <- new.env()
  source(file.path(DIR_DEWARP,"jagsdata.txt"),local=new_env)
  bugs.dat <- as.list(new_env)
  rm(new_env)
  res_dewarp <- coda::read.coda(file.path(DIR_DEWARP,"CODAchain1.txt"),
                                file.path(DIR_DEWARP,"CODAindex.txt"),
                                quiet=TRUE)

  res <- list(bugs.dat = bugs.dat,
              res_dewarp = res_dewarp)
}

