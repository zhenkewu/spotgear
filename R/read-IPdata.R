#' Read Immuno-Precipitation (IP) data
#'
#' Read autoradiographic intensities from multiple lanes in a single gel. NB: add
#' specifics about the technology.
#'
#' @param dir_raw Path to IP data file (\code{.txt})
#' @param dir_raw_marked Default is \code{NULL}. Can supply a path to hand-picked bands from
#' the IP autoradiography (\code{.txt}).
#'
#' @return \itemize{
#' \item \code{dat_raw} A data frame for raw autoragraphical intensities, one row per
#' location along the original gel scan. The first column records the location along
#' the gel that linearly increase from 0 to 1; the second column and beyond are
#' radio-intensity readings obtained from the tested samples. The Lane 1 (in the second column)
#' intensities are observed from the reference sample comprising 7 known
#' molecules of distinct weights (200, 116, 97, 66, 45, 31, 21.5 kDa).
#'
#' \item \code{dat_raw_marked} A data frame with one hand-picked band per row.
#' At least three columns are expected:
#' \itemize{
#' \item Lane Index (from 1 to the total no. of tested samples),
#' \item Location along the gel (a value between 0 and 1), and
#' \item Intensity value (positive numeric value)
#' }
#' }
#'
#' @export
read_IPdata <- function(dir_raw, dir_raw_marked=NULL){
  # read in raw band data:
  dat_raw <- utils::read.table(dir_raw,skip=1,sep="\t",header=TRUE)
  dat_raw_marked <- NA
  if (!is.null(dir_raw_marked)){
    # read in marked data:
    dat_raw_marked <- utils::read.table(dir_raw_marked,skip=6,sep="",header=!TRUE,nrow=82,comment.char = "=")
  }
  make_list(dat_raw,dat_raw_marked)
}
