#' spotgear: \strong{S}ubset \strong{P}rofiling and \strong{O}rganizing \strong{T}ools 
#' for \strong{G}el \strong{E}lectrophoresis \strong{A}utoradiography in \strong{R}
#' 
#' An R Package for Fitting Bayesian Two-Dimensional Image Dewarping Models and 
#' Estimating Disease Subsets and Signatures
#'
#' Autoimmune diseases, e.g., scleroderma, are human immune system’s
#' responses to autoantigens in which the body produces specific autoantibodies
#' that target these autoantigens but also cause tissue damage. The autoantibody
#' composition is strikingly different among patients, as indicated by the many
#' different radiolabeled patterns obtained from the mass-spectrometry-based
#' technology - gel electrophoresis autoradiograms (GEA). Human recognition of
#' patterns is not optimal when the patterns are composite/complex, or patterns
#' are scattered across a large number of sam- ples. However, multiple sources of
#' error (including irrelevant intensity differences across gels and warping of 
#' the gels) have traditionally precluded automation of pattern discovery using
#' autoradiograms. In this package, we overcome these limitations by novel initial
#' gel data preprocessing (Bayesian two dimensional image dewarping/registration) and then provide methods to
#  detect disease subgroups.
#' 
#' @seealso 
#' \itemize{
#' \item \url{https://github.com/zhenkewu/spotgear} for the source code
#' and system/software requirements to use \code{spotgear} for your data.
#' }
#'
#' @import rjags
#'
#' @section spotgear functions:
#'\code{\link{dewarp2d}}
#'
#' @docType package
#' @name spotgear
NULL
#> NULL

