#' Takes any number of R objects as arguments and returns a list whose names are
#' derived from the names of the R objects.
#'
#' Roger Peng's listlabeling challenge from
#' \url{http://simplystatistics.tumblr.com/post/11988685443/computing-on-the-language}.
#' Code copied from \url{https://gist.github.com/ajdamico/1329117/0134148987859856fcecbe4446cfd37e500e4272}
#'
#' @param ... any R objects
#'
#' @return a list as described above
#'
#' @examples
#' #create three example variables for a list
#' x <- 1
#' y <- 2
#' z <- "hello"
#' #display the results
#' make_list( x , y , z )
#' @export
make_list <- function(...) {
  #put all values into a list
  argument_values <- list(...)

  #save all argument names into another list
  argument_names <- as.list(sys.call())

  #cycle through the first list and label with the second, ignoring the function itself
  for (i in 2:length(argument_names)) {
    names(argument_values)[i - 1] <- argument_names[i]
  }
  #return the newly-labeled function
  argument_values
}




