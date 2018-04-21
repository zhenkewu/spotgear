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


#' Compare two partitions (could be of different lengths)
#' @param x cluster indicators for the first partition
#' @param y cluster indicators for the second partition (truth should be put here if any)
#'
#' @return a list with adjusted Rand index (large = agreement), sensitivity and specificity, and Variation
#' of information (0 means perfect agreement) 
#'
#' 
#' @export

myARI <- function (x, y) {
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y)) 
    stop("arguments must be vectors of the same length")
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1))) 
    return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + 
                                                     a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  RI  <- (a+d)/(a+b+c+d)
  
  vi_dist <- function (cl1, cl2, parts = FALSE, base = 2) 
  {
    if (length(cl1) != length(cl2)) 
      stop("cl1 and cl2 must have same length")
    ent <- function(cl) {
      n <- length(cl)
      p <- table(cl)/n
      -sum(p * log(p, base = base))
    }
    mi <- function(cl1, cl2) {
      p12 <- table(cl1, cl2)/length(cl1)
      p1p2 <- outer(table(cl1)/length(cl1), table(cl2)/length(cl2))
      sum(p12[p12 > 0] * log(p12[p12 > 0]/p1p2[p12 > 0], base = base))
    }
    if (!parts) 
      return(ent(cl1) + ent(cl2) - 2 * mi(cl1, cl2))
    ent1 <- ent(cl1)
    ent2 <- ent(cl2)
    mi12 <- mi(cl1, cl2)
    c(vi = ent1 + ent2 - 2 * mi12, `H(1|2)` = ent1 - mi12, `H(2|1)` = ent2 - 
        mi12)
  }
  
  list(ARI=ARI, RI=RI,sensitivity=a/(a+c),specificity = d/(b+d),VI =vi_dist(x,y))
}




#' This function convert ids: (will need to move to the spotgear package.)
#'
#' @param N_g_vec 
#' @param set_id 
#' @param lane_id 
#' @param include_orig_lane1 
#'
#' @return
#' @export
#'
#' @examples
#' N_g_vec     <- c(20,18,20,20,20,20,   20,20,20,9,   20,20,20) # <-- number of lanes per gel.
#' index_start <- c(1,7,11) # the 1st, 7th, 11th gel that represents the first of each cancer group.
#' set_no_vec  <- c(1:6,1:4,1:3) 
#' # index_start <- c(1,7,11,14)
#' # set_no_vec <- c(1:6,1:4,1:3,1:4)
#' PMSCL_machine_lane_id <- c(20,18,16,20,2,20,   20,20,20,NA,   2,20,2,NA,NA,NA,NA) # NB: 2nd one has deleted two empty lanes.
#' # ?? no ref_lanes? Note that we call the first lane to be the reference lanes, sometimes we call them lane 1s.
#' # controls do not have identical bands: (lesson - clear documentation of the variable names after the first test of the model.)
#' # need two ids to identify the control lanes (set id, lane id - prior to removal of the first lane).
#' control_set_id  <- c(6,6,6,6,6,6,6,6,6,6, 
#'                      10,10,10,10,
#'                      13,13,13,13
#' )
#' control_lane_id <- c(10,11,12,13,14,15,16,17,18,19,
#'                      3,4,5,6,
#'                      17,18,19,20
#' )
#' 
#' # the control id if all the lanes are stacked together:
#' ctrl_in_id <- pair2id(N_g_vec,control_set_id,control_lane_id) #ctrl_in_id <- c(102:111+6,171:174+10,231:234+13) <-- including the original lane 1s.
#' 
pair2id <- function(N_g_vec,set_id,lane_id,include_orig_lane1=TRUE){
  if (length(set_id)!=length(lane_id)){stop("==[spotgear] lengths of 'set_id' and 'lane_id' differ. Make them equal. ==")}
  if (length(paste(set_id,lane_id,sep="."))>length(unique(paste(set_id,lane_id,sep=".")))){
    warning("==[spotgear] duplicated (set_id, lane_id) pairs. ==")
    }
  if (include_orig_lane1){
    cat("==[spotgear]ID conversions include the original reference lanes. ==\n")
    res <- c(0,cumsum(N_g_vec))[set_id]+lane_id
  } else{
    cat("==[spotgear]ID conversions DO NOT include the original reference lanes. ==\n")
    res <- cumsum(N_g_vec-1)[set_id] - (N_g_vec[set_id]-lane_id)
  }
  res
} 




