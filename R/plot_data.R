# Define some variables shared within the package:
package_env <- new.env(parent = emptyenv())
package_env$fancy_xlab <- expression(paste("Normalized gel location (",t[j]^{(g)},")"))
package_env$colors     <- rep(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                              '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'),2)
package_env$ref_wt     <- c(21.5,31,45,66,97,116,200) # The known molecular weights for the reference lane.

#' Make grayscale heatmap for one gel
#'
#'
#' @param dat_raw Raw autoradiointensity data as obtained from \code{\link{read_IPdata}}
#' @param band_at A vector of integers selected from 1 to the total number of lanes to highlight selected lanes
#' @param ... other graphical parameters
#'
#' @return A heatmap visualizing the autoradiointensity data obtained from a single gel.
#' @export
image_gel <- function(dat_raw,band_at=NULL,...){# dat_raw has the first column being the gel locations.
  # grayscale heatmap:
  ref <- dat_raw[,1]
  n_g <- ncol(dat_raw)-1
  graphics::image(ref,1:n_g,as.matrix(dat_raw[,-1]),
        #axes= TRUE,
        xlab= package_env$fancy_xlab,
        ylab="Lane",
        col = rev(grDevices::grey(seq(0, 1, length = 256))),yaxt="n",...)
  graphics::axis(2,at=1:n_g,labels=1:n_g,las=1,cex.axis=2)
  if (!is.null(band_at)){
    graphics::abline(h=band_at+.5,col="blue",lty=2)
    graphics::abline(h=band_at-.5,col="blue",lty=2)
    graphics::mtext(band_at,4,cex=2,col="blue",las=1)
  }
}



#' Plot intensity of lanes (for raw data prior to cropping)
#'
#' Can be used for binned data
#'
#' @param dat The data for plotting
#' @param mark The handpicked bands
#' @param lanes The index of lanes to be plotted
#' @param threshold Default is 0. Plot data above this threshold value.
#' @param add_mark Default is \code{TRUE}. Add handpicked marks into the plot.
#' @param ... Other graphical parameters.
#' @return  plot with the lanes stacked together that could be marked with
#' hand-picked bands
#' @export
plot_lane <- function(dat,mark,lanes,threshold=0,add_mark=TRUE,...){
  curr_dat0 <- dat
  # plot data:
  curr_dat <- curr_dat0[,-1]
  curr_dat[curr_dat<threshold] <- 0
  if (max(lanes) > ncol(curr_dat)) {
    stop(paste0("==Please remove the lane number: ",
                lanes[which(lanes>ncol(curr_dat))],"==\n"))
  }
  ct <- 0
  for (lane in lanes){
    ct <- ct+1
    ref <- curr_dat0$rf
    if (ct==1){
      graphics::matplot(ref,(as.matrix(curr_dat[,lane])),type="l",
              xlab=package_env$fancy_xlab,
              ylab="Intensity",
              #ylim=c(0,1.4),
              col = package_env$colors[10],
              #col=package_env$colors[lane],
              lwd=2,...)
    }else{
      graphics::matplot(ref,(as.matrix(curr_dat[,lane])),type="l",add=TRUE,
                        #col=package_env$colors[lane],
                        col=package_env$colors[10],
              lwd=2)
    }

    if (add_mark){
      # add marked data:
      curr_marked <- mark
      curr_ind <- which(curr_marked[,1]==lane)
      graphics::points(curr_marked[curr_ind,3],curr_marked[curr_ind,5],type="h",
             #col=package_env$colors[lane],
             col=package_env$colors[10])
    }
  }
  # add molelular weight landmarks:
  #points(ref,(as.matrix(curr_dat0[,2])),type="l",col="grey")
  #legend("topleft",paste("Lane",lanes,sep=": "),
  #       col=colors[lanes],lwd=2,bty="n",horiz=TRUE)
}
