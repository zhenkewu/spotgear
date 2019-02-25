#' Identify peaks based on the ridges in 2-D CWT coefficient matrix
#'
#' Identify the peaks based on the ridge list (returned by getRidge)
#' in 2-D CWT coefficient matrix and estimated Signal to Noise Ratio
#' (SNR)
#'
#'
#' @param ms	The autoradiointensity data
#' @param ridgeList Returned by \code{\link[MassSpecWavelet]{getRidge}}
#' @param wCoefs 2-D CWT coefficients
#' @param scales The scales of the wavelet function to be matched to the signal.
#' Scales of CWT, by default it is the colnames of wCoefs
#' @param SNR.Th Threshold of SNR
#' @param peakScaleRange The CWT scale range of the peak.
#' @param ridgeLength The maximum ridge scale of the major peaks.
#' @param nearbyPeak Determine whether to include the small peaks close to large major peaks
#' @param nearbyWinSize The window size to determine the nearby peaks. Only effective when nearbyPeak is true.
#' @param winSize.noise The local window size to estimate the noise level.
#' @param SNR.method Method to estimate noise level. Currently, only 95 percentage quantile is supported.
#' @param minNoiseLevel The minimum noise level used in calculating SNR, i.e.,
#' if the estimated noise level is less than "minNoiseLevel", it will use "minNoiseLevel" instead. If the noise level is less than 0.5, it will be treated as the ratio to the maximum amplitude of the spectrum.
#' @param shape_ratio A peak should have peak area at least shape_ratio of the maximum one. Because the
#' same sample could be run on multiple gels with different exposure lengths, we
#' pick peaks relative to the maximum peak values (using shape_ratio).
#' @param minPeakArea The minimum area for a peak if to call a location a peak.
#' An important peak should have non-ignorable peak area
#'
#' @return Return a list with following elements:
#' \itemize{
#' \item \code{peakIndex} The indexes of the identified peaks
#' \item \code{peakValue} For identified peaks, we record the peak values to be used as peak shape information.
#' We will combine this crucial information with other landmark locations' "peak areas"
#' to form the basis for signature estimation in the Bayesian clustering model.
#' \item \code{peakCenterIndex}	 The indexes of peak centers, which correspond to the maximum on the ridge. peakCenterIndex includes all the peaks, not just the identified major peaks.
#' \item \code{peakCenterValue} The CWT coefficients (the maximum on the ridge) corresponding to peakCenterIndex
#' \item \code{peakSNR} The SNR of the peak, which is the ratio of peakCenterValue and noise level
#' \item \code{peakScale} The estimated scale of the peak, which corresponds to the \code{peakCenterIndex}
#' \item \code{potentialPeakIndex} The indexes of all potential peaks, which satisfy all requirements of a peak without considering its SNR. Useful, if you want to change to a lower SNR threshold later.
#' \item \code{allPeakIndex} The indexes of all the peaks, whose order is the same as peakCenterIndex, peakCenterValue, peakSNR and peakScale.
#' }
#' All of these return elements have peak names, which are the same as the corresponding peak ridges. see getRidge for details.
#' @export
spotgear_identifyMajorPeaks <- function (ms, ridgeList, wCoefs,
                                   scales = as.numeric(colnames(wCoefs)),
                                   SNR.Th = 3,
                                   peakScaleRange = c(5,64),
                                   ridgeLength = 32,
                                   nearbyPeak = TRUE,
                                   nearbyWinSize = 100,
                                   winSize.noise = 10,
                                   SNR.method = "quantile",
                                   minNoiseLevel = 0.001,
                                   shape_ratio = 0.01,
                                   minPeakArea = 0.05){
  if (is.null(scales)) {
    scales <- 1:ncol(wCoefs)
    colnames(wCoefs) <- scales
  }else if (is.character(scales)) {
    scales <- as.numeric(scales)
  }
  if (ridgeLength > max(scales))
    ridgeLength <- max(scales)
  if (length(peakScaleRange) == 1) {
    peakScaleRange <- scales[scales >= peakScaleRange]
  }else {
    peakScaleRange <- scales[scales >= peakScaleRange[1] &
                               scales <= peakScaleRange[2]]
  }
  if (minNoiseLevel >= 1)
    names(minNoiseLevel) <- "fixed"
  if (is.null(minNoiseLevel)) {
    minNoiseLevel <- 0
  }else {
    if (is.null(names(minNoiseLevel))) {
      minNoiseLevel <- max(wCoefs) * minNoiseLevel
      #minNoiseLevel <- minNoiseLevel#max(wCoefs) * minNoiseLevel
    }else if (names(minNoiseLevel) != "fixed") {
      minNoiseLevel <- max(wCoefs) * minNoiseLevel
      #minNoiseLevel <- minNoiseLevel#max(wCoefs) * minNoiseLevel
    }
  }
  ridgeLen <- sapply(ridgeList, length)
  ridgeName <- names(ridgeList)
  ridgeInfo <- matrix(as.numeric(unlist(strsplit(ridgeName,
                                                 "_"))), nrow = 2)
  ridgeLevel <- ridgeInfo[1, ]
  notnull <- sapply(ridgeList, function(x) {
    !is.null(x[1])
  })
  mzInd <- sapply(ridgeList[notnull], function(x) {
    x[1]
  })
  ord <- order(mzInd)
  ridgeName <- ridgeName[ord]
  ridgeLen <- ridgeLen[ord]
  ridgeLevel <- ridgeLevel[ord]
  ridgeList <- ridgeList[ord]
  mzInd <- mzInd[ord]
  peakScale <- NULL
  peakCenterInd <- NULL
  peakValue <- NULL
  for (i in 1:length(ridgeList)) {
    ridge.i <- ridgeList[[i]]
    level.i <- ridgeLevel[i]
    levels.i <- level.i:(level.i + ridgeLen[i] - 1)
    scales.i <- scales[levels.i]
    selInd.i <- which(scales.i %in% peakScaleRange)
    if (length(selInd.i) == 0) {
      peakScale <- c(peakScale, scales.i[1])
      peakCenterInd <- c(peakCenterInd, ridge.i[1])
      peakValue <- c(peakValue, 0)
      next
    }
    levels.i <- levels.i[selInd.i]
    scales.i <- scales.i[selInd.i]
    ridge.i <- ridge.i[selInd.i]
    if (scales.i[1] == 0) {
      ind.i <- cbind(ridge.i[-1], levels.i[-1])
    }
    else {
      ind.i <- cbind(ridge.i, levels.i)
    }
    ridgeValue.i <- wCoefs[ind.i]
    maxInd.i <- which.max(ridgeValue.i)
    peakScale <- c(peakScale, scales.i[maxInd.i])
    peakCenterInd <- c(peakCenterInd, ridge.i[maxInd.i])
    peakValue <- c(peakValue, ridgeValue.i[maxInd.i])
  }
  noise <- abs(wCoefs[, "1"])
  peakSNR <- NULL
  nMz <- nrow(wCoefs)
  for (k in 1:length(ridgeList)) {
    ind.k <- mzInd[k]
    start.k <- ifelse(ind.k - winSize.noise < 1, 1, ind.k -
                        winSize.noise)
    end.k <- ifelse(ind.k + winSize.noise > nMz, nMz, ind.k +
                      winSize.noise)
    ms.int <- ms[start.k:end.k]
    noiseLevel.k <- switch(SNR.method, quantile = stats::quantile(noise[start.k:end.k],
                                                           probs = 0.95),
                           sd =stats::sd(noise[start.k:end.k]), mad = stats::mad(noise[start.k:end.k],
                                                                                                                   center = 0), data.mean = mean(ms.int), data.mean.quant = mean(ms.int[ms.int <
                                                                                                                                                                                          stats::quantile(ms.int, probs = 0.95)]))
    if (noiseLevel.k < minNoiseLevel)
      noiseLevel.k <- minNoiseLevel
    peakSNR <- c(peakSNR, peakValue[k]/noiseLevel.k)
  }
  selInd1 <- ((scales[ridgeLevel + ridgeLen - 1] >= ridgeLength) &  (peakValue > max(max(peakValue)*shape_ratio,minPeakArea))) # <-- major peaks should have long ridge AND non-ignorable peak area.
  if (nearbyPeak) {
    selInd1 <- which(selInd1)
    index <- 1:length(mzInd)
    #nearbyWinSize <- 150
    tempInd <- NULL
    for (ind.i in selInd1) {
      tempInd <- c(tempInd, index[mzInd >= mzInd[ind.i] -
                                    nearbyWinSize & mzInd <= mzInd[ind.i] + nearbyWinSize])
    }
    selInd1 <- (index %in% tempInd)
  }
  nearbyWinSize_zw <- 10 # <- control edge effect.
  selInd2 <- (peakSNR > SNR.Th)
  selInd3 <- !(mzInd %in% c(1:(nearbyWinSize_zw/2), (nrow(wCoefs) -
                                                       (nearbyWinSize_zw/2) + 1):nrow(wCoefs)))
  #selInd4 <- (peakValue > max(max(peakValue)*shape_ratio,minPeakArea))  # <---- added shape ratio constraint, so that a peak could only be called when its at least 5% of the largest shape.
  selInd <- (selInd1 & selInd2 & selInd3)
  names(peakSNR) <- names(peakScale) <- names(peakCenterInd) <- names(peakValue) <- names(mzInd) <- ridgeName

  #print("-------")
  #print(ridgeLen)
  return(list(peakIndex = mzInd[selInd], peakValue = peakValue,
              peakCenterIndex = peakCenterInd, peakSNR = peakSNR,
              peakScale = peakScale, potentialPeakIndex = mzInd[selInd1 &
                                                                  selInd3], allPeakIndex = mzInd))
}
