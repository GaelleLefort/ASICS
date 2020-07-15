# adaptation of alignment functions of speaq package

## Perform alignment for one spectra with clustering algorithm
.doAlignment <- function(tarSpec, tarPeakList, refSpec, refPeakList,
                              maxShift = 65, acceptLostPeak = TRUE){

  # List of peaks from both spectra
  peakList <- rbind(refPeakList, tarPeakList)
  peakList$peakLabel <- c(rep(1, nrow(refPeakList)), rep(0, nrow(tarPeakList)))

  # Alignment
  res <- .alignAlgo(refSpec, tarSpec, peakList,
                            1, length(tarSpec), maxShift = maxShift,
                            acceptLostPeak = acceptLostPeak)

  return(res$tarSpec)
}

## Clustering algorithm for spectrum alignment
#' @importFrom stats hclust dist cutree
.alignAlgo <- function(refSpec, tarSpec, peakList, startP,
                               endP, maxShift = 65, acceptLostPeak = FALSE){

  tarSpec_align <- tarSpec
  peakList_align <- peakList

  # max shift
  if (length(maxShift) == 1) maxShift <- rep(maxShift, 2) * c(-1, 1)

  # segment to align
  minpeakList <- min(peakList$peak)
  maxpeakList <- max(peakList$peak)

  startCheckP <- min(peakList$bInf[peakList$peak == minpeakList])
  if (is.na(startCheckP) | startCheckP < 1) startCheckP <- startP

  endCheckP <- max(peakList$bSup[peakList$peak == maxpeakList])
  if (is.na(endCheckP) | endCheckP > length(tarSpec)) endCheckP <- endP

  # do alignment only if length of segment > 2
  if (endCheckP - startCheckP >= 2) {
    maxShift_seg <-
      c(-min(abs(maxShift[1]), endCheckP - startCheckP - 2,
             min(peakList_align$peak[peakList_align$peakLabel == 0]) -
               startCheckP - 2),
        min(maxShift[2], endCheckP - startCheckP - 2,
            endCheckP -
              max(peakList_align$peak[peakList_align$peakLabel == 0]) - 2))
    # shift
    adj <- .findBestShift(refSpec[startCheckP:endCheckP],
                                   tarSpec[startCheckP:endCheckP],
                                   maxShift = maxShift_seg)$stepAdj

    # do shift
    if (adj != 0) {
      if (acceptLostPeak || ((adj < 0 && adj + minpeakList >= startCheckP) ||
                             (adj > 0 && adj + maxpeakList <= endCheckP))) {
        tarSpec_align[startCheckP:endCheckP] <-
          .doShift(tarSpec[startCheckP:endCheckP], adj,
                         tarSpec[startCheckP - 1], tarSpec[endCheckP + 1])
        maxShift <- maxShift - adj

        peakList_align$peak[peakList_align$peakLabel == 0] <-
          peakList_align$peak[peakList_align$peakLabel == 0] + adj
        peakList_align$bInf[peakList_align$peakLabel == 0] <-
          peakList_align$bInf[peakList_align$peakLabel == 0] + adj
        peakList_align$bSup[peakList_align$peakLabel == 0] <-
          peakList_align$bSup[peakList_align$peakLabel == 0] + adj

        keepPeaks <- which(peakList_align$peak > 0 &
                             peakList_align$peak <= length(tarSpec))
        peakList_align <- peakList_align[keepPeaks, ]
      }
    }

    # continue only if it still remains more than three peaks
    if (nrow(peakList_align) >= 3) {
      # clustering
      hc <- hclust(dist(peakList_align$peak), method = "average")
      clusterLabel <- cutree(hc, h = hc$height[length(hc$height) - 1])

      # if more than one cluster
      if (length(unique(clusterLabel)) > 1){

        # cluster 1 need to be before cluster 2 on spectra
        if (max(peakList_align$peak[clusterLabel == 1]) >
            min(peakList_align$peak[clusterLabel == 2])) {
          clusterLabel <- ifelse(clusterLabel == 1, 2, 1)
        }

        maxsubData1 <- max(peakList_align$bSup[clusterLabel == 1])
        minsubData2 <- min(peakList_align$bInf[clusterLabel == 2])

        # cluster bounds
        endP1 <- maxsubData1
        if (is.na(endP1) | endP1 > length(tarSpec_align)) endP1 <- maxsubData1

        startP2 <- minsubData2
        # cluster 1
        if (length(unique(peakList_align$peakLabel[clusterLabel == 1])) > 1) {
          res1 <- .alignAlgo(refSpec = refSpec, tarSpec = tarSpec_align,
                             peakList = peakList_align[clusterLabel == 1, ],
                             startP = startCheckP, endP = endP1,
                             maxShift = maxShift,
                             acceptLostPeak = acceptLostPeak)
          tarSpec_align <- res1$tarSpec

          if (nrow(peakList_align[clusterLabel == 1, ]) ==
              nrow(res1$peakList)) {
            peakList_align[clusterLabel == 1, ] <- res1$peakList
          } else {
            toConcat <- matrix(rep(res1$peakList[1, ],
                                   nrow(peakList_align[clusterLabel == 1]) -
                                     nrow(res1$peakList)),
                               ncol = 4)
            res1$peakList <- rbind(res1$peakList, toConcat)
            peakList_align[clusterLabel == 1, ] <- res1$peakList
          }
        }

        # cluster 2
        if (length(unique(peakList_align$peakLabel[clusterLabel == 2])) > 1) {
          res2 <- .alignAlgo(refSpec = refSpec, tarSpec = tarSpec_align,
                             peakList = peakList_align[clusterLabel == 2, ],
                             startP = startP2, endP = endCheckP,
                             maxShift = maxShift, acceptLostPeak=acceptLostPeak)
          tarSpec_align <- res2$tarSpec

          if (nrow(peakList_align[clusterLabel == 2, ]) ==
              nrow(res2$peakList)) {
            peakList_align[clusterLabel == 2, ] <- res2$peakList
          } else {
            toConcat <- matrix(rep(res2$peakList[1, ],
                                   nrow(peakList_align[clusterLabel == 2]) -
                                     nrow(res2$peakList)),
                               ncol = 4)
            res2$peakList <- rbind(res2$peakList, toConcat)
            peakList_align[clusterLabel == 2, ] <- res2$peakList
          }
        }
      }
    }
  }

  return(list(tarSpec = tarSpec_align, peakList = peakList_align))
}


#' @importFrom stats median fft
.findBestShift <- function (refSpec, tarSpec, maxShift = 0) {

  # max shift
  maxShift[maxShift > length(refSpec)] <- length(refSpec)
  if (length(maxShift) == 1) maxShift <- rep(maxShift, 2) * c(-1, 1)
  if (maxShift[1] > 0) maxShift[1] <- 0
  if (maxShift[2] < 0) maxShift[2] <- 0

  # compute fft
  M <- length(refSpec)
  zeroAdd <- 2^ceiling(log2(M)) - M
  M <- M + zeroAdd
  R <- fft(c(refSpec * 1e6, double(zeroAdd))) *
    Conj(fft(c(tarSpec * 1e6, double(zeroAdd)))) / M
  vals <- Re(fft(R, inverse = TRUE) / length(R))

  corValue <- -1
  maxpos <- 1
  lenVals <- length(vals)

  # find best shift
  if (anyNA(vals)) {
    stepAdj <- 0
  } else {
    maxVals <- c(ifelse(maxShift[1] != 0,
                        max(rev(vals)[seq_len(abs(maxShift[1]))], na.rm = TRUE),
                        NA), max(vals[seq_len(maxShift[2] + 1)], na.rm = TRUE))
    corValue <- max(maxVals, na.rm = TRUE)

    if (corValue < 0.1) {
      stepAdj <- 0
    } else if (which.max(maxVals) == 1) {
      stepAdj <- which.max(rev(vals)[seq_len(abs(maxShift[1]))]) * -1
    } else {
      stepAdj <- which.max(vals[seq_len(maxShift[2] + 1)]) - 1
    }
  }

  return(list(corValue = corValue, stepAdj = stepAdj))
}


# Perform shift of one segment
#' @importFrom stats approx
.doShift <- function(specSeg, shiftStep, bInf, bSup){

  newSegment <- specSeg
  if (shiftStep < 0) {
    if (is.na(bSup)) bSup <- 0
    approxEnd <- approx(x = seq_len(2), y = c(specSeg[length(specSeg)], bSup),
                        n = abs(shiftStep) + 2)$y
    newSegment <- c(specSeg[(abs(shiftStep) + 1):length(specSeg)],
                    approxEnd[2:(length(approxEnd) - 1)])
  } else if (shiftStep > 0) {
    if (is.na(bInf)) bInf <- 0
    approxStart <- approx(x = seq_len(2), y = c(bInf, specSeg[1]),
                          n = abs(shiftStep) + 2)$y
    newSegment <- c(approxStart[2:(length(approxStart) - 1)],
                    specSeg[seq_len(length(specSeg) - shiftStep)])
  }

  return(newSegment)
}


## Find the spectrum of reference
#' @importFrom BiocParallel bplapply
#' @importFrom plyr llply
.findReference <- function (spectra, ncores = 1, verbose) {

  # similarity
  if (verbose) cat("Compute FFT correlations \n")
  simi_matrix <- bplapply(as.list(1:ncol(spectra)),
     function(x) vapply(1:ncol(spectra),
                        function(y) ifelse(x <= y, 0, .computeFFT(spectra[,y],
                                                                  spectra[,x])),
                                             FUN.VALUE = numeric(1)),
                          BPPARAM = .createEnv(ncores, ncol(spectra), verbose))

  simi_matrix <- do.call("cbind", simi_matrix)

  simi_matrix <- simi_matrix + t(simi_matrix)
  diag(simi_matrix) <- 1

  return(which.max(rowSums(simi_matrix)))
}


## Compute LCSS similarity
.computeFFT <- function (refSpec, tarSpec) {

  M <- length(refSpec)
  R <- fft((refSpec - mean(refSpec))/sd(refSpec)) *
    Conj(fft((tarSpec - mean(tarSpec))/sd(tarSpec))) / (M^2)
  vals <- (2^ceiling(log2(M))/length(refSpec)) * Re(fft(R, inverse = TRUE))

  return(vals[1])
}


## Obtain  a  set  of  standard deviation SDset from intensity with 2w + 1 point
#sliding windows
#' @importFrom stats sd
#' @keywords internal
.getSD <- function(intensity, sliding_windows){
  SDset <- numeric(length(intensity))
  for(i in seq_along(intensity)){
    SDset[i] <- sd(intensity[max(1, i - sliding_windows):
                               min(i + sliding_windows, length(intensity))])
  }
  return(SDset)
}

## Calculate the median m1 from SDset. Exclude the elements greater than 2m1
#from SDset, and recalculate the median m2 from SDset. Repeat this routine until
#m2/m1 converges and set m2 as the expected value of the noise standard
#deviation
#' @importFrom stats median
#' @keywords internal
.findNoiseSD <- function(SDset, ratio){
  m1 <- median(SDset)
  SDset <- SDset[SDset < 2 * m1]
  m2 <- median(SDset)

  while(m2 / m1 < ratio){
    m1 <- m2
    SDset <- SDset[SDset < 2 * m1]
    m2 <- median(SDset) + 1e-5
  }
  return(m2)
}

## Use the noise standard deviation sigma to determine whether each data point
#is signal or noise
.isSignal <- function(sigma, SDset, windows, threshold){
  SNvector <- SDset * 0

  for(i in seq_along(SNvector)){
    if(SDset[i] > sigma * threshold){
      SNvector[max(1, i - windows):min(i + windows, length(SNvector))] <- 1
    }
  }
  return(SNvector)
}

## Find peaks to perform a clustering on them
.findPeaks <- function(idx, spectra) {
  sliding_windows <- 10
  convergence_ratio <- 0.999
  nb_central_points <- 5
  cutoff_threshold <- 3
  optimum_windows <- 3

  spectrum <- spectra[, idx]

  # 1 - calculating the expected value of noise standard deviation
  SDset <- .getSD(spectrum, sliding_windows)
  sigma <- .findNoiseSD(SDset, convergence_ratio)

  # 2 & 3 - classification of windows and spectral points
  SNvector <- .isSignal(sigma, SDset, nb_central_points, cutoff_threshold)

  # 4 - find local optimum
  optimum <- .localOptimum(spectrum, windows = optimum_windows)
  optimum <- lapply(optimum, function(x) x[x %in% which(SNvector == 1)])
  optimum_df <- data.frame(idx = as.numeric(unlist(optimum)),
                           opti = c(rep("min", length(optimum$minima)),
                                    rep("max", length(optimum$maxima))))

  # 5 - signal extremities
  signal_ext <- which(SNvector == 1)

  min_extremities <- signal_ext[!((signal_ext - 1) %in% signal_ext)]
  max_extremities <- signal_ext[!((signal_ext + 1) %in% signal_ext)]
  peaks_extremities <- cbind(min_extremities, max_extremities)

  # 6 - clean each signal segment
  for (i in seq_len(nrow(peaks_extremities))) {
    # optimum corresponding to this segment
    opt_seg <- optimum_df[optimum_df$idx >= peaks_extremities[i, 1] &
                            optimum_df$idx <= peaks_extremities[i, 2], ]

    #at least one max is nedded
    if (!("max" %in% opt_seg$opti)) {
      optimum_df <- optimum_df[!(optimum_df$idx %in% opt_seg$idx), ]
      next
    }

    # condition 1: first and last optimum need to be minima
    if (opt_seg$opti[which.min(opt_seg$idx)] != "min") {
      if (opt_seg$idx[which.min(opt_seg$idx)] == peaks_extremities[i, 1]) {
        optimum_df <- rbind(optimum_df, c(peaks_extremities[i, 1] - 1, "min"))
        peaks_extremities[i, 1] <- peaks_extremities[i, 1] - 1
      } else {
        optimum_df <- rbind(optimum_df, c(peaks_extremities[i, 1], "min"))
      }
    }
    if (opt_seg$opti[which.max(opt_seg$idx)] != "min") {
      if (opt_seg$idx[which.max(opt_seg$idx)] == peaks_extremities[i, 2]) {
        optimum_df <- rbind(optimum_df, c(peaks_extremities[i, 2] + 1, "min"))
        peaks_extremities[i, 2] <- peaks_extremities[i, 2] + 1
      } else {
        optimum_df <- rbind(optimum_df, c(peaks_extremities[i, 2], "min"))
      }
    }
    optimum_df$idx <- as.numeric(optimum_df$idx)
    optimum_df <- optimum_df[order(optimum_df$idx), ]
    opt_seg <- optimum_df[optimum_df$idx >= peaks_extremities[i, 1] &
                            optimum_df$idx <= peaks_extremities[i, 2], ]

    # condition 2: difference between consecutive optimum need to be large
    #enough
    for (j in 2:(nrow(opt_seg) - 1)) {
      if (abs(spectrum[opt_seg$idx[j - 1]] - spectrum[opt_seg$idx[j]]) <
          abs(max(spectrum[optimum_df$idx]) -
              min(spectrum[optimum_df$idx])) * 0.001 &
          abs(spectrum[opt_seg$idx[j]] - spectrum[opt_seg$idx[j + 1]]) <
          abs(max(spectrum[optimum_df$idx]) -
              min(spectrum[optimum_df$idx])) * 0.001) {
        optimum_df <- optimum_df[optimum_df$idx != opt_seg$idx[j], ]
      }
    }
    opt_seg <- optimum_df[optimum_df$idx >= peaks_extremities[i, 1] &
                            optimum_df$idx <= peaks_extremities[i, 2], ]

    #at least one max is nedded
    if (!("max" %in% opt_seg$opti)) {
      optimum_df <- optimum_df[!(optimum_df$idx %in% opt_seg$idx), ]
      next
    }

    # condition3: we should have min, max, min, max, ...
    for (j in seq_len(nrow(opt_seg) - 1)) {
      if (opt_seg$opti[j] == opt_seg$opti[j + 1]) {
        if (opt_seg$opti[j] == "min") {
          toRemove <- which.max(c(spectrum[opt_seg$idx[j]],
                                  spectrum[opt_seg$idx[j + 1]])) - 1
          optimum_df <- optimum_df[optimum_df$idx != opt_seg$idx[j +
                                                                   toRemove], ]
        } else {
          toRemove <- which.min(c(spectrum[opt_seg$idx[j]],
                                  spectrum[opt_seg$idx[j + 1]])) - 1
          optimum_df <- optimum_df[optimum_df$idx != opt_seg$idx[j +
                                                                   toRemove], ]
        }
      }
    }
  }

  # 7 - peak data-frame
  peak_table <- data.frame(peak = optimum_df$idx[optimum_df$opti == "max"])
  peak_table <-
    cbind(peak_table,
          bInf = apply(peak_table, 1,
                       function(x) max(optimum_df$idx[optimum_df$idx < x])))
  peak_table <-
    cbind(peak_table,
          bSup = apply(peak_table, 1,
                       function(x) min(optimum_df$idx[optimum_df$idx > x[1]])))

  return(peak_table)
}


## Find local optimum
.localOptimum <- function(x, windows = 1){
  up <- vapply(seq_len(windows), function(n) c(x[-(seq(n))], rep(NA, n)),
               FUN.VALUE = numeric(length(x)))
  down <- vapply(-seq_len(windows),
                 function(n) c(rep(NA,abs(n)),
                               x[-seq(length(x), length(x) - abs(n) + 1)]),
                 FUN.VALUE = numeric(length(x)))
  a <- cbind(x, up, down)
  list(minima = which(apply(a, 1, min) == a[,1] &
                        (up[, 1] != 0 | down[, 1] != 0)),
       maxima = which(apply(a, 1, max) == a[,1]  &
                        (up[, 1] != 0 | down[, 1] != 0)))
}
