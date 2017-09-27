## Implemented directly from pseudo code in "Wang, K.C., Wang, S.Y., Kuo, C.H.,
#& Tseng Y.J. (2013). Distribution-based classification method for baseline
#correction of metabolomic 1D proton nuclear magnetic resonance spectra.
#Analytical Chemistry, 85(2), 1231-1239."

baseline_corrector <- function(intensity){
  sliding_windows <- 1
  convergence_ratio <- 0.999
  nb_central_points <- 1
  cutoff_threshold <- 1.5
  nb_interpolation_points <- 3
  nb_smooth_points <- 10

  #1 - Calculating the expected value of noise standard deviation
  SDset <- getSD(intensity, sliding_windows)
  sigma <- findNoiseSD(SDset, convergence_ratio)

  #2 & 3 - Classification of windows and spectral points
  SNvector <- isSignal(sigma, SDset, nb_central_points, cutoff_threshold)

  #4 - Baseline fitting
  sStart <- getSignalStart(SNvector)
  sEnd <- sort(1 + length(intensity) - getSignalStart(rev(SNvector)))
  tempBaseline <- getTempBaseline(intensity, sStart, sEnd,
                                  nb_interpolation_points)
  baseline <- smoothBaseline(tempBaseline, nb_smooth_points)

  return(intensity - baseline)
}


## Obtain  a  set  of  standard deviation SDset from intensity with 2w + 1 point
#sliding windows
#' @importFrom stats sd
#' @keywords internal
getSD <- function(intensity, sliding_windows){
  SDset <- numeric(length(intensity))
  for(i in 1:length(intensity)){
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
findNoiseSD <- function(SDset, ratio){
  m1 <- median(SDset)
  SDset <- SDset[SDset < 2 * m1]
  m2 <- median(SDset)

  while(m2 / m1 < ratio){
    m1 <- m2
    SDset <- SDset[SDset < 2 * m1]
    m2 <- median(SDset)
  }
  return(m2)
}

## Use the noise standard deviation sigma to determine whether each data point
#is signal or noise
isSignal <- function(sigma, SDset, windows, threshold){
  SNvector <- SDset * 0

  for(i in 1: length(SNvector)){
    if(SDset[i] > sigma * threshold){
      SNvector[max(1, i - windows):min(i + windows, length(SNvector))] <- 1
    }
  }
  return(SNvector)
}

## Obtain the start point of each signal-fragment
getSignalStart <- function(SNvector){
  sStart <- numeric()

  for(i in 2:length(SNvector)){
    if(SNvector[i] - SNvector[i - 1] > 0){
      sStart <- c(sStart, i)
    }
  }

  if(SNvector[1] == 1){
    sStart <- c(1, sStart)
  }

  return(sStart)
}

## Linear interpolation of each signal segment according to the boundary
getTempBaseline <- function(intensity, sStart, sEnd, windows){
  tempBaseline <- intensity

  for(i in 1:length(sStart)){
    tempBaseline[sStart[i]] <- mean(intensity[max(1, sStart[i] - windows):
                                                min(sStart[i] + windows,
                                                    length(tempBaseline))])

    tempBaseline[sEnd[i]] <- mean(intensity[max(1, sEnd[i] - windows):
                                              min(sEnd[i] + windows,
                                                  length(tempBaseline))])

    for(j in sStart[i]:sEnd[i]){
      tempBaseline[j] <- (intensity[sEnd[i]] - intensity[sStart[i]]) *
        (j - sStart[i]) / (sEnd[i] - sStart[i]) + intensity[sStart[i]]
    }
  }
  return(tempBaseline)
}

## Use 2 * windows + 1 point mean filter over tempBaseline to obtain the
#baseline
smoothBaseline <- function(tempBaseline, windows){
  baseline <- tempBaseline

  for(i in 1:length(tempBaseline)){
    baseline[i] <- mean(tempBaseline[max(1, i - windows):
                                       min(i + windows, length(tempBaseline))])
  }

  return(baseline)
}
