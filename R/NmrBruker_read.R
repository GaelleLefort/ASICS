## RAW BRUKER FILE READING FUNCTION
NmRBrucker_read <- function(DataDir,SampleSpectrum)
{
  
  bruker.get_param <- function (ACQ,paramStr)
  {
    regexpStr <- paste("^...",paramStr,"=",sep="")
    as.numeric(gsub("^[^=]+= ","" ,ACQ[which(simplify2array(regexpr(regexpStr,ACQ))>0)]))
  }
  
  ACQFILE <- file.path(DataDir,"acqus")
  SPECFILE <- file.path(DataDir,"pdata/1/1r")
  PROCFILE <- file.path(DataDir,"pdata/1/procs")
  
  ACQ <- readLines(ACQFILE)
  TD      <- bruker.get_param(ACQ,"TD")
  SW      <- bruker.get_param(ACQ,"SW")
  SWH     <- bruker.get_param(ACQ,"SW_h")
  DTYPA   <- bruker.get_param(ACQ,"DTYPA")
  BYTORDA <- bruker.get_param(ACQ,"BYTORDA")
  #ENDIAN = ifelse( BYTORDA==0, "little", "big")
  ENDIAN <- "little"
  SIZE = ifelse( DTYPA==0, 4, 8)
  
  PROC <- readLines(PROCFILE)
  OFFSET <- bruker.get_param(PROC,"OFFSET")
  SI <- bruker.get_param(PROC,"SI")
  
  to.read = file(SPECFILE,"rb")
  maxTDSI = max(TD,SI)
  #  signal<-rev(readBin(to.read, what="int",size=SIZE, n=TD, signed = TRUE, endian = ENDIAN))
  signal <- rev(readBin(to.read, what="int",size=SIZE, n=maxTDSI, signed = TRUE, endian = ENDIAN))
  close(to.read)
  
  td <- length(signal)
  
  #  dppm <- SW/(TD-1)
  dppm <- SW/(td-1)
  pmax <- OFFSET
  pmin <- OFFSET - SW
  ppmseq <- seq(from=pmin, to=pmax, by=dppm)
  signal <- 100*signal/max(signal)
  
  SampleSpectrum <- cbind(ppmseq,signal)
  return(SampleSpectrum)
}

#SampleDir <- "K:/PROJETS/Metabohub/Bruker/Tlse_BPASourisCerveau/BPA_c21_aq_121-control/1/"
#setwd(SampleDir)
#DataDir <- "pdata/1"
#rawSpectrum <- NmRBrucker_read(DataDir,rawSpectrum)
