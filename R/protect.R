#' @noRd
mrglib <- function(){
  system.file("mrgsolve", package = "hmctdmr")

}

#' @noRd
stanLib <- function(){
  system.file("stan", package="hmctdmr")

}

#' @noRd
sampleLib <- function(){
  system.file("sample", package="hmctdmr")

}

stanD <- function(data){ 
  nt <- nrow(data)
  iObs <- with(data, (1:nrow(data))[evid == 0])
  nObs <- length(iObs)
  cObs <- with(data, DV[iObs])
  # time <- with(data, time[time != 0])

  return (with(data, list( nt=nt, nObs=nObs, iObs=as.array(iObs), cObs=as.array(cObs), 
                          cmt=as.array(cmt), time=as.array(time), amt=as.array(amt), rate= as.array(rate), 
                          evid=as.array(evid), ii=as.array(ii), addl=as.array(addl), ss=as.array(ss), 
                          WT=as.array(WT), CLCR=as.array(CLCR), LBW=as.array(LBW))) )
}

stanInit <- function(chains) list( ln_VD_NR=log(0.8), 
                                       ln_Km=log(5.0),
                                       ln_Vmax=log(500))