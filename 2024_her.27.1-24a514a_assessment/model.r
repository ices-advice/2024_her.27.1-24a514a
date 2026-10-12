# ----------------------------------------------------------------------#
# Remove everything from memory to avoide conflicts
#-------------------------------------------------------------------------------#
rm(list=ls())


#-------------------------------------------------------------------------------#
# Set the working directory to location of the script. 
# (it requieres a new rstudio version)
#-------------------------------------------------------------------------------#
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))



################################################################################
#Prepare SAM config
################################################################################
library(stockassessment)
icesTAF::mkdir('model')
load('input/data.Rda')
conf = loadConf(dat,"bootstrap/model/model.cfg")
par<-defpar(dat,conf)
par$logSdLogN = c(-0.35, -5)
map = list(logSdLogN = as.factor(c(0,NA)))
fit <- sam.fit(dat,conf,par=par,map=map)
save(fit,file='model/fit_official.Rda')


#Run a 5 year retro
ret <- retro(fit,year = 5)
save(ret,file = 'model/retro.Rda')



#Run the residual
res <- resid(fit)
save(res,file='model/residual_official.Rda')


#Due to the assessment uses fixed parameters, and map we need to construct our 
#own leaveout
leaveOutXSAM<-function(fit){
  fleet = as.list(2:fit$data$noFleets)
  runs0 <- lapply(fleet, function(f) runwithout(fit, fleet = f))
  runs <- list()
  for(i in 1:length(runs0)){
    dat0<-runs0[[i]]$data
    conf0<-runs0[[i]]$conf
    par0<-defpar(dat0,conf0)
    par0$logSdLogN = c(-0.35, -5)
    map = list(logSdLogN = as.factor(c(0,NA)))
    fit0 <- sam.fit(dat0,conf0,par=par0,map=map)
    runs[[i]]<-fit0
  }
  names(runs) <- paste0("w.o. ", lapply(fleet, function(x) paste(attr(fit$data, 
                                                                      "fleetNames")[x], collapse = " and ")))
  attr(runs, "fit") <- fit
  class(runs) <- "samset"
  return(runs)
}
LO<-leaveOutXSAM(fit)
save(LO,file='model/lo_official.Rda')





