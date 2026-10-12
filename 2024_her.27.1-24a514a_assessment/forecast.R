
#-------------------------------------------------------------------------------#
# Remove everything from memory to avoide conflicts
#-------------------------------------------------------------------------------#
rm(list=ls())


icesTAF::mkdir('output')
icesTAF::mkdir('output/tables')
icesTAF::mkdir('output/figures')


# #-------------------------------------------------------------------------------#
# # Set the working directory to location of the script. 
# # (it requieres a new rstudio version)
# #-------------------------------------------------------------------------------#
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# 
# management strategy  340489.4  F=0.14  0.099(0.08,0.121) 3295.978(2617.798,4125.677)
# risc ssbchange TACchange Catchchange
# 1 0.021    14.581     -12.7       -35.8


#load packages
source('bootstrap/rscript/SAM_forecast_from_GIT.R')
library(ggplot2)
load('model/fit_official_2023.Rda')
fit_old <- fit
load('model/fit_official.Rda')
load('model/forecast_last.Rda')
forecast_last <- tt_deterministic
source('bootstrap/rscript/forecast_utils.R')
source('bootstrap/model/utils.R')
library('stockassessment')




################################################################################
#some user input
################################################################################
TAC <- 390010#386297                 #The last year advice (TAC).  This is used to compute relative changes in TAC
nation_quota = 446.928  #439.757#451.038#701.353        #The sum of national quota. This is used in the forecast and to compute changes
set.seed(1234)
rec.years = fit$data$years      #Specify the range to pull recruitment from 
rec_from_par = T              #Compute recruitment from parameter ( maybe the above is redundant)
switch_year_range = T         #Flag to initialize a small hack in the forecast function, so it fits the structure of the data
ave.years <- tail(fit$data$years,3) # which years to compute the means





################################################################################
#Reference points
################################################################################
blim <- 2500
bpa <- 3184
flim = 0.291
fmsy <- 0.157
fmp <- 0.14



#------------------------------------------------------------#
#Define Maturity for the next years. This is set by the working group
#------------------------------------------------------------#
fit2 <- fit
# mo <- rbind(c(0,0.0,0.4,0.8,1,1.0,1, 1,  1,  1,  1),#Next year
#             c(0,0.0,0.4,0.8,1,1.0,1, 1,  1,  1,  1),#two years
#             c(0,0.0,0.4,0.8,1,1.0,1, 1,  1,  1,  1))#three years
mo <- rbind(c(0,0.0,0.1,0.8,1,1.0,1, 1,  1,  1,  1),#Next year
            c(0,0.0,0.1,0.6,1,1.0,1, 1,  1,  1,  1),#two years
            c(0,0.0,0.4,0.6,0.9,1.0,1, 1,  1,  1,  1))#three years
rownames(mo)<-c((as.numeric(max(rownames(fit2$data$propM)))+1):(as.numeric(max(rownames(fit2$data$propM)))+nrow(mo)))
fit2$data$propMat<-rbind(fit2$data$propMat,mo)




#------------------------------------------------------------#
#A hack as there is something strange with the dimension in the fit object
#------------------------------------------------------------#
fit2$data$catchMeanWeight<-fit2$data$catchMeanWeight[,,]
fit2$data$landFrac<-fit2$data$landFrac[,,]
fit2$data$landMeanWeight<-fit2$data$landMeanWeight[,,]
fit2$data$propF<-fit2$data$propF[,,]





#------------------------------------------------------------#
#Define MP, to be able to reduce F when SSB is below btrigger
#------------------------------------------------------------#
scenario <- function(blim,bpa,ftarget){
MP<- c()
MP$ssb<-c(seq(0,blim),seq(blim,bpa),seq(bpa,12000))
MP$F <-c(rep(0.05,length(seq(0,blim))),
         (seq(blim,bpa)-blim)*(ftarget-0.05)/(bpa-blim)+0.05,
         rep(ftarget,length(seq(bpa,12000))))
MP <- data.frame(MP)
return(MP)
}


#Write input table
tmp <- data.frame(
  age=min(fit2$data$minAgePerFleet):max(fit2$data$maxAgePerFleet),
  stockno = as.numeric(tail(ntable(fit2),1)),
  natmort=as.numeric(tail(fit2$data$natMor,1)),
  maturity=as.numeric(tail(fit$data$propMat,1)),
  PropM=as.numeric(tail(fit2$data$propM,1)),
  PropF=as.numeric(tail(fit2$data$propF,1)),
  SWeight=as.numeric(round(tail(fit2$data$stockMeanWeight,1),digits = 3)),
  Explot = as.numeric(round(tail(faytable(fit),1),3)),
  CWeight=round(colMeans(tail(fit2$data$catchMeanWeight,3)),digits = 3)
)
write.csv2(tmp,
file='output/tables/4.8.1.1a.txt'
)

tmp$stockno[-1]<-NA
tmp$maturity<-mo[1,]
write.csv2(
  tmp,
  file='output/tables/4.8.1.1b.txt'
)






#------------------------------------------------------------#
#For bookkeeping
#------------------------------------------------------------#
out_data <- c()
out_table <- c()


#------------------------------------------------------------#
#Compute the TAC for next year. 
#The structure is written to mimic the XSAM structure of forecast
#------------------------------------------------------------#
MP<-scenario(blim,bpa,ftarget = fmp)
tt<-forecast_XSAM(fit2,fBarW =  c(NA,fmp,fmp,fmp),
                    catchval.exact = c(nation_quota,NA,NA,NA),ave.years = ave.years,
                    rec.years = rec.years,deterministic = T,
                    rec_from_par = rec_from_par,switch_year_range = switch_year_range)
new_F<-MP$F[MP$ssb==round(median(tt[2][[1]]$ssb))]
tt<-forecast_XSAM(fit2,fBarW =  c(NA,new_F,fmp,fmp),
                    catchval.exact = c(nation_quota,NA,NA,NA),ave.years = ave.years,
                    rec.years = rec.years,deterministic = T,
                    rec_from_par = rec_from_par,switch_year_range = switch_year_range)
new_F2<-MP$F[MP$ssb==round(median(tt[3][[1]]$ssb))]
tt<-forecast_XSAM(fit2,fBarW =  c(NA,new_F,new_F2,fmp),
                    catchval.exact = c(nation_quota,NA,NA,NA),ave.years = ave.years,
                    rec.years = rec.years,deterministic = T,
                    rec_from_par = rec_from_par,switch_year_range = switch_year_range)
new_F3<-MP$F[MP$ssb==round(median(tt[4][[1]]$ssb))]
tt<-forecast_XSAM(fit2,fBarW =  c(NA,new_F,new_F2,new_F3),
                    catchval.exact = c(nation_quota,NA,NA,NA),ave.years = ave.years,
                    rec.years = rec.years,
                    deterministic = T,
                    rec_from_par = rec_from_par,switch_year_range = switch_year_range)
tt_stochastic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                               catchval.exact  =  c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),
                                                    mean(tt[[3]]$catch),mean(tt[[4]]$catch)),
                               ave.years = ave.years,
                               deterministic = F,
                               rec.years = rec.years,
                               rec_from_par = rec_from_par,switch_year_range = switch_year_range,bias_correct = T)
tt_deterministic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                                  catchval.exact  =  c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),
                                                       mean(tt[[3]]$catch),mean(tt[[4]]$catch)),
                                  ave.years = ave.years,
                                  deterministic = T,
                                  rec.years = rec.years,
                                  rec_from_par = rec_from_par,switch_year_range = switch_year_range)




#Handle table 4.8.2.1
write.csv2(
  data.frame(
  'SSB_intermediat' = round(median(tt_deterministic[[1]]$ssb)),
  catch_intermediat = round(median(tt_deterministic[[1]]$catch)*1000),
  'SSB_forecast'=round(median(tt_deterministic[[2]]$ssb)),
  'FbarW_intermediat'=round(median(tt_deterministic[[1]]$fbar),digits = 3),
  'Recruitment'=paste(round(c(mean(tt_deterministic[[1]]$rec),mean(tt_deterministic[[2]]$rec),mean(tt_deterministic[[3]]$rec))/1000,digits = 3),sep = ',')
  ),
  file = 'output/tables/tab4.8.2.1a.txt')


sam_forecast <- SAMforecast_output(tt_stochastic)
sam_forecast<-merge(sam_forecast,
                   data.frame('year'=2023:2026,'CatchSAM'=
                                c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),
                                  mean(tt[[3]]$catch),mean(tt[[4]]$catch)))
                   )
sam_forecast$rationale <- 'MP'
sam_forecast$assessment_year <- max(fit2$data$years)
out_data<-rbind(out_data,sam_forecast)


fbarw <- round(mean(tt_stochastic[[2]]$fbar),digits = 3)
fbarw_q <- as.numeric(round(quantile(tt_stochastic[[2]]$fbar, probs = c(0.05, 0.95)),digits=3))
ssb <- round(mean(tt_stochastic[[3]]$ssb),digits = 3)
ssb_q <- as.numeric(round(quantile(tt_stochastic[[3]]$ssb, probs = c(0.05, 0.95)),digits=3))

out_table <- rbind(out_table,
data.frame(Rationale = 'management strategy',
           Catches = mean(tt_deterministic[[2]]$catch)*1000,
           Basis = 'F=0.14',
           Finn = round(new_F,3),
           FW  = paste0(fbarw,'(',fbarw_q[1],',',fbarw_q[2],')'),
           ssb  = paste0(ssb,'(',ssb_q[1],',',ssb_q[2],')'),
           risc = mean(tt_stochastic[[3]]$ssb<blim),
           ssbchange = round((round(mean(tt_stochastic[[3]]$ssb),digits = 3)-round(mean(tt_stochastic[[2]]$ssb),digits = 3))/round(mean(tt_stochastic[[2]]$ssb),digits = 3)*100,digits = 3),
           TACchange = round((mean(tt_deterministic[[2]]$catch)*1000-TAC)/TAC*100,digits = 1),
           Catchchange = round((mean(tt_deterministic[[2]]$catch)-nation_quota)/nation_quota*100,digits = 1)
           )
)



save(tt_deterministic,file = 'model/forecast.Rda')
save(tt_stochastic,file = 'model/stochastic.Rda')


sam_forecast <- SAMforecast_output(tt_stochastic)
sam_forecast<-merge(sam_forecast,
                    data.frame('year'=2023:2026,'CatchSAM'=
                                 c(mean(tt_deterministic[[1]]$catch),mean(tt_deterministic[[2]]$catch),mean(tt_deterministic[[3]]$catch),mean(tt_deterministic[[4]]$catch))))
sam_forecast$assessment_year <- max(fit$data$years)


#------------------------------------------------------------#
#MSY approach
#------------------------------------------------------------#
fbarW <- c(NA,fmsy,fmsy,fmsy)
for(i in 2:4){
  tt<-forecast_XSAM(fit2,fBarW =  fbarW,
                    catchval.exact = c(nation_quota,NA,NA,NA),ave.years = ave.years,
                    rec.years = rec.years,deterministic = T,
                    rec_from_par = rec_from_par,switch_year_range = switch_year_range)
  fbarW[i]<-fbarW[i]*mean(tt[[i]]$ssb)/bpa
  if(fbarW[i]>fmsy)fbarW[i]<-fmsy
}


tt_stochastic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                             catchval.exact  =  c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),
                                                  mean(tt[[3]]$catch),mean(tt[[4]]$catch)),
                             ave.years = ave.years,
                             deterministic = F,
                             rec.years = rec.years,
                             rec_from_par = rec_from_par,switch_year_range = switch_year_range,bias_correct = T)
tt_deterministic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                                catchval.exact  =  c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),
                                                     mean(tt[[3]]$catch),mean(tt[[4]]$catch)),
                                ave.years = ave.years,
                                deterministic = T,
                                rec.years = rec.years,
                                rec_from_par = rec_from_par,switch_year_range = switch_year_range)


sam_forecast <- SAMforecast_output(tt_stochastic)
sam_forecast<-merge(sam_forecast,
                   data.frame('year'=2023:2026,'CatchSAM'=
                                c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),mean(tt[[3]]$catch),mean(tt[[4]]$catch))))
sam_forecast$rationale <- 'F=0.157'
sam_forecast$assessment_year <- max(fit2$data$years)


fbarw <- round(median(tt_stochastic[[2]]$fbar),digits = 3)
fbarw_q <- as.numeric(round(quantile(tt_stochastic[[2]]$fbar, probs = c(0.05, 0.95)),digits=3))
ssb <- round(median(tt_stochastic[[3]]$ssb),digits = 3)
ssb_q <- as.numeric(round(quantile(tt_stochastic[[3]]$ssb, probs = c(0.05, 0.95)),digits=3))


out_data<-rbind(out_data,sam_forecast)
out_table <- rbind(out_table,
data.frame(Rationale = 'MSY',
           Catches = mean(tt_deterministic[[2]]$catch)*1000,
           Basis = 'F=0.157',
           Finn=round(fbarW[2],3),
           FW  = paste0(fbarw,'(',fbarw_q[1],',',fbarw_q[2],')'),
           ssb  = paste0(ssb,'(',ssb_q[1],',',ssb_q[2],')'),
           risc = mean(tt_stochastic[[3]]$ssb<blim),
           ssbchange = round((round(mean(tt_stochastic[[3]]$ssb),digits = 3)-round(mean(tt_stochastic[[2]]$ssb),digits = 3))/round(mean(tt_stochastic[[2]]$ssb),digits = 3)*100,digits = 3),
           TACchange = round((mean(tt_deterministic[[2]]$catch)*1000-TAC)/TAC*100,digits = 1),
           Catchchange = round((mean(tt_deterministic[[2]]$catch)-nation_quota)/nation_quota*100,digits = 1)
))




sam_forecast <- SAMforecast_output(tt_stochastic)
sam_forecast<-merge(sam_forecast,
                    data.frame('year'=2023:2026,'CatchSAM'=
                                 c(mean(tt_deterministic[[1]]$catch),mean(tt_deterministic[[2]]$catch),mean(tt_deterministic[[3]]$catch),mean(tt_deterministic[[4]]$catch))))
sam_forecast$assessment_year <- max(fit$data$years)
ssb <- as.data.frame(stockassessment::ssbtable(fit2))
ssb$year <- as.numeric(rownames(ssb))
names(ssb)<-c('ssb','ssblow','ssbhigh','year')
ssb<-rbind(ssb,sam_forecast[,c('year','ssb','ssblow','ssbhigh')])


sam_forecast <- SAMforecast_output(tt)
sam_forecast<-merge(sam_forecast,
                    data.frame('year'=2023:2026,'CatchSAM'=
                                 c(mean(tt_deterministic[[1]]$catch),mean(tt_deterministic[[2]]$catch),mean(tt_deterministic[[3]]$catch),mean(tt_deterministic[[4]]$catch))))
sam_forecast$assessment_year <- max(fit$data$years)

  
  
  
#------------------------------------------------------------#
#Zero catch
#Cannot run the forecast with 0 catch, therefore use an F that is close to 0
#------------------------------------------------------------#
tt<-forecast_XSAM(fit2,
                  catchval.exact = c(nation_quota,NA,NA,NA),ave.years = ave.years,
                  fscale = c(NA,0.00000000001,0.00000000001,0.0000000001),
                  rec.years = rec.years,deterministic = F,
                  rec_from_par = rec_from_par,switch_year_range = switch_year_range)
sam_forecast <- SAMforecast_output(tt)
sam_forecast<-merge(sam_forecast,
                   data.frame('year'=2023:2026,'CatchSAM'=
                                c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),mean(tt[[3]]$catch),mean(tt[[4]]$catch))))
sam_forecast$rationale <- 'F=0'
sam_forecast$assessment_year <- max(fit2$data$years)


fbarw <- round(median(tt[[2]]$fbar),digits = 3)
fbarw_q <- as.numeric(round(quantile(tt[[2]]$fbar, probs = c(0.05, 0.95)),digits=3))
ssb <- round(median(tt[[3]]$ssb),digits = 3)
ssb_q <- as.numeric(round(quantile(tt[[3]]$ssb, probs = c(0.05, 0.95)),digits=3))

out_data<-rbind(out_data,sam_forecast)
out_table <- rbind(out_table,
data.frame(Rationale = 'Zero Catch',
           Catches = round(mean(tt[[2]]$catch)),
           Basis = 'F=0',Finn=0,
           FW  = paste0(fbarw,'(',fbarw_q[1],',',fbarw_q[2],')'),
           ssb  = paste0(ssb,'(',ssb_q[1],',',ssb_q[2],')'),
           risc = mean(tt[[3]]$ssb<blim),
           ssbchange = round((round(mean(tt[[3]]$ssb),digits = 3)-round(mean(tt[[2]]$ssb),digits = 3))/round(mean(tt[[2]]$ssb),digits = 3)*100,digits = 3),
           TACchange = round((mean(tt[[2]]$catch)*1000-TAC)/TAC*100,digits = 1),
           Catchchange = round((mean(tt[[2]]$catch)-nation_quota)/nation_quota*100,digits = 1)
)
)






#------------------------------------------------------------#
# Fpa
#------------------------------------------------------------#
tt<-forecast_XSAM(fit2,fBarW =  c(NA,fmsy,fmsy,fmsy),
                  catchval.exact = c(nation_quota,NA,NA,NA),ave.years = ave.years,
                  rec.years = rec.years,deterministic = F,
                  rec_from_par = rec_from_par,switch_year_range = switch_year_range)
tt_stochastic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                             catchval.exact  =  c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),
                                                  mean(tt[[3]]$catch),mean(tt[[4]]$catch)),
                             ave.years = ave.years,
                             deterministic = F,
                             rec.years = rec.years,
                             rec_from_par = rec_from_par,switch_year_range = switch_year_range,bias_correct = T)
tt_deterministic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                                catchval.exact  =  c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),
                                                     mean(tt[[3]]$catch),mean(tt[[4]]$catch)),
                                ave.years = ave.years,
                                deterministic = T,
                                rec.years = rec.years,
                                rec_from_par = rec_from_par,switch_year_range = switch_year_range)


sam_forecast <- SAMforecast_output(tt_stochastic)
sam_forecast<-merge(sam_forecast,
                   data.frame('year'=2023:2026,'CatchSAM'=
                                c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),mean(tt[[3]]$catch),mean(tt[[4]]$catch))))
sam_forecast$rationale <- 'F=Fpa'
sam_forecast$assessment_year <- max(fit2$data$years)


fbarw <- round(median(tt_stochastic[[2]]$fbar),digits = 3)
fbarw_q <- as.numeric(round(quantile(tt_stochastic[[2]]$fbar, probs = c(0.05, 0.95)),digits=3))
ssb <- round(median(tt_stochastic[[3]]$ssb),digits = 3)
ssb_q <- as.numeric(round(quantile(tt_stochastic[[3]]$ssb, probs = c(0.05, 0.95)),digits=3))

out_data<-rbind(out_data,sam_forecast)
out_table <- rbind(out_table,data.frame(Rationale = 'Fpa',
           Catches = mean(tt_stochastic[[2]]$catch)*1000,
           Basis = 'F=Fpa',
           Finn=round(fmsy,3),
           FW  = paste0(fbarw,'(',fbarw_q[1],',',fbarw_q[2],')'),
           ssb  = paste0(ssb,'(',ssb_q[1],',',ssb_q[2],')'),
           risc = mean(tt_stochastic[[3]]$ssb<blim),
           ssbchange = round((round(mean(tt_stochastic[[3]]$ssb),digits = 3)-round(mean(tt_stochastic[[2]]$ssb),digits = 3))/round(mean(tt_stochastic[[2]]$ssb),digits = 3)*100,digits = 3),
           TACchange = round((mean(tt_stochastic[[2]]$catch)*1000-TAC)/TAC*100,digits = 1),
           Catchchange = round((mean(tt_stochastic[[2]]$catch)-nation_quota)/nation_quota*100,digits = 1)
))



###Flim

tt<-forecast_XSAM(fit2,fBarW =  c(NA,flim,flim,flim),
                  catchval.exact = c(nation_quota,NA,NA,NA),ave.years = ave.years,
                  rec.years = rec.years,deterministic = F,
                  rec_from_par = rec_from_par,switch_year_range = switch_year_range,)
tt_stochastic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                             catchval.exact  =  c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),
                                                  mean(tt[[3]]$catch),mean(tt[[4]]$catch)),
                             ave.years = ave.years,
                             deterministic = F,
                             rec.years = rec.years,
                             rec_from_par = rec_from_par,switch_year_range = switch_year_range,bias_correct = T)
tt_deterministic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                                catchval.exact  =  c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),
                                                     mean(tt[[3]]$catch),mean(tt[[4]]$catch)),
                                ave.years = ave.years,
                                deterministic = T,
                                rec.years = rec.years,
                                rec_from_par = rec_from_par,switch_year_range = switch_year_range)

sam_forecast <- SAMforecast_output(tt_stochastic)
sam_forecast<-merge(sam_forecast,
                   data.frame('year'=2023:2026,'CatchSAM'=
                                c(mean(tt_stochastic[[1]]$catch),mean(tt_stochastic[[2]]$catch),
                                  mean(tt_stochastic[[3]]$catch),mean(tt_stochastic[[4]]$catch))))
sam_forecast$rationale <- 'F=Flim'
sam_forecast$assessment_year <- max(fit2$data$years)


fbarw <- round(median(tt_stochastic[[2]]$fbar),digits = 3)
fbarw_q <- as.numeric(round(quantile(tt_stochastic[[2]]$fbar, probs = c(0.05, 0.95)),digits=3))
ssb <- round(median(tt_stochastic[[3]]$ssb),digits = 3)
ssb_q <- as.numeric(round(quantile(tt_stochastic[[3]]$ssb, probs = c(0.05, 0.95)),digits=3))


out_data <- rbind(out_data,sam_forecast)
out_table <- rbind(out_table,
data.frame(Rationale = 'Flim',
           Catches = mean(tt_stochastic[[2]]$catch)*1000,
           Basis = 'F=0.291',
           Finn= flim,
           FW  = paste0(fbarw,'(',fbarw_q[1],',',fbarw_q[2],')'),
           ssb  = paste0(ssb,'(',ssb_q[1],',',ssb_q[2],')'),
           risc = mean(tt_stochastic[[3]]$ssb<blim),
           ssbchange = round((round(mean(tt_stochastic[[3]]$ssb),digits = 3)-round(mean(tt_stochastic[[2]]$ssb),digits = 3))/round(mean(tt_stochastic[[2]]$ssb),digits = 3)*100,digits = 3),
           TACchange = round((mean(tt_stochastic[[2]]$catch)*1000-TAC)/TAC*100,digits = 1),
           Catchchange = round((mean(tt_stochastic[[2]]$catch)-nation_quota)/nation_quota*100,digits = 1)
))






#SSB = Blim
progress <- list()
objective_function <- function(f_target,b_target) {
  tt <- forecast_XSAM(fit2, fBarW = c(NA, f_target, f_target, f_target),
                      catchval.exact = c(nation_quota, NA, NA, NA),
                      ave.years = ave.years,
                      rec.years = rec.years,
                      deterministic = T,
                      rec_from_par = rec_from_par,
                      switch_year_range = switch_year_range)
  mean_ssb <- median(tt[[3]]$ssb)
  diff <- abs(b_target - mean_ssb)
  progress <<- append(progress, list(list(f_target = f_target, diff = diff)))
  return(diff)
}

# Initial guess for f_target
initial_f_target <- 0.14
result <- optim(par = initial_f_target, 
                fn = objective_function, 
                method = "Brent", 
                lower = 0, upper = 1,
                b_target=blim)  # Assuming f_target is between 0 and 1
optimal_f_target <- result$par

f_target <- round(optimal_f_target,3)
# f_target<-0.447
tt<-forecast_XSAM(fit2,fBarW =  c(NA,f_target,f_target,f_target),
                  catchval.exact = c(nation_quota,NA,NA,NA),ave.years = ave.years,
                  rec.years = rec.years,deterministic = T,
                  rec_from_par = rec_from_par,switch_year_range = switch_year_range,)
tt_stochastic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                             catchval.exact  =  c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),
                                                  mean(tt[[3]]$catch),mean(tt[[4]]$catch)),
                             ave.years = ave.years,
                             deterministic = F,
                             rec.years = rec.years,
                             rec_from_par = rec_from_par,switch_year_range = switch_year_range,bias_correct = T)

# median(tt_stochastic[[3]]$ssb)-blim

tt_deterministic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                                catchval.exact  =  c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),
                                                     mean(tt[[3]]$catch),mean(tt[[4]]$catch)),
                                ave.years = ave.years,
                                deterministic = T,
                                rec.years = rec.years,
                                rec_from_par = rec_from_par,switch_year_range = switch_year_range)

sam_forecast <- SAMforecast_output(tt_stochastic)
sam_forecast<-merge(sam_forecast,
                   data.frame('year'=2023:2026,'CatchSAM'=
                                c(mean(tt_stochastic[[1]]$catch),
                                  mean(tt_stochastic[[2]]$catch),
                                  mean(tt_stochastic[[3]]$catch),
                                  mean(tt_stochastic[[4]]$catch))))
sam_forecast$rationale <- 'SSB=Blim'
sam_forecast$assessment_year <- max(fit2$data$years)


fbarw <- round(median(tt_stochastic[[2]]$fbar),digits = 3)
fbarw_q <- as.numeric(round(quantile(tt_stochastic[[2]]$fbar, probs = c(0.05, 0.95)),digits=3))
ssb <- round(median(tt_stochastic[[3]]$ssb),digits = 3)
ssb_q <- as.numeric(round(quantile(tt_stochastic[[3]]$ssb, probs = c(0.05, 0.95)),digits=3))


out_data <- rbind(out_data,sam_forecast)
out_table <- rbind(out_table,
                   data.frame(Rationale = 'SSB=Blim',
                              Catches = mean(tt[[2]]$catch)*1000,
                              Basis = paste0('F=',fbarw),
                              FW  = paste0(fbarw,'(',fbarw_q[1],',',fbarw_q[2],')'),
                              Finn=f_target,
                              ssb  = paste0(ssb,'(',ssb_q[1],',',ssb_q[2],')'),
                              risc = mean(tt_stochastic[[3]]$ssb<blim),
                              ssbchange = round((round(mean(tt_stochastic[[3]]$ssb),digits = 3)-round(mean(tt_stochastic[[2]]$ssb),digits = 3))/round(mean(tt_stochastic[[2]]$ssb),digits = 3)*100,digits = 3),
                              TACchange = round((mean(tt_stochastic[[2]]$catch)*1000-TAC)/TAC*100,digits = 1),
                              Catchchange = round((mean(tt_stochastic[[2]]$catch)-nation_quota)/nation_quota*100,digits = 1)
                   ))

###SSB = Bpa
# Initial guess for f_target
initial_f_target <- 0.14
result <- optim(par = initial_f_target, 
                fn = objective_function, 
                method = "Brent", 
                lower = 0, upper = 1,b_target=bpa)  # Assuming f_target is between 0 and 1
optimal_f_target <- result$par

f_target <- round(optimal_f_target,3)
tt<-forecast_XSAM(fit2,fBarW =  c(NA,f_target,f_target,f_target),
                  catchval.exact = c(nation_quota,NA,NA,NA),ave.years = ave.years,
                  rec.years = rec.years,deterministic = T,
                  rec_from_par = rec_from_par,switch_year_range = switch_year_range,)
# tt<-forecast_XSAM(fit2,nextssb = c(NA,bpa,bpa,bpa),
#                   catchval.exact = c(nation_quota,NA,NA,NA),ave.years = ave.years,
#                   rec.years = rec.years,deterministic = F,
#                   rec_from_par = rec_from_par,switch_year_range = switch_year_range,)
tt_stochastic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                             catchval.exact  =  c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),
                                                  mean(tt[[3]]$catch),mean(tt[[4]]$catch)),
                             ave.years = ave.years,
                             deterministic = F,
                             rec.years = rec.years,
                             rec_from_par = rec_from_par,switch_year_range = switch_year_range,bias_correct = T)
tt_deterministic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                                catchval.exact  =  c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),
                                                     mean(tt[[3]]$catch),mean(tt[[4]]$catch)),
                                ave.years = ave.years,
                                deterministic = T,
                                rec.years = rec.years,
                                rec_from_par = rec_from_par,switch_year_range = switch_year_range)


sam_forecast <- SAMforecast_output(tt_stochastic)
sam_forecast<-merge(sam_forecast,
                   data.frame('year'=2023:2026,'CatchSAM'=
                                c(mean(tt_stochastic[[1]]$catch),
                                  mean(tt_stochastic[[2]]$catch),
                                  mean(tt_stochastic[[3]]$catch),
                                  mean(tt_stochastic[[4]]$catch))))
sam_forecast$rationale <- 'SSB=Bpa'
sam_forecast$assessment_year <- max(fit2$data$years)


fbarw <- round(median(tt_stochastic[[2]]$fbar),digits = 3)
fbarw_q <- as.numeric(round(quantile(tt_stochastic[[2]]$fbar, probs = c(0.05, 0.95)),digits=3))
ssb <- round(median(tt_stochastic[[3]]$ssb),digits = 3)
ssb_q <- as.numeric(round(quantile(tt_stochastic[[3]]$ssb, probs = c(0.05, 0.95)),digits=3))


out_data <- rbind(out_data,sam_forecast)
out_table <- rbind(out_table,
                   data.frame(Rationale = 'SSB=Bpa',
                              Catches = mean(tt_stochastic[[2]]$catch)*1000,
                              Basis = paste0('F=',fbarw),
                              Finn=f_target,
                              FW  = paste0(fbarw,'(',fbarw_q[1],',',fbarw_q[2],')'),
                              ssb  = paste0(ssb,'(',ssb_q[1],',',ssb_q[2],')'),
                              risc = mean(tt_stochastic[[3]]$ssb<blim),
                              ssbchange = round((round(mean(tt_stochastic[[3]]$ssb),digits = 3)-round(mean(tt_stochastic[[2]]$ssb),digits = 3))/round(mean(tt_stochastic[[2]]$ssb),digits = 3)*100,digits = 3),
                              TACchange = round((mean(tt_stochastic[[2]]$catch)*1000-TAC)/TAC*100,digits = 1),
                              Catchchange = round((mean(tt_stochastic[[2]]$catch)-nation_quota)/nation_quota*100,digits = 1)
                   ))



#Status quo
#SSB = Bpa


tt<-forecast_XSAM(fit2,fBarW =  c(NA,fmp,fmp,fmp),
                  catchval.exact = c(nation_quota,NA,NA,NA),ave.years = ave.years,
                  rec.years = rec.years,deterministic = T,
                  rec_from_par = rec_from_par,switch_year_range = switch_year_range)

f_sq <- head(as.data.frame(fbartable(tt))$Estimate,1)
tt<-forecast_XSAM(fit2,fBarW =  c(NA,f_sq,f_sq,f_sq),
                  catchval.exact = c(nation_quota,NA,NA,NA),ave.years = ave.years,
                  rec.years = rec.years,deterministic = F,
                  rec_from_par = rec_from_par,switch_year_range = switch_year_range)
tt_stochastic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                             catchval.exact  =  c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),
                                                  mean(tt[[3]]$catch),mean(tt[[4]]$catch)),
                             ave.years = ave.years,
                             deterministic = F,
                             rec.years = rec.years,
                             rec_from_par = rec_from_par,switch_year_range = switch_year_range,bias_correct = T)
tt_deterministic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                                catchval.exact  =  c(mean(tt[[1]]$catch),mean(tt[[2]]$catch),
                                                     mean(tt[[3]]$catch),mean(tt[[4]]$catch)),
                                ave.years = ave.years,
                                deterministic = T,
                                rec.years = rec.years,
                                rec_from_par = rec_from_par,switch_year_range = switch_year_range)



sam_forecast <- SAMforecast_output(tt_stochastic)
sam_forecast<-merge(sam_forecast,
                   data.frame('year'=2023:2026,'CatchSAM'=
                                c(mean(tt_stochastic[[1]]$catch),
                                  mean(tt_stochastic[[2]]$catch),
                                  mean(tt_stochastic[[3]]$catch),mean(tt_stochastic[[4]]$catch))))
sam_forecast$rationale <- 'SQ'
sam_forecast$assessment_year <- max(fit2$data$years)


fbarw <- round(median(tt_stochastic[[2]]$fbar),digits = 3)
fbarw_q <- as.numeric(round(quantile(tt_stochastic[[2]]$fbar, probs = c(0.05, 0.95)),digits=3))
ssb <- round(median(tt_stochastic[[3]]$ssb),digits = 3)
ssb_q <- as.numeric(round(quantile(tt_stochastic[[3]]$ssb, probs = c(0.05, 0.95)),digits=3))


out_data <- rbind(out_data,sam_forecast)
out_table <- rbind(out_table,
                   data.frame(Rationale = 'SQ',
                              Catches = mean(tt[[2]]$catch)*1000,
                              Basis = paste0('F=',fbarw),
                              Finn = f_sq,
                              FW  = paste0(fbarw,'(',fbarw_q[1],',',fbarw_q[2],')'),
                              ssb  = paste0(ssb,'(',ssb_q[1],',',ssb_q[2],')'),
                              risc = mean(tt_stochastic[[3]]$ssb<blim),
                              ssbchange = round((round(mean(tt_stochastic[[3]]$ssb),digits = 3)-round(mean(tt_stochastic[[2]]$ssb),digits = 3))/round(mean(tt_stochastic[[2]]$ssb),digits = 3)*100,digits = 3),
                              TACchange = round((mean(tt_stochastic[[2]]$catch)*1000-TAC)/TAC*100,digits = 1),
                              Catchchange = round((mean(tt_stochastic[[2]]$catch)-nation_quota)/nation_quota*100,digits = 1)
                   ))


#Replace basis table with actual input of F values
out_table$Basis <- out_table$Finn
out_table$Finn<-NULL

#print the output tables
print(out_table)
write.csv2(out_table,file='output/tables/tab4.8.2.1b.txt')
save(out_table,file='model/forecast_options.Rda')


