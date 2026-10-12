# Assessment and forecast script for the Norwegian Spring Spawning Herring
# Year: 2024
#
#
#Author. Sindre Vatnehol


# -----------------------------------------------------------------------------#
# Remove everything from memory to avoid conflicts
#------------------------------------------------------------------------------#
rm(list=ls())


#------------------------------------------------------------------------------#
# Set the working directory to location of the script. 
# (it requieres a new rstudio version)
#------------------------------------------------------------------------------#
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))



#------------------------------------------------------------------------------#
#load packages
#------------------------------------------------------------------------------#
library(stockassessment)
library(icesTAF)



#------------------------------------------------------------------------------#
#load and prepare data for SAM 
#------------------------------------------------------------------------------#
# cn<-read.ices("bootstrap/data/caa.dat")
cn <- read.ices("bootstrap/data/caa.txt") #To be updated
cw<-read.ices("bootstrap/data/cw.txt") # to be updated
write.ices(cw,"bootstrap/data/dw.txt")
write.ices(cw,"bootstrap/data/lw.txt")
dw<-read.ices("bootstrap/data/dw.txt")
lf<-read.ices("bootstrap/data/lf.txt")
lw<-read.ices("bootstrap/data/lw.txt")

#The maturity file includes all data from 2 to age 15. but only 2 to 12 is used in the assessment
#Also the file includes data from 1950, but only data from 1988 is used
mo <- read.ices("bootstrap/data/matprop.txt")
mo <- mo[as.numeric(rownames(mo))>=1988,]
mo <- mo[,as.numeric(colnames(mo))%in%c(2:12)]

nm<-read.ices("bootstrap/data/nm.txt")
pf<-read.ices("bootstrap/data/pf.txt")
pm<-read.ices("bootstrap/data/pm.txt")
#read weight in stock
#this dataset is from age 0 to 15+, where this is converted into age 2 to 12+
#in the pluss group the mean values for age >=12 are used
sw<-read.ices("bootstrap/data/west.txt")
sw <- sw[as.numeric(rownames(sw))>=1988,]
sw<-sw[,as.numeric(colnames(sw))>=2]
sw[,as.numeric(colnames(sw))==12]<-rowMeans(sw[,as.numeric(colnames(sw))>=12])
sw <- sw[,as.numeric(colnames(sw))<=12]

surveys<-read.ices("bootstrap/data/survey.dat")
# 

#Source the script to prepare the variance data
source('bootstrap/rscript/prepare_data.r')
load('bootstrap/data/varCAA.dat')
sd_I[is.na(sd_I)]<-1
attributes(cn)$weight = 1/(sd_I**2)
load('bootstrap/data/varFleet1.dat')
sd_I[is.na(sd_I)]<-1
attributes(surveys[[1]])$weight = 1/(sd_I**2)

load('bootstrap/data/varFleet4.dat')
sd_I[is.na(sd_I)]<-1
attributes(surveys[[2]])$weight = 1/(sd_I**2)
load('bootstrap/data/varFleet5.dat')
sd_I2 <- sd_I
sd_I[is.na(sd_I)]<-1
attributes(surveys[[3]])$weight = 1/(sd_I**2)


dat<-setup.sam.data(surveys=surveys,
                    residual.fleet = cn,
                    prop.mature = mo,
                    stock.mean.weight = sw,
                    catch.mean.weight = cw,
                    dis.mean.weight = dw,
                    land.mean.weight = lw,
                    prop.f=pf,
                    prop.m = pm,
                    natural.mortality = nm,
                    land.frac = lf)

icesTAF::mkdir('input')
save(dat,file='input/data.Rda')

# 
# 
# stox_fleet1_var <- read.csv('bootstrap/original_data/caa_var.txt')
# stox_fleet1_var<-stox_fleet1_var[c('year','age','mean','v')]
# stox_fleet1_var$last <- F
# stox_fleet1_var$last[stox_fleet1_var$year==max(stox_fleet1_var$year)]<-T
# stox_fleet1_var$slast <- F
# stox_fleet1_var$slast[stox_fleet1_var$year==max(stox_fleet1_var$year)-1]<-T
# stox_fleet1_var$fleet <- 'Fleet 1'
# out <-stox_fleet1_var
# 
# library(ggplot2)
# ggplot(out,aes(x=log(mean),y=log(v)))+geom_point()+
#   geom_point(data=out[out$last==T,],colour='red')+
#   geom_point(data=out[out$slast==T,],colour='green')
# ggsave('terst.png')
# 
# tt<-read.csv('bootstrap/2023_caa/NbyAge.txt',sep='\t')
# tt$Nor <- tt$CatchAtAge/0.7
# tt[tt$Age==12,]$Nor<-sum(tt[tt$Age>=12,]$Nor)
# tt<-tt[tt$Age<=12,]
# 
# tt<-read.csv('bootstrap/2023_caa/NbyAge.txt',sep='\t')
# tt$SD**2

# tt<-read.csv2(file = 'bootstrap/original_data/ReportBootstrapDataSpawn.txt',sep='\t')
# 
# as.numeric(tt$Abundance_sum_mean)*as.numeric(tt$Abundance_sum_cv)
# as.numeric(tt$Abundance_sum_sd)**2
# 
# dd <- data.frame(en=1,year=2023,age=tt$Age,mean=as.numeric(tt$CatchAtAge),sd=as.numeric(tt$SD))
# dd$cv = dd$sd/dd$mean
# dd$nor <- dd$mean/1000
# (dd$nor*dd$cv)**2
# 
# dd$v = dd$
# write.csv2(data.frame(en=1,year=2023,age=tt$Age,mean=as.numeric(tt$CatchAtAge/1000),sd=as.numeric(tt$SD/1000)**2),
#            file = 'cv_caa.txt')

