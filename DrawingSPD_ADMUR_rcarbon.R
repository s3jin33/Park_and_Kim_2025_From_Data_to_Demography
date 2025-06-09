library(rcarbon)
library(dplyr)
library(RColorBrewer)
library(ADMUR)
library(DEoptimR)
library(dplyr)
library(purrr)
library(ggplot2)
library(maptools)
library(rworldmap)
library(spatstat) 
library("sf")
library(gridExtra)
library(readr)
library(readxl)

getwd()
setwd("/Volumes/SAMSUNG/23-2학회/BHDC/Revision")
setwd("D:/23-2학회/BHDC/분석")


# Reading Data(After Combine ver)
Data<-read.csv("After_Combine_250604.csv")

# Excluding outliers
Data_clean <- Data %>% filter(is.na(Mismatch_Flag))

# Checking no. of site codes-> total 228(regardless of period categorization)
summary(unique(Data$Site_Code))



################# ADMUR: STEP1 Calculating CPL For Entire Han River Basin(Figure 4a) ###########
newData <-Data_clean %>% rename("age"="BP","sd"="error","site"="Site_Code") 
newData<-as.data.frame(newData)


###PD default width = 200. matching with rcarbon
CalArray <- makeCalArray(calcurve=intcal20, calrange = c(1200,3500))
SPD<- summedPhaseCalibrator(data = newData, calcurve=intcal20, calrange=c(1200,3500))
PD <- phaseCalibrator(data=newData, CalArray, remove.external = TRUE)

exp <- JDEoptim(lower=-0.01, upper=0.01, fn=objectiveFunction, PDarray=PD, type='exp', trace=T, NP=20)

fn <- objectiveFunction
CPL1 <- JDEoptim(lower=rep(0,1), upper=rep(1,1), fn, PDarray=PD, type='CPL',trace=T,NP=20)
CPL2 <- JDEoptim(lower=rep(0,3), upper=rep(1,3), fn, PDarray=PD, type='CPL',trace=T,NP=60)
CPL3 <- JDEoptim(lower=rep(0,5), upper=rep(1,5), fn, PDarray=PD, type='CPL',trace=T,NP=100)
CPL4 <- JDEoptim(lower=rep(0,7), upper=rep(1,7), fn, PDarray=PD, type='CPL',trace=T,NP=140)
CPL5 <- JDEoptim(lower=rep(0,9), upper=rep(1,9), fn, PDarray=PD, type='CPL',trace=T,NP=180)
CPL6 <- JDEoptim(lower=rep(0,11),upper=rep(1,11),fn, PDarray=PD, type='CPL',trace=T,NP=220)

save(SPD, PD, exp, CPL1, CPL2, CPL3, CPL4, CPL5, CPL6, file='results.RData',version=2)

load('results.RData')



# Calculate BICs for all six models
# name of each model
model <- c('EXP','1-CPL','2-CPL','3-CPL','4-CPL','5-CPL','6-CPL')

# extract  log likelihoods for each model
loglik <- c(-exp$value, -CPL1$value, -CPL2$value, -CPL3$value, -CPL4$value, -CPL5$value, -CPL6$value)

# extract effective sample sizes
N <- c(rep(ncol(PD),7))

# number of parameters for each model
K <- c(1, 1, 3, 5, 7, 9, 11)

# calculate BIC for each model
BICs <- log(N)*K - 2*loglik

# model comparison. Figure 4-a
png('BICvalueHanRiverBasin.png',width=4000,height=3000,res=500)
par(mar=c(8, 8, 2, 1))
red <- 'firebrick'
blue <- 'steelblue'
col <- rep('grey35',7); col[which(BICs==min(BICs))] <- red
plot(BICs,xlab='',ylab='',xaxt='n', pch=20,cex=5,col=col,main='',las=1,cex.axis=2, font.axis=2)
labels <- c('EXP','1-CPL','2-CPL','3-CPL','4-CPL','5-CPL','6-CPL')
axis(side=1, at=1:7, las=2, labels=labels, cex.axis=2, font=2)
mtext(side=2, at=mean(BICs),text='BIC',las=0,line=6, cex=2, font=2)
dev.off()


## Drawing SPD with ADMUR and rcarbon
CPL1 <- convertPars(pars=CPL1$par, years=1200:3500, type='CPL') 
CPL2 <- convertPars(pars=CPL2$par, years=1200:3500, type='CPL')  
CPL3 <- convertPars(pars=CPL3$par, years=1200:3500, type='CPL')  
CPL4 <- convertPars(pars=CPL4$par, years=1200:3500, type='CPL')  
CPL5 <- convertPars(pars=CPL5$par, years=1200:3500, type='CPL')  
CPL6 <- convertPars(pars=CPL6$par, years=1200:3500, type='CPL')  
EXP <- convertPars(pars=exp$par, years=1200:3500, type='exp')  


#Calibration with rcarbon for rolling mean
Data.caldates=calibrate(x=newData$age, errors=newData$sd, calCurves='intcal20')
Data.bins = rcarbon::binPrep(sites=newData$site,ages=newData$age, h=200)
length(unique(Data.bins)) #453
Data.spd.bins = spd(Data.caldates,bins=Data.bins,timeRange=c(3500,1200), spdnormalised=TRUE)

png('SPD_HanRiverBasin.png',width=9000,height=3000,res=500)
par(mar=c(5,7,1,1), oma=c(4,6,1,1))
years <- as.numeric(row.names(SPD))
plot(NULL,xlim=rev(range(years)), ylim=range(SPD),
     type='l',xaxt='n',ylab='',xlab='',las=1,cex.axis=3,cex.lab=2.5) 
axis(1, at=seq(3500, 1200, by=-100), labels=seq(35, 12, by=-1), cex.axis=2.5, padj=1.2, tck=-0.04, lwd=3, font=2)
axis(2,las=1, cex.axis=3, tck=-0.02, lwd=3, font=2)
mtext(side=1, text='kyr cal BP',line=6.5, cex=2.5, font=2)
mtext(side=2, at=max(SPD[,1])/2,text='PD',las=0,cex=2.5, line=9,font=2)
polygon(c(min(years),years,max(years)),c(0,SPD[,1],0),col="light grey",border="darkgrey")
lines(CPL6$year,CPL6$pdf,col="#5D3A9B",lwd=6)
plot(Data.spd.bins,xaxt='n',runm=100,add=TRUE,type="simple",col="#E66100",lwd=6,lty=2, border="black")
legend("topleft",legend=c("ADMUR", "rcarbon"),
     col=c("#5D3A9B", "#E66100"),lty=c(1,2),lwd=c(6,6), bty='n', cex=2.5)
x_center <- mean(range(years))
y_top <- max(SPD[,1]) 
text(x = x_center, y = y_top * 0.80, labels = "(a) Han River Basin", 
     cex = 2.5, font = 2, pos = 3)
text(x = x_center, y = y_top * 0.65, labels = "N = 2150, bins = 453", 
     cex = 2.3, font = 2, pos = 3)
dev.off()



################# ADMUR: STEP2 Calculating CPL For Northern Gyeonggi(Figure 4b) ###########
NG <- newData %>% filter(Sub_Region == "Northern Gyeonggi") #162

CalArray_NG <- makeCalArray(calcurve=intcal20, calrange = c(1200, 3500))
SPD_NG <- summedPhaseCalibrator(data = NG, calcurve = intcal20, calrange = c(1200, 3500))
PD_NG <- phaseCalibrator(data = NG, CalArray_NG, remove.external = TRUE)

# Model Fitting
exp_NG <- JDEoptim(lower = -0.01, upper = 0.01, fn = objectiveFunction, PDarray = PD_NG, type = 'exp', trace = TRUE, NP = 20)
fn <- objectiveFunction

CPL1_NG <- JDEoptim(lower=rep(0,1), upper=rep(1,1), fn, PDarray=PD_NG, type='CPL', trace=TRUE, NP=20)
CPL2_NG <- JDEoptim(lower=rep(0,3), upper=rep(1,3), fn, PDarray=PD_NG, type='CPL', trace=TRUE, NP=60)
CPL3_NG <- JDEoptim(lower=rep(0,5), upper=rep(1,5), fn, PDarray=PD_NG, type='CPL', trace=TRUE, NP=100)
CPL4_NG <- JDEoptim(lower=rep(0,7), upper=rep(1,7), fn, PDarray=PD_NG, type='CPL', trace=TRUE, NP=140)
CPL5_NG <- JDEoptim(lower=rep(0,9), upper=rep(1,9), fn, PDarray=PD_NG, type='CPL', trace=TRUE, NP=180)
CPL6_NG <- JDEoptim(lower=rep(0,11), upper=rep(1,11), fn, PDarray=PD_NG, type='CPL', trace=TRUE, NP=220)


# Saving results
save(SPD_NG, PD_NG, exp_NG, CPL1_NG, CPL2_NG, CPL3_NG, CPL4_NG, CPL5_NG, CPL6_NG,
     file = 'results_NG.RData', version = 2)


load('results_NG.RData')

# Calculate BICs for all six models (Northern Gyeonggi)
model_NG <- c('EXP','1-CPL','2-CPL','3-CPL','4-CPL','5-CPL','6-CPL')

# extract log likelihoods for each model
loglik_NG <- c(-exp_NG$value, -CPL1_NG$value, -CPL2_NG$value, -CPL3_NG$value, 
               -CPL4_NG$value, -CPL5_NG$value, -CPL6_NG$value)

# extract effective sample sizes
N_NG <- c(rep(ncol(PD_NG), 7))

# number of parameters for each model
K_NG <- c(1, 1, 3, 5, 7, 9, 11)

# calculate BIC for each model
BICs_NG <- log(N_NG) * K_NG - 2 * loglik_NG

# model comparison plot for Northern Gyeonggi Figrue 4-b
png('BICvalue_NorthernGyeonggi.png', width=4000, height=3000, res=500)
par(mar=c(8, 8, 2, 1))
red <- 'firebrick'
col <- rep('grey35', 7)
col[which(BICs_NG == min(BICs_NG))] <- red
plot(BICs_NG, xlab='', ylab='', xaxt='n', pch=20, cex=5, col=col, main='', las=1,
     cex.axis=2, font.axis=2)
labels <- model_NG
axis(side=1, at=1:7, las=2, labels=labels, cex.axis=2, font=2)
mtext(side=2, at=mean(BICs_NG), text='BIC', las=0, line=6, cex=2, font=2)
dev.off()


## Drawing SPD with ADMUR and rcarbon
CPL1_NG <- convertPars(pars=CPL1_NG$par, years=1200:3500, type='CPL') 
CPL2_NG <- convertPars(pars=CPL2_NG$par, years=1200:3500, type='CPL')  
CPL3_NG <- convertPars(pars=CPL3_NG$par, years=1200:3500, type='CPL')  
CPL4_NG <- convertPars(pars=CPL4_NG$par, years=1200:3500, type='CPL')  
CPL5_NG <- convertPars(pars=CPL5_NG$par, years=1200:3500, type='CPL')  
CPL6_NG <- convertPars(pars=CPL6_NG$par, years=1200:3500, type='CPL')  
EXP_NG <- convertPars(pars=exp_NG$par, years=1200:3500, type='exp')  


#Calibration with rcarbon for rolling mean
NG.caldates=calibrate(x=NG$age, errors=NG$sd, calCurves='intcal20')
NG.bins = rcarbon::binPrep(sites=NG$site,ages=NG$age, h=200)
length(unique(NG.bins)) #44
NG.spd.bins = spd(NG.caldates,bins=NG.bins,timeRange=c(3500,1200), spdnormalised=TRUE)

png('SPD_NG_CPL4.png',width=9000,height=3000,res=500)
par(mar=c(5,7,1,1), oma=c(4,6,1,1))
years <- as.numeric(row.names(SPD_NG))
plot(NULL,xlim=rev(range(years)), ylim=range(SPD_NG),
     type='l',xaxt='n',ylab='',xlab='',las=1,cex.axis=3,cex.lab=2.5) 
axis(1, at=seq(3500, 1200, by=-100), labels=seq(35, 12, by=-1), cex.axis=2.5, padj=1.2, tck=-0.04, lwd=3, font=2)
axis(2,las=1, cex.axis=3, tck=-0.02, lwd=3, font=2)
mtext(side=1, text='kyr cal BP',line=6.5, cex=2.5, font=2)
mtext(side=2, at=max(SPD_NG[,1])/2,text='PD',las=0,cex=2.5, line=9,font=2)
polygon(c(min(years),years,max(years)),c(0,SPD_NG[,1],0),col="light grey",border="darkgrey")
lines(CPL4_NG$year,CPL4_NG$pdf,col="#5D3A9B",lwd=6)
plot(NG.spd.bins,xaxt='n',runm=100,add=TRUE,type="simple",col="#E66100",lwd=6,lty=2, border="black")
legend("topleft",legend=c("ADMUR", "rcarbon"),
       col=c("#5D3A9B", "#E66100"),lty=c(1,2),lwd=c(6,6), bty='n', cex=2.5)
text(x = min(years) + 100, y = max(SPD_NG[,1]) * 0.80, labels = "bins = 44", cex = 2.5, font = 2)
dev.off()

png('SPD_NG_CPL3.png',width=9000,height=3000,res=500)
par(mar=c(5,7,1,1), oma=c(4,6,1,1))
years <- as.numeric(row.names(SPD_NG))
plot(NULL,xlim=rev(range(years)), ylim=range(SPD_NG),
     type='l',xaxt='n',ylab='',xlab='',las=1,cex.axis=3,cex.lab=2.5) 
axis(1, at=seq(3500, 1200, by=-100), labels=seq(35, 12, by=-1), cex.axis=2.5, padj=1.2, tck=-0.04, lwd=3, font=2)
axis(2,las=1, cex.axis=3, tck=-0.02, lwd=3, font=2)
mtext(side=1, text='kyr cal BP',line=6.5, cex=2.5, font=2)
mtext(side=2, at=max(SPD_NG[,1])/2,text='PD',las=0,cex=2.5, line=9,font=2)
polygon(c(min(years),years,max(years)),c(0,SPD_NG[,1],0),col="light grey",border="darkgrey")
lines(CPL3_NG$year,CPL3_NG$pdf,col="#5D3A9B",lwd=6)
plot(NG.spd.bins,xaxt='n',runm=100,add=TRUE,type="simple",col="#E66100",lwd=6,lty=2, border="black")
legend("topleft",legend=c("ADMUR", "rcarbon"),
       col=c("#5D3A9B", "#E66100"),lty=c(1,2),lwd=c(6,6), bty='n', cex=2.5)
x_center <- mean(range(years))
y_top <- max(SPD_NG[,1])      
text(x = x_center, y = y_top * 0.80, labels = "(b) Northern Gyeonggi", 
     cex = 2.5, font = 2, pos = 3)
text(x = x_center, y = y_top * 0.65, labels = "N = 162, bins = 44", 
     cex = 2.3, font = 2, pos = 3)
dev.off()

dev.set(dev.next())

################# ADMUR: STEP3 Calculating CPL For Southern Gyeonggi(Figure 4c) ###########
#Filter Data
SG <- newData %>% filter(Sub_Region == "Southern Gyeonggi") #555

CalArray_SG <- makeCalArray(calcurve=intcal20, calrange = c(1200, 3500))
SPD_SG <- summedPhaseCalibrator(data = SG, calcurve = intcal20, calrange = c(1200, 3500))
PD_SG <- phaseCalibrator(data = SG, CalArray_SG, remove.external = TRUE)

exp_SG <- JDEoptim(lower = -0.01, upper = 0.01, fn = objectiveFunction, PDarray = PD_SG, type = 'exp', trace = TRUE, NP = 20)
fn <- objectiveFunction

CPL1_SG <- JDEoptim(lower=rep(0,1), upper=rep(1,1), fn, PDarray=PD_SG, type='CPL', trace=TRUE, NP=20)
CPL2_SG <- JDEoptim(lower=rep(0,3), upper=rep(1,3), fn, PDarray=PD_SG, type='CPL', trace=TRUE, NP=60)
CPL3_SG <- JDEoptim(lower=rep(0,5), upper=rep(1,5), fn, PDarray=PD_SG, type='CPL', trace=TRUE, NP=100)
CPL4_SG <- JDEoptim(lower=rep(0,7), upper=rep(1,7), fn, PDarray=PD_SG, type='CPL', trace=TRUE, NP=140)
CPL5_SG <- JDEoptim(lower=rep(0,9), upper=rep(1,9), fn, PDarray=PD_SG, type='CPL', trace=TRUE, NP=180)
CPL6_SG <- JDEoptim(lower=rep(0,11), upper=rep(1,11), fn, PDarray=PD_SG, type='CPL', trace=TRUE, NP=220)

save(SPD_SG, PD_SG, exp_SG, CPL1_SG, CPL2_SG, CPL3_SG, CPL4_SG, CPL5_SG, CPL6_SG,
     file = 'results_SG.RData', version = 2)


load('results_SG.RData')

# Calculate BICs for all six models (Northern Gyeonggi)
model_SG <- c('EXP','1-CPL','2-CPL','3-CPL','4-CPL','5-CPL','6-CPL')

# extract log likelihoods for each model
loglik_SG <- c(-exp_SG$value, -CPL1_SG$value, -CPL2_SG$value, -CPL3_SG$value, 
               -CPL4_SG$value, -CPL5_SG$value, -CPL6_SG$value)

# extract effective sample sizes
N_SG <- c(rep(ncol(PD_SG), 7))

# number of parameters for each model
K_SG <- c(1, 1, 3, 5, 7, 9, 11)

# calculate BIC for each model
BICs_SG <- log(N_SG) * K_SG - 2 * loglik_SG

png('BICvalue_SouthernGyeonggi.png', width=4000, height=3000, res=500)
par(mar=c(8, 8, 2, 1))
red <- 'firebrick'
col <- rep('grey35', 7)
col[which(BICs_SG == min(BICs_SG))] <- red
plot(BICs_SG, xlab='', ylab='', xaxt='n', pch=20, cex=5, col=col, main='', las=1,
     cex.axis=2, font.axis=2)
axis(side=1, at=1:7, las=2, labels=model_SG, cex.axis=2, font=2)
mtext(side=2, at=mean(BICs_SG), text='BIC', las=0, line=6, cex=2, font=2)
dev.off()


## Drawing SPD with ADMUR and rcarbon
CPL1_SG <- convertPars(pars=CPL1_SG$par, years=1200:3500, type='CPL') 
CPL2_SG <- convertPars(pars=CPL2_SG$par, years=1200:3500, type='CPL')  
CPL3_SG <- convertPars(pars=CPL3_SG$par, years=1200:3500, type='CPL')  
CPL4_SG <- convertPars(pars=CPL4_SG$par, years=1200:3500, type='CPL')  
CPL5_SG <- convertPars(pars=CPL5_SG$par, years=1200:3500, type='CPL')  
CPL6_SG <- convertPars(pars=CPL6_SG$par, years=1200:3500, type='CPL')  
EXP_SG <- convertPars(pars=exp_SG$par, years=1200:3500, type='exp')  


#Calibration with rcarbon for rolling mean
SG.caldates=calibrate(x=SG$age, errors=SG$sd, calCurves='intcal20')
SG.bins = rcarbon::binPrep(sites=SG$site,ages=SG$age, h=200)
length(unique(SG.bins)) #145
SG.spd.bins = spd(SG.caldates,bins=SG.bins,timeRange=c(3500,1200), spdnormalised=TRUE)


png('SPD_SG_6CPL.png',width=9000,height=3000,res=500)
par(mar=c(5,7,1,1), oma=c(4,6,1,1))
years <- as.numeric(row.names(SPD_SG))
plot(NULL,xlim=rev(range(years)), ylim=range(SPD_SG),
     type='l',xaxt='n',ylab='',xlab='',las=1,cex.axis=3,cex.lab=2.5) 
axis(1, at=seq(3500, 1200, by=-100), labels=seq(35, 12, by=-1), cex.axis=2.5, padj=1.2, tck=-0.04, lwd=3, font=2)
axis(2,las=1, cex.axis=3, tck=-0.02, lwd=3, font=2)
mtext(side=1, text='kyr cal BP',line=6.5, cex=2.5, font=2)
mtext(side=2, at=max(SPD_SG[,1])/2,text='PD',las=0,cex=2.5, line=9,font=2)
polygon(c(min(years),years,max(years)),c(0,SPD_SG[,1],0),col="light grey",border="darkgrey")
lines(CPL6_SG$year,CPL6_SG$pdf,col="#5D3A9B",lwd=6)
plot(SG.spd.bins,xaxt='n',runm=100,add=TRUE,type="simple",col="#E66100",lwd=6,lty=2, border="black")
legend("topleft",legend=c("ADMUR", "rcarbon"),
       col=c("#5D3A9B", "#E66100"),lty=c(1,2),lwd=c(6,6), bty='n', cex=2.5)
x_center <- mean(range(years))
y_top <- max(SPD_SG[,1])      
text(x = x_center, y = y_top * 0.80, labels = "(c) Southern Gyeonggi", 
     cex = 2.5, font = 2, pos = 3)
text(x = x_center, y = y_top * 0.65, labels = "N = 555, bins = 145", 
     cex = 2.3, font = 2, pos = 3)
dev.off()

################# ADMUR: STEP4 Calculating CPL For North Han River Basin(Figure 4d) ###########
# Filter Data
NHRB <- newData %>% filter(Sub_Region == "North Han River Basin") #929

CalArray_NHRB <- makeCalArray(calcurve = intcal20, calrange = c(1200, 3500))
SPD_NHRB <- summedPhaseCalibrator(data = NHRB, calcurve = intcal20, calrange = c(1200, 3500))
PD_NHRB <- phaseCalibrator(data = NHRB, CalArray_NHRB, remove.external = TRUE)

exp_NHRB <- JDEoptim(lower = -0.01, upper = 0.01, fn = objectiveFunction, PDarray = PD_NHRB, type = 'exp', trace = TRUE, NP = 20)
fn <- objectiveFunction

CPL1_NHRB <- JDEoptim(lower = rep(0, 1), upper = rep(1, 1), fn, PDarray = PD_NHRB, type = 'CPL', trace = TRUE, NP = 20)
CPL2_NHRB <- JDEoptim(lower = rep(0, 3), upper = rep(1, 3), fn, PDarray = PD_NHRB, type = 'CPL', trace = TRUE, NP = 60)
CPL3_NHRB <- JDEoptim(lower = rep(0, 5), upper = rep(1, 5), fn, PDarray = PD_NHRB, type = 'CPL', trace = TRUE, NP = 100)
CPL4_NHRB <- JDEoptim(lower = rep(0, 7), upper = rep(1, 7), fn, PDarray = PD_NHRB, type = 'CPL', trace = TRUE, NP = 140)
CPL5_NHRB <- JDEoptim(lower = rep(0, 9), upper = rep(1, 9), fn, PDarray = PD_NHRB, type = 'CPL', trace = TRUE, NP = 180)
CPL6_NHRB <- JDEoptim(lower = rep(0, 11), upper = rep(1, 11), fn, PDarray = PD_NHRB, type = 'CPL', trace = TRUE, NP = 220)

save(SPD_NHRB, PD_NHRB, exp_NHRB, CPL1_NHRB, CPL2_NHRB, CPL3_NHRB, CPL4_NHRB, CPL5_NHRB, CPL6_NHRB,
     file = 'results_NHRB.RData', version = 2)

load('results_NHRB.RData')

# Calculate BICs for all six models (North Han River Basin)
model_NHRB <- c('EXP','1-CPL','2-CPL','3-CPL','4-CPL','5-CPL','6-CPL')

loglik_NHRB <- c(-exp_NHRB$value, -CPL1_NHRB$value, -CPL2_NHRB$value, -CPL3_NHRB$value,
                 -CPL4_NHRB$value, -CPL5_NHRB$value, -CPL6_NHRB$value)

N_NHRB <- c(rep(ncol(PD_NHRB), 7))

K_NHRB <- c(1, 1, 3, 5, 7, 9, 11)

BICs_NHRB <- log(N_NHRB) * K_NHRB - 2 * loglik_NHRB

png('BICvalue_NorthHanRiverBasin.png', width = 4000, height = 3000, res = 500)
par(mar = c(11, 8, 2, 1))
red <- 'firebrick'
col <- rep('grey35', 7)
col[which(BICs_NHRB == min(BICs_NHRB))] <- red
plot(BICs_NHRB, xlab = '', ylab = '', xaxt = 'n', pch = 20, cex = 5, col = col, main = '', las = 1,
     cex.axis = 2, font.axis = 2)
axis(side = 1, at = 1:7, las = 2, labels = model_NHRB, cex.axis = 2, font = 2)
mtext(side = 2, at = mean(BICs_NHRB), text = 'BIC', las = 0, line = 6, cex = 2, font = 2)
dev.off()



## Drawing SPD with ADMUR and rcarbon
CPL1_NHRB <- convertPars(pars=CPL1_NHRB$par, years=1200:3500, type='CPL') 
CPL2_NHRB <- convertPars(pars=CPL2_NHRB$par, years=1200:3500, type='CPL')  
CPL3_NHRB <- convertPars(pars=CPL3_NHRB$par, years=1200:3500, type='CPL')  
CPL4_NHRB <- convertPars(pars=CPL4_NHRB$par, years=1200:3500, type='CPL')  
CPL5_NHRB <- convertPars(pars=CPL5_NHRB$par, years=1200:3500, type='CPL')  
CPL6_NHRB <- convertPars(pars=CPL6_NHRB$par, years=1200:3500, type='CPL')  
EXP_NHRB <- convertPars(pars=exp_NHRB$par, years=1200:3500, type='exp')  


#Calibration with rcarbon for rolling mean
NHRB.caldates=calibrate(x=NHRB$age, errors=NHRB$sd, calCurves='intcal20')
NHRB.bins = rcarbon::binPrep(sites=NHRB$site,ages=NHRB$age, h=200)
length(unique(NHRB.bins)) #126
NHRB.spd.bins = spd(NHRB.caldates,bins=NHRB.bins,timeRange=c(3500,1200), spdnormalised=TRUE)



png('SPD_NHRB_6CPL.png',width=9000,height=3000,res=500)
par(mar=c(5,7,1,1), oma=c(4,6,1,1))
years <- as.numeric(row.names(SPD_NHRB))
plot(NULL,xlim=rev(range(years)), ylim=range(SPD_NHRB),
     type='l',xaxt='n',ylab='',xlab='',las=1,cex.axis=3,cex.lab=2.5) 
axis(1, at=seq(3500, 1200, by=-100), labels=seq(35, 12, by=-1), cex.axis=2.5, padj=1.2, tck=-0.04, lwd=3, font=2)
axis(2,las=1, cex.axis=3, tck=-0.02, lwd=3, font=2)
mtext(side=1, text='kyr cal BP',line=6.5, cex=2.5, font=2)
mtext(side=2, at=max(SPD_NHRB[,1])/2,text='PD',las=0,cex=2.5, line=9,font=2)
polygon(c(min(years),years,max(years)),c(0,SPD_NHRB[,1],0),col="light grey",border="darkgrey")
lines(CPL6_NHRB$year,CPL6_NHRB$pdf,col="#5D3A9B",lwd=6)
plot(NHRB.spd.bins,xaxt='n',runm=100,add=TRUE,type="simple",col="#E66100",lwd=6,lty=2, border="black")
legend("topleft",legend=c("ADMUR", "rcarbon"),
       col=c("#5D3A9B", "#E66100"),lty=c(1,2),lwd=c(6,6), bty='n', cex=2.5)
x_center <- mean(range(years))
y_top <- max(SPD_NHRB[,1])      
text(x = x_center, y = y_top * 0.80, labels = "(d) North Han River Basin", 
     cex = 2.5, font = 2, pos = 3)
text(x = x_center, y = y_top * 0.65, labels = "N = 929, bins = 126", 
     cex = 2.3, font = 2, pos = 3)
dev.off()



################# ADMUR: STEP4 Calculating CPL For South Han River Basin(Figure 4e) ###########
# Filter Data
SHRB <- newData %>% filter(Sub_Region == "South Han River Basin") #230

CalArray_SHRB <- makeCalArray(calcurve = intcal20, calrange = c(1200, 3500))
SPD_SHRB <- summedPhaseCalibrator(data = SHRB, calcurve = intcal20, calrange = c(1200, 3500))
PD_SHRB <- phaseCalibrator(data = SHRB, CalArray_SHRB, remove.external = TRUE)

exp_SHRB <- JDEoptim(lower = -0.01, upper = 0.01, fn = objectiveFunction, PDarray = PD_SHRB, type = 'exp', trace = TRUE, NP = 20)
fn <- objectiveFunction

CPL1_SHRB <- JDEoptim(lower = rep(0, 1), upper = rep(1, 1), fn, PDarray = PD_SHRB, type = 'CPL', trace = TRUE, NP = 20)
CPL2_SHRB <- JDEoptim(lower = rep(0, 3), upper = rep(1, 3), fn, PDarray = PD_SHRB, type = 'CPL', trace = TRUE, NP = 60)
CPL3_SHRB <- JDEoptim(lower = rep(0, 5), upper = rep(1, 5), fn, PDarray = PD_SHRB, type = 'CPL', trace = TRUE, NP = 100)
CPL4_SHRB <- JDEoptim(lower = rep(0, 7), upper = rep(1, 7), fn, PDarray = PD_SHRB, type = 'CPL', trace = TRUE, NP = 140)
CPL5_SHRB <- JDEoptim(lower = rep(0, 9), upper = rep(1, 9), fn, PDarray = PD_SHRB, type = 'CPL', trace = TRUE, NP = 180)
CPL6_SHRB <- JDEoptim(lower = rep(0, 11), upper = rep(1, 11), fn, PDarray = PD_SHRB, type = 'CPL', trace = TRUE, NP = 220)

save(SPD_SHRB, PD_SHRB, exp_SHRB, CPL1_SHRB, CPL2_SHRB, CPL3_SHRB, CPL4_SHRB, CPL5_SHRB, CPL6_SHRB,
     file = 'results_SHRB.RData', version = 2)

load('results_SHRB.RData')

# Calculate BICs for all six models (South Han River Basin)
model_SHRB <- c('EXP','1-CPL','2-CPL','3-CPL','4-CPL','5-CPL','6-CPL')

loglik_SHRB <- c(-exp_SHRB$value, -CPL1_SHRB$value, -CPL2_SHRB$value, -CPL3_SHRB$value,
                 -CPL4_SHRB$value, -CPL5_SHRB$value, -CPL6_SHRB$value)

N_SHRB <- c(rep(ncol(PD_SHRB), 7))

K_SHRB <- c(1, 1, 3, 5, 7, 9, 11)

BICs_SHRB <- log(N_SHRB) * K_SHRB - 2 * loglik_SHRB

png('BICvalue_SouthHanRiverBasin.png', width = 4000, height = 3000, res = 500)
par(mar = c(8, 8, 2, 1))
red <- 'firebrick'
col <- rep('grey35', 7)
col[which(BICs_SHRB == min(BICs_SHRB))] <- red
plot(BICs_SHRB, xlab = '', ylab = '', xaxt = 'n', pch = 20, cex = 5, col = col, main = '', las = 1,
     cex.axis = 2, font.axis = 2)
axis(side = 1, at = 1:7, las = 2, labels = model_SHRB, cex.axis = 2, font = 2)
mtext(side = 2, at = mean(BICs_SHRB), text = 'BIC', las = 0, line = 6, cex = 2, font = 2)
dev.off()


## Drawing SPD with ADMUR and rcarbon
CPL1_SHRB <- convertPars(pars=CPL1_SHRB$par, years=1200:3500, type='CPL') 
CPL2_SHRB <- convertPars(pars=CPL2_SHRB$par, years=1200:3500, type='CPL')  
CPL3_SHRB <- convertPars(pars=CPL3_SHRB$par, years=1200:3500, type='CPL')  
CPL4_SHRB <- convertPars(pars=CPL4_SHRB$par, years=1200:3500, type='CPL')  
CPL5_SHRB <- convertPars(pars=CPL5_SHRB$par, years=1200:3500, type='CPL')  
CPL6_SHRB <- convertPars(pars=CPL6_SHRB$par, years=1200:3500, type='CPL')  
EXP_SHRB <- convertPars(pars=exp_SHRB$par, years=1200:3500, type='exp')  


#Calibration with rcarbon for rolling mean
SHRB.caldates=calibrate(x=SHRB$age, errors=SHRB$sd, calCurves='intcal20')
SHRB.bins = rcarbon::binPrep(sites=SHRB$site,ages=SHRB$age, h=200)
length(unique(SHRB.bins)) #70
SHRB.spd.bins = spd(SHRB.caldates,bins=SHRB.bins,timeRange=c(3500,1200), spdnormalised=TRUE)


png('SPD_SHRB.png',width=9000,height=3000,res=500)
par(mar=c(5,7,1,1), oma=c(4,6,1,1))
years <- as.numeric(row.names(SPD_SHRB))
plot(NULL,xlim=rev(range(years)), ylim=range(SPD_SHRB),
     type='l',xaxt='n',ylab='',xlab='',las=1,cex.axis=3,cex.lab=2.5) 
axis(1, at=seq(3500, 1200, by=-100), labels=seq(35, 12, by=-1), cex.axis=2.5, padj=1.2, tck=-0.04, lwd=3, font=2)
axis(2,las=1, cex.axis=3, tck=-0.02, lwd=3, font=2)
mtext(side=1, text='kyr cal BP',line=6.5, cex=2.5, font=2)
mtext(side=2, at=max(SPD_SHRB[,1])/2,text='PD',las=0,cex=2.5, line=9,font=2)
polygon(c(min(years),years,max(years)),c(0,SPD_SHRB[,1],0),col="light grey",border="darkgrey")
lines(CPL5_SHRB$year,CPL5_SHRB$pdf,col="#5D3A9B",lwd=6)
plot(SHRB.spd.bins,xaxt='n',runm=100,add=TRUE,type="simple",col="#E66100",lwd=6,lty=2, border="black")
legend("topleft",legend=c("ADMUR", "rcarbon"),
       col=c("#5D3A9B", "#E66100"),lty=c(1,2),lwd=c(6,6), bty='n', cex=2.5)
x_center <- mean(range(years))
y_top <- max(SPD_SHRB[,1])      
text(x = x_center, y = y_top * 0.80, labels = "(e) South Han River Basin", 
     cex = 2.5, font = 2, pos = 3)
text(x = x_center, y = y_top * 0.65, labels = "N = 230, bins = 70", 
     cex = 2.3, font = 2, pos = 3)
dev.off()

################# ADMUR: STEP5 Calculating CPL For West Coast (Figure 4f) ###########
# Filter Data
WC <- newData %>% filter(Sub_Region == "West Coast") #274

CalArray_WC <- makeCalArray(calcurve = intcal20, calrange = c(1200, 3500))
SPD_WC <- summedPhaseCalibrator(data = WC, calcurve = intcal20, calrange = c(1200, 3500))
PD_WC <- phaseCalibrator(data = WC, CalArray_WC, remove.external = TRUE)

exp_WC <- JDEoptim(lower = -0.01, upper = 0.01, fn = objectiveFunction, PDarray = PD_WC, type = 'exp', trace = TRUE, NP = 20)
fn <- objectiveFunction

CPL1_WC <- JDEoptim(lower = rep(0, 1), upper = rep(1, 1), fn, PDarray = PD_WC, type = 'CPL', trace = TRUE, NP = 20)
CPL2_WC <- JDEoptim(lower = rep(0, 3), upper = rep(1, 3), fn, PDarray = PD_WC, type = 'CPL', trace = TRUE, NP = 60)
CPL3_WC <- JDEoptim(lower = rep(0, 5), upper = rep(1, 5), fn, PDarray = PD_WC, type = 'CPL', trace = TRUE, NP = 100)
CPL4_WC <- JDEoptim(lower = rep(0, 7), upper = rep(1, 7), fn, PDarray = PD_WC, type = 'CPL', trace = TRUE, NP = 140)
CPL5_WC <- JDEoptim(lower = rep(0, 9), upper = rep(1, 9), fn, PDarray = PD_WC, type = 'CPL', trace = TRUE, NP = 180)
CPL6_WC <- JDEoptim(lower = rep(0, 11), upper = rep(1, 11), fn, PDarray = PD_WC, type = 'CPL', trace = TRUE, NP = 220)

save(SPD_WC, PD_WC, exp_WC, CPL1_WC, CPL2_WC, CPL3_WC, CPL4_WC, CPL5_WC, CPL6_WC,
     file = 'results_WC.RData', version = 2)

load('results_WC.RData')

# Calculate BICs for all six models (West Coast)
model_WC <- c('EXP','1-CPL','2-CPL','3-CPL','4-CPL','5-CPL','6-CPL')

loglik_WC <- c(-exp_WC$value, -CPL1_WC$value, -CPL2_WC$value, -CPL3_WC$value,
               -CPL4_WC$value, -CPL5_WC$value, -CPL6_WC$value)

N_WC <- c(rep(ncol(PD_WC), 7))

K_WC <- c(1, 1, 3, 5, 7, 9, 11)

BICs_WC <- log(N_WC) * K_WC - 2 * loglik_WC

png('BICvalue_WestCoast.png', width = 4000, height = 3000, res = 500)
par(mar = c(8, 8, 2, 1))
red <- 'firebrick'
col <- rep('grey35', 7)
col[which(BICs_WC == min(BICs_WC))] <- red
plot(BICs_WC, xlab = '', ylab = '', xaxt = 'n', pch = 20, cex = 5, col = col, main = '', las = 1,
     cex.axis = 2, font.axis = 2)
axis(side = 1, at = 1:7, las = 2, labels = model_WC, cex.axis = 2, font = 2)
mtext(side = 2, at = mean(BICs_WC), text = 'BIC', las = 0, line = 6, cex = 2, font = 2)
dev.off()


## Drawing SPD with ADMUR and rcarbon
CPL1_WC <- convertPars(pars=CPL1_WC$par, years=1200:3500, type='CPL') 
CPL2_WC <- convertPars(pars=CPL2_WC$par, years=1200:3500, type='CPL')  
CPL3_WC <- convertPars(pars=CPL3_WC$par, years=1200:3500, type='CPL')  
CPL4_WC <- convertPars(pars=CPL4_WC$par, years=1200:3500, type='CPL')  
CPL5_WC <- convertPars(pars=CPL5_WC$par, years=1200:3500, type='CPL')  
CPL6_WC <- convertPars(pars=CPL6_WC$par, years=1200:3500, type='CPL')  
EXP_WC <- convertPars(pars=exp_WC$par, years=1200:3500, type='exp')  


#Calibration with rcarbon for rolling mean
WC.caldates=calibrate(x=WC$age, errors=WC$sd, calCurves='intcal20')
WC.bins = rcarbon::binPrep(sites=WC$site,ages=WC$age, h=200)
length(unique(WC.bins)) #69
WC.spd.bins = spd(WC.caldates,bins=WC.bins,timeRange=c(3500,1200), spdnormalised=TRUE)


png('SPD_WC.png',width=9000,height=3000,res=500)
par(mar=c(5,7,1,1), oma=c(4,6,1,1))
years <- as.numeric(row.names(SPD_WC))
plot(NULL,xlim=rev(range(years)), ylim=range(SPD_WC),
     type='l',xaxt='n',ylab='',xlab='',las=1,cex.axis=3,cex.lab=2.5) 
axis(1, at=seq(3500, 1200, by=-100), labels=seq(35, 12, by=-1), cex.axis=2.5, padj=1.2, tck=-0.04, lwd=3, font=2)
axis(2,las=1, cex.axis=3, tck=-0.02, lwd=3, font=2)
mtext(side=1, text='kyr cal BP',line=6.5, cex=2.5, font=2)
mtext(side=2, at=max(SPD_WC[,1])/2,text='PD',las=0,cex=2.5, line=9,font=2)
polygon(c(min(years),years,max(years)),c(0,SPD_WC[,1],0),col="light grey",border="darkgrey")
lines(CPL4_WC$year,CPL4_WC$pdf,col="#5D3A9B",lwd=6)
plot(WC.spd.bins,xaxt='n',runm=100,add=TRUE,type="simple",col="#E66100",lwd=6,lty=2, border="black")
legend("topleft",legend=c("ADMUR", "rcarbon"),
       col=c("#5D3A9B", "#E66100"),lty=c(1,2),lwd=c(6,6), bty='n', cex=2.5)
x_center <- mean(range(years))
y_top <- max(SPD_NHRB[,1])      
text(x = x_center, y = y_top * 0.7, labels = "(f) West Coast", 
     cex = 2.5, font = 2, pos = 3)
text(x = x_center, y = y_top * 0.55, labels = "N = 274, bins = 69", 
     cex = 2.3, font = 2, pos = 3)
dev.off()
