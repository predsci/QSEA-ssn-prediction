# SSN prediction using Q-SEA model
# This code is intended to reproduce Figures 1-7 in the 
# 2023 solar physics paper: 
# "On the Strength and Duration of Solar Cycle 25: 
# A Novel Quantile-based Superposed-Epoch Analysis"
# As such, it's not the most elegantly designed script. 
# It's designed to generate each figure without having to 
# change parameters interactively. 
# An analysis version of this code provides more flexibility
# but requires choosing/setting a number of parameters to 
# get to the results in the paper. 
#
# Please note that some of the later figures require that 
# variables used in earlier figures be previously 
# declared and populated with values. 
# 
# It can - and hopefully will - be further developed 
# to be a more robust package. However, to be 
# generally useful to the community, it should probably
# be converted to a Python package. 

# If you would like to use the logic/approach outlined 
# in this code, please feel free to do so. I'd be 
# grateful if you sent me an email (pete@predsci.com)
# letting me know that you are using it, and please 
# don't hesitate to contact me if you have any questions 
# about the code. 

# The code relies on two generally purpose libraries, 
# both of which can be downloaded via CRAN. 

library(aTSA) # load to test if time series is stationary
library(Metrics) # needed for the MAE calculation 

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# read in the data

# get the data from here
# https://www.sidc.be/silso/datafiles

dir = "/Users/pete/Dropbox/research/ssn-prediction/"
file1 = "SN_m_tot_V2.0.csv"
arr1  = read.csv(paste(dir,file1,sep=""),header=F,sep=";")
dyear = arr1$V3
ssn_m = arr1$V4
file2 = "SN_ms_tot_V2.0.csv"
arr2  = read.csv(paste(dir,file2,sep=""),header=F,sep=";")
ssn_ms = arr2$V4

# identify all the solar minima from the 13-month running average
# start with the rough month of the minimum, then finesse it.

tmin = c(1755+2/12.,
         1766+6/12.,
         1775+6/12.,
         1784+9/12.,
         1798+4/12.,
         1810+8/12.,
         1823+5/12.,
         1833+11/12.,
         1843+7/12.,
         1855+12/12.,
         1867+3/12.,
         1878+12/12.,
         1890+3/12.,
         1902+1/12.,
         1913+7/12.,
         1923+8/12.,
         1933+9/12.,
         1944+2/12.,
         1954+4/12.,
         1964+10/12.,
         1976+3/12.,
         1986+9/12.,
         1996+8/12.,
         2008+12/12.,
         2019+12/12.)

nmin = length(tmin)

# create an interpolated time-series for making Q-SEA model later
dyear1 = 0.1
yearMax = 11.0
time_dec = seq(0,yearMax,dyear1) 
ssn_dec  = matrix(nrow=nmin,ncol=length(time_dec))

# choose which figure from the manuscript to plot
# 1: figure 1 - continuous time series.
# 2: etc, etc.

ssn = ssn_ms

whPlot = 6
#
#---------------------------------------------------------------
#
if (whPlot == 1) {
ssn = ssn_m
plot(dyear,ssn,ylab="SSN",xlab="Time (Years)",type="l")
good13 = (ssn_ms > 0)
lines(dyear[good13],ssn_ms[good13],col="green")

# add vertical lines for cycle boundaries and number 
# each of the cycles

for (i in 1:nmin) {
  lines(c(0,0)+tmin[i],c(0,500),col="red")
  text(tmin[i]+5,350,i)
}
}
#
#---------------------------------------------------------------
#
# make a superposed epoch analysis based on these minima:
#
if (whPlot > 1) {
  ssn_spe = ssn_ms
  tmax = 12
   plot_25 = 'n'
   fit_rise = 'n'
   plot_sea = 'y'
   lw24 = 1


if (whPlot == 2) {
  lw24 = 4
}

if (whPlot == 3) {
  tmax= 4
  ssn_spe = ssn
  fit_rise = 'y'
}

if (whPlot == 4) {
  ssn_spe = ssn_m
  plot_25 = 'y'
  lw24 = 4
}

if (whPlot == 5) {
  plot_25 = 'y'
  lw24 = 4
}

plot(dyear-tmin[1],ssn_spe,xlim=c(0,tmax),ylim=c(0,300),type="l",xlab="Time (Years)",ylab="SSN",col=1) #,
   #  main="13-month running average Sunspot number (except cycle 25, which is monthly)")

col = c("red","blue","green","")
for (i in 1:nmin) {
  lwd = 1
#  if (i != nmin) {
    lwd = 1
    if (i==24) lwd=lw24
    if (i != 25) lines(dyear-tmin[i],ssn_spe,col=i,type="l",lwd=lwd)
#  } else {
    lwd=4
    if (i==25) lines(dyear-tmin[i],ssn_m,col=nmin,type="l",lwd=lwd)
#  }
  ssn_dec[i,] = approx(dyear-tmin[i],ssn_spe,time_dec)$y
}  

# correct interpolation when only one point inside bounds
ltzero = ssn_dec[25,] < 0

ssn_dec[25,ltzero] = NA

legend(0, 300, legend=1:nmin,col=1:nmin, lty=rep(1,nmin), cex=0.65)

# test for stationarity...doesn't appear to be
# dickey_fuller_test = adf.test(ssn,nlag=170)

# Now do the analysis on the relative position of the current cycle: 

rank_dec  = 0.0*ssn_dec

for (j in 1:length(time_dec)) {
   rank_dec[,j] = rank(ssn_dec[,j])
}

# The following statement prints out the relative ranking of solar cycle 
# 25 with respect to the others. All points with "25" are put at the end 
# because they are "NA"...the data hasn't been seen yet.
# the final non-25 value of 2 is also an artefact of the interpolation 
# routine. This number will depend on the number resolution of the 
# interpolation, which is a free varaible in this analysis. 

rank_dec[25,]

# create the predicted time series for cycle 25

ssn25   = 0.0*time_dec
ssn_sea = 0.0*time_dec

for (j in 1:length(time_dec)) {
  quant = quantile(ssn_dec[1:(nmin-1),j])
  ssn25[j]   = data.matrix(quant)[2]
  ssn_sea[j] = mean(ssn_dec[1:(nmin-1),j])
}

#---------------------------------------------------------------------------

# plot the 25% curve
if (plot_25 == 'y') {
  lines(time_dec,ssn25,col="red",type="l",lwd=4)
}

#---------------------------------------------------------------------------

# plot the SEA curve
if (plot_sea == 'y') {
  lines(time_dec,ssn_sea,col="green",type="l",lwd=4)
}

#---------------------------------------------------------------------------

if (fit_rise == 'y') {
fit_model = 2

min_t = 0.5
max_t = 2.75

if (fit_model == 1) {
# fit an exponential to each curve and print out the parameters
for (i in 1:nmin) {
  it_good = (time_dec > min_t) & (time_dec < max_t)
  myModel <- lm(log(ssn_dec[i,it_good]) ~ time_dec[it_good])
  #summary(model)
  model_y = exp(myModel$coefficients[1])*exp(myModel$coefficients[2]*time_dec[it_good])
  lines(time_dec[it_good],model_y,type="l",col="blue",lwd=2)
  print(paste("cycle/growth-rate/gradient:",i, myModel$coefficients[2], exp(myModel$coefficients[1])))
}
}

if (fit_model ==2 ) {
# fit a linear model to each curve and print out the parameters
for (i in 1:nmin) {
  it_good = (time_dec > min_t) & (time_dec < max_t)
  myModel <- lm(ssn_dec[i,it_good] ~ time_dec[it_good])
  #summary(model)
  model_y = myModel$coefficients[1] + myModel$coefficients[2]*time_dec[it_good]
  lines(time_dec[it_good],model_y,type="l",col="blue",lwd=2)
  print(paste("cycle/intercept/gradient:",i, myModel$coefficients[1], myModel$coefficients[2]))
}
}

}

#---------------------------------------------------------------------------

# plot 6 and 7

if (whPlot >= 6) {

# rank of each cycle for the first 2.75 years (modal value)
# btw, a more general code would make the cycle number be 
# a variable. But this code specifically analyses cycle 25
# so I'm hard-wiring that point in on purpose. 

  # need to recalcuate the rankings since we're not 
  # including cycle 25
  
rank_dec  = 0.0*ssn_dec[1:24,]
ssn_dec   = ssn_dec[1:24,]

for (j in 1:length(time_dec)) {
  rank_dec[,j] = rank(ssn_dec[,j])
}
  
cycle_rank = rep(0.0, 24)
for (i in 1:24) {
  cycle_rank[i] = mode(rank_dec[i,13:24])
}

ssn_spe = ssn_ms
tmax = 10.5
fit_rise = 'n'

ssn_sea = 0.0*time_dec
ssn_QSEA = 0.0*time_dec
col = c("red","blue","green","")
lwd = 1
  
plot(dyear-tmin[1],ssn_spe,xlim=c(0,tmax),ylim=c(0,300),type="l",xlab="Time (Years)",ylab="SSN",col=1) #,
mae_QSEA = rep(0.0,24)
mae_mean = rep(0.0,24)
for (i in 1:24) {
  lines(dyear-tmin[i],ssn_ms,col=i,type="l",lwd=lwd)
  ssn_dec[i,] = approx(dyear-tmin[i],ssn_spe,time_dec)$y
  for (j in 1:length(time_dec)) {
    quant = quantile(ssn_dec[1:(nmin-1),j],probs=as.double(cycle_rank[i])/24.)
    ssn_QSEA[j]   = data.matrix(quant)[1]
  }
  lines(time_dec,ssn_QSEA,col=i,type="l",lwd=lwd,lty=2)
  mae_QSEA[i] = mae(ssn_dec[i,],ssn_QSEA)
  mae_mean[i] = mae(ssn_dec[i,],ssn_sea)
  #print("cycle/MAE_QSEA/MAE_mean: ",i, mae_QSEA[i],mae_mean)
  colour = i
  adjustcolor(colour, alpha.f=0.5)
  polygon(c(dyear-tmin[i], rev(time_dec)), c(ssn_ms, rev(ssn_QSEA)),
            col = rgb(col2rgb(i)[1]/255.,col2rgb(i)[2]/255.,col2rgb(i)[3]/255.,0.25), lty = 0)
}  
  
legend(0, 300, legend=1:nmin,col=1:nmin, lty=rep(1,nmin), cex=0.65)
  
for (j in 1:length(time_dec)) {
  ssn_sea[j] = mean(ssn_dec[1:(nmin-1),j])
}

if (whPlot ==7) {
# plot 7 (requires plot 6 analysis)
# finally, compute the MAEs for the predicted curves
# as compared with the average curve
cycle1_24 = 1:24

# This plots the values of MAE and compares with the ratio...
# not very illuminating
par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
plot(cycle1_24,mae_QSEA,xlim=c(0,25),ylim=c(0,50),type="l",xlab="Cycle Number",ylab="MAE (ssn^2)",col="blue") #,
lines(cycle1_24,mae_mean,col="red")  
legend(0,10, legend=c("QSEA","Average"),col=c("blue","red"), lty=rep(1,2), cex=1)
par(new = TRUE)
plot(cycle1_24,mae_QSEA/mae_mean, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(mae_QSEA/mae_mean)))
mtext("MAE(QSEA)/MAE(AVERAGE)", side=4, line=3)

hist(mae_QSEA/mae_mean) # this is a simple histogram. OK. 

# This is the best representation of the data
den <- density(mae_QSEA/mae_mean)
plot(den, frame = FALSE, col = "blue",main = "MAE(QSEA) / MAE(AVERAGE)")
abline(v=1, col="red",lty=1)

}

}

}

