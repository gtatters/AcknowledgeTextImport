############################## TABLE OF CONTENTS #######################################
#0.  Define all libraries and functions
#1.  Human Inputs - manually enter file name as variable f and file date
#2.  File handling - import file into dataframe called 'alldata'
#3.  Estimate time from file start or first time point and construct a Time Variable
#4.  Define the inputted variables, assign col.names, re-sample to 1 samp/sec
#5.  Lag time detection for WVD/O2/CO2 (Human Input here)
#6.  Align channels based on lag values into new dataframe called 'align.data'
#7.  Show the aligned data, compensated for the lag
#8.  Find the baseline regions (Human Input here)
#9.  Calculating the EWL, VH2O, VO2, VCO2
#10. Output finalised data to new file in variable f.out 

# READ ME BEFORE RUNNING!
# This script works with text files exported from Windaq's Acknowledge 1.8.1
# Importing Text files generated from Windaq Acknowledge
#    Two ways to generate txt files:
#    One without header (recommended).
#    Choose tab as the delimiter when exporting, since comma delimiter leads to excessive
#    coding
#    One with header, which creates a number of initial rows that need to be 
#    removed after import 
#    This script is not written to accomodate the header from Acknowledge files
# Baseline Information:   
#    This script is intended to work with a data set where the files starts on the baseline
#    I have not written the script to accomodate starting on the chamber
#    If you happen to have a file that does start on the chamber, you will have to edit the
#    text file in order to produce a few rows of cal data.  
# File info:
#    You will have to input directly your file name below. 
#    A new file will be created called 'filename_output.csv' which will contain the 
#    EWL and metabolic rate data with time stamps and all the necessary information for filtering
#    For completeness, please enter the date that the file was recorded.  Unfortunately,
#    Biopac does not save date information in the txt file.


#0. Define all libraries used and functions ####
library("Thermimage")
library("lattice")

rm(list = setdiff(ls(), lsf.str()))

#1. Human Inputs - Filename and date of actual experiment ####
mainDir<-"~/Dropbox/R/MyProjects/AcknowledgeImport/"
#mainDir<-"/Volumes/TattersallSD/ZebrafinchBrock/Analysis/TNZ/Metabolic/TXT/"
setwd(mainDir)
getwd() #show working directory

l.files<-list.files(mainDir,all.files=FALSE,full.names=FALSE,recursive=TRUE,
                    include.dirs=TRUE, pattern = "\\.txt$")
l.files
(f<-l.files[3])


filedate<-"2018-02-01"
f.out<-paste(gsub(".txt", "", f), "_output.csv", sep="")
f.pdf<-paste(gsub(".txt", "", f), "_plot.pdf", sep="")

#2. File Handling - manually enter file name as variable f ####
# define the file output names and locations
#  must give a filename ("in_file")

finfo = file.info(f)
#alldata<-read.csv(f, header=TRUE, skip=36) 
  # use this if you save with header and use comma as delimiter
#alldata<-read.delim(f, header=TRUE, sep="", skip=35) 
  # use this if you save with header and use tab as delimiter
  # skip 35 creates proper number of coolumns
# Recommend this method for robust input 
alldata<-read.delim(f, header=FALSE, sep="", skip=26)
  # Do not export a header
  # Save the timescale when exporting 
  # use tab as the delimter
rownames(alldata)<-NULL
str(alldata)
alldata<-na.omit(alldata)
  # remove NA values in alldata dataframe in case bad records are imported

# if it imports an extra column :
#alldata<-alldata[,-ncol(alldata)]
  # remove last column (ncol(alldata)) - error during export/import from WDQ leads to an 
  # extra column of non data
#str(alldata)
# if it imports a row with # of samples in file:
#alldata<-alldata[-1,]
  # remove first row (contains the total # of samples) from all data 
#str(alldata)

#3. Estimate time from file start or first time point ####
(timestart<-alldata[1,1]) 
  # this is how many seconds have elapsed since midnight
  # should be able to conver this to a real time value
hrs.decimal<-timestart/3600
hrs<-trunc(hrs.decimal)
mins<-trunc((hrs.decimal-hrs)*60)
secs<-round(((hrs.decimal-hrs)*60-mins)*60,digits=0)
start.time<-paste(filedate, " ", sprintf("%02d",hrs),":",sprintf("%02d",mins),":",
                  sprintf("%02d",secs), sep="")
  # pad values with a zero in front for single digit h,m,s
start.time<-as.POSIXct(start.time, format="%Y-%m-%d %H:%M:%OS")
start.time
  # populate vector variable times with the elapsed time channel from alldata
Fs<-mean(diff(alldata[,1]))  # sample rate in seconds
cat('Sample rate (in samples/s) is: ', Fs)
cat('\n')
TimeVar<-seq(start.time,start.time+Fs*(nrow(alldata)-1),Fs)
alldata<-data.frame(TimeVar,alldata)
rm(TimeVar)

#4. Define the Inputted Variables ####  

# Starting Units of Inputted Variables
# Time: Time variable in POSIXct units (real time in format="%Y-%m-%d %H:%M:%OS")
# ElapTime: Elapsed Time in seconds
# WVD: Water vapour Density in ug H2O/mL
# FO2:  Oxygen in fractional (in ACQ is saved as percent. Divide by 100 below to get F)
# FCO2: CO2 in fractional (in ACQ is saved as percent. Divide by 100 below to get F)
# Flow rate: mL/min (flow standardised to standard temperature and pressure)
# Inc Temp: Incubator temperature (attempted temperature)
# Resp Temp: Respirometer temperature (temperature around bird)
# Foxbox temp: Temperature of the area around the foxbox
# Pressure: differential pressure waveform (volts: arbitrary measurement)
# Barometric Pressure: kPa
# DewPt: Dew Point of the inspired air
# WVD, O2, CO2, FR, Ti, Tr
col.names<-c("Time", "ElapTime","WVD","FO2","FCO2","FR","Ti","Tr","Pres", "Tfox", "Baseline", "BP","DewP")
colnames(alldata)<-col.names
rownames(alldata)<-NULL

alldata$"ElapTime"<-alldata$"ElapTime"-alldata$"ElapTime"[1]
alldata$"FO2"<-alldata$"FO2"/100
alldata$"FCO2"<-alldata$"FCO2"/100
alldata$"Baseline"[which(alldata$"Baseline"<=0.2)]<-0
alldata$"Baseline"[which(alldata$"Baseline">0.2)]<-1

# Elapsed time now starts at 0 seconds, gas percentages now expressed as fractions

# Keeping this for posterity in case the alldata variable was sampled at a different sample rate
# this is commented out from when we had recorded at > 1 samp/s
#alldata1s<-apply(alldata,2,meanEveryN,50,2,showsamples=FALSE)
  # resample all data from native sample rate down by 50
  # this will convert the time intervals to be ~ every second
alldata1s<-alldata
# populate alldata1s with alldata and remove alldata
# only doing this to retain easier ability to re-write script in case files are sampled differently
rm(alldata)
alldata1s<-data.frame(alldata1s,row.names=NULL)
colnames(alldata1s)<-col.names

#5. Lag detection for WVD/O2/CO2 & optimisation (Human Input) ####
# since the gas channels are not necessarily aligned, this section will attempt to find the 
# time lag which corresponds to the highest correlation between each gas channel and the 
# Baseline channel (~0 volts = Baseline, ~5 volts = Chamber)
maxlag<-60  ##### 
  # 60 seconds is deemed to be the max # of seconds required to correct a lag
len<-nrow(alldata1s)
  # Order of gas flow measurements: WVD-CO2-O2
maxlen<-nrow(alldata1s)-maxlag+1
#w1<-alldata1s$"WVD"[1:maxlen]-5
baselinedata<-alldata1s$"Baseline"[1:maxlen]

corelw1<-NULL;corelc1<-NULL;corelo1<-NULL
par(mfrow=c(3,1),mar=c(3,3,3,3))
for (i in maxlag:1)  
{
  j<-i+len-maxlag
  lmw1<-lm(alldata1s$"WVD"[i:j]~baselinedata)
  lmc1<-lm(alldata1s$"FCO2"[i:j]~baselinedata)
  lmo1<-lm(alldata1s$"FO2"[i:j]~baselinedata)
  
  rsw1<-summary(lmw1)$r.squared
  rsc1<-summary(lmc1)$r.squared
  rso1<-summary(lmo1)$r.squared
  # 
  # plot(baselinedata,alldata1s$"WVD"[i:j], type="l", lty=3, ylab="", bty="n",
  #     main=paste("WVD Lag",i,"r²=",round(rsw1,digits=3)))
  # abline(lmw1$coef[1],lmw1$coef[2], col="red", lwd=5)
  # plot(baselinedata,alldata1s$"FCO2"[i:j], type="l", lty=3, ylab="",  bty="n",
  #      main=paste("FCO2 Lag",i,"r²=",round(rsc1,digits=3)))
  # abline(lmc1$coef[1],lmc1$coef[2], col="red", lwd=5)
  # plot(baselinedata,alldata1s$"FO2"[i:j], type="l", lty=3, ylab="", bty="n",
  #      main=paste("FO2 Lag",i,"r²=",round(rso1,digits=3)))
  # abline(lmo1$coef[1],lmo1$coef[2], col="red", lwd=5)
  # 
  corelw1<-c(corelw1,rsw1)
  corelc1<-c(corelc1,rsc1)
  corelo1<-c(corelo1,rso1)
  
  #Sys.sleep(0.1)
}

# Lag estimates in seconds for wvd and co2
w1.lag<-which(rev(corelw1)==max(corelw1, na.rm=TRUE))  
c1.lag<-which(rev(corelc1)==max(corelc1, na.rm=TRUE))  
o1.lag<-which(rev(corelo1)==max(corelo1, na.rm=TRUE))


#6. Align Channels based on lag ######
# Align Channels considering lags
# First need to remember that the oxygen channels are all laggin behind the real time channels
# This was determined in the lab by measuring how long it took the oxygen meter to react
# a series of NA values will have to be added to the end of any channel with a lag
# w1.lag<-8 # redundant now we have a Baseline channel.
# store the length of the dataset prior to lag correction in variable len
len<-nrow(alldata1s)

# Define start and end values for the lagged channels
  #w1.lag.end<-len-w1.lag
  #c1.lag<-c1.lag+w1.lag
  #o1.lag<-o1.lag+w1.lag
# populate temporary variables o1, c1 and w1 with lagged values
  w1<-c(alldata1s$"WVD"[w1.lag:len],rep(NA,w1.lag-1))
  c1<-c(alldata1s$"FCO2"[c1.lag:len],rep(NA,c1.lag-1))
  o1<-c(alldata1s$"FO2"[o1.lag:len],rep(NA,o1.lag-1))

  align.data<-data.frame(alldata1s$"Time", alldata1s$"ElapTime",w1,o1,c1,alldata1s$"FR",alldata1s$"Ti",
                       alldata1s$"Tr",alldata1s$"Pres",alldata1s$"Tfox", alldata1s$"Baseline",
                       alldata1s$"BP", alldata1s$"DewP")
  colnames(align.data)<-col.names
  #align.data<-na.omit(align.data)
  align.data<-align.data[-c(seq(nrow(align.data)-o1.lag,nrow(align.data),1)),]
  # remove the Na values at the end that resulted from the lag correction

#7. Show the aligned data ####
par(mfrow=c(2,1), mar= c(2, 2, 2, 2))
attach(alldata1s)
plot(ElapTime, Baseline, type="l", lty=3, xlab="", ylab="", xaxt="n", yaxt="n",
     col="black", axes=FALSE)
par(new = TRUE)
plot(ElapTime, WVD, type="l", xlab="", ylab="", xaxt="n", yaxt="n",
     col="blue", axes=FALSE, main="Plotted without lag")
par(new = TRUE)
plot(ElapTime, FCO2, type="l", xlab="", ylab="", xaxt="n", yaxt="n",
     col="green", axes=FALSE)
par(new = TRUE)
plot(ElapTime, FO2, type="l", xlab="", ylab="", xaxt="n", yaxt="n", 
     col="red", axes=FALSE)
detach(alldata1s)

cat(paste("WVD Lag",w1.lag,"samples"))
cat('\n')
cat(paste("FCO2 Lag",c1.lag,"samples"))
cat('\n')
cat(paste("FO2 Lag", o1.lag,"samples"))
cat('\n')
par(mar= c(2, 2, 2, 2))

attach(align.data)
plot(ElapTime, Baseline, type="l", lty=3, xlab="", ylab="", xaxt="n", yaxt="n",
     col="black", axes=FALSE)
par(new = TRUE)
plot(ElapTime, WVD, type="l", xlab="", ylab="", xaxt="n", yaxt="n",
     col="blue", axes=FALSE, main=paste("WVD Lag",w1.lag,"samples"))
par(new = TRUE)
plot(ElapTime, FCO2, type="l", xlab="", ylab="", xaxt="n", yaxt="n",
     col="green", axes=FALSE)
par(new = TRUE)
plot(ElapTime, FO2, type="l", xlab="", ylab="", xaxt="n", yaxt="n", 
     col="red", axes=FALSE)
detach(align.data)


#8. Baseline determination & delta determination ####
# Assumptions: 
# Baseline variable has been changed so that Baseline = 0 means baseline
# Baseline = 1 means reading from chamber
attach(align.data)
# Define first the starting points for the baselines and chambers switches
diff.base<-diff(Baseline)
diff.base<-c(diff.base[1], diff(Baseline))
# construct a variable the same length as Baseline which is the difference in
# adjacent points in time.  
# diff.base will help to locate where the transition points are between 0 and 1
# 1-0 = 1, meaning chamber start point
# 0-1 = -1, meaning baseline start point
if(Baseline[1]==1) 
  {
  chamber.starts<-c(1, which(diff.base==1))
  base.starts<-which(diff.base==-1)
}
if(Baseline[1]==0)
{
  chamber.starts<-which(diff.base==1)
  base.starts<-c(1,which(diff.base==-1))
}
# Define baseline points which will be trustworthy for calibration
# Assumptions: 
#     The respirometry system starts with a baseline measurement?  
#     There will be some delay following the baseline switch before data can be trusted   
if(length(chamber.starts)>length(base.starts)) base.starts<-c(base.starts,length(Baseline))
if(length(base.starts)>length(chamber.starts)) chamber.starts<-c(chamber.starts,length(Baseline))
chamber.delay<-40   #####
# delay in samples (use this to ignore periods of time after switching to chamber)
if(base.starts[1]<chamber.starts[1])
{
  calib.points<-Map(':', base.starts+chamber.delay, chamber.starts)
  # create a list that consists of the index values corresponding to each baseline period
  chamber.points<-Map(':', chamber.starts[1:(length(chamber.starts)-1)]+chamber.delay, 
                    base.starts[2:(length(chamber.starts))]-chamber.delay)
  
}
chamber.starts.n1<-c(chamber.starts[2:length(chamber.starts)],nrow(align.data))
if(base.starts[1]>chamber.starts[1])
{
  calib.points<-Map(':', base.starts+chamber.delay, chamber.starts.n1)
  # create a list that consists of the index values corresponding to each baseline period
  chamber.points<-Map(':', chamber.starts[1:(length(chamber.starts))]+chamber.delay, 
                      base.starts[1:(length(chamber.starts))]-chamber.delay)
}
index<-seq(1, nrow(align.data), 1)
chamber.index<-unlist(chamber.points)
# use chamber.index to filer out non-measurement periods
index.na<-rep(NA, nrow(align.data))
non.chamber.index<-index[-chamber.index]
index.na[chamber.index]<-1
# use index.na when saving the final file.  1 = periods measurement periods

# the following section will define the regions of the data where baselines occur
# calib, refers to the value within a baseline period that best corresponds to the value 
# to which the baseline would be assigned/calibrated to.
calib.WVD<-vector("list", length(calib.points))
baseline.WVD<-vector("list", length(calib.points))
calib.FO2<-vector("list", length(calib.points))
baseline.FO2<-vector("list", length(calib.points))
calib.FCO2<-vector("list", length(calib.points))
baseline.FCO2<-vector("list", length(calib.points))
# populate .na variable with NA to the same length as the respirometry data
# these .na variables will be built up with actual baseline data.  The reason to populate with
# NA is in case an error is made, R can handle NAs
WVD.na<-rep(NA, nrow(align.data))
FO2.na<-rep(NA, nrow(align.data))
FCO2.na<-rep(NA, nrow(align.data))
par(mfrow=c(2,round(length(calib.points)/2)), mar=c(4,4,4,2))
for (i in 1:length(calib.points))
{
  j<-i+1
  x<-unlist(calib.points[i])
  calib.WVD[i]<-min(WVD[x])
  calib.FO2[i]<-quantile(FO2[x], probs=seq(0,1,0.1), na.rm=TRUE)[9]
  #calib.FCO2[i]<-quantile(FCO2[x], probs=seq(0,1,0.1), na.rm=TRUE)[2]
  calib.FCO2[i]<-median(FCO2[x], na.rm=TRUE)
  #calib.FCO2[i]<-mean(FCO2[which(FCO2[x]<quantile(FCO2[x], probs=seq(0,1,0.1), na.rm=TRUE)[2])])
  # define the stretches of baselines in a list variable
  baseline.WVD[[i]]<-rep(unlist(calib.WVD[i]), length(x))
  baseline.FO2[[i]]<-rep(unlist(calib.FO2[i]), length(x))
  baseline.FCO2[[i]]<-rep(unlist(calib.FCO2[i]), length(x)) 
  plot(ElapTime[x], FO2[x], type="p", main=paste("Baseline # ", i))
  title(main=paste("FO2 =", calib.FO2[i]), line=0.5, cex=0.8)
  lines(c(ElapTime[min(x)], ElapTime[max(x)]), c(calib.FO2[i], calib.FO2[i]), col="red", lwd=4)
  plot(ElapTime[x], FCO2[x], type="p", main=paste("Baseline # ", i))
  title(main=paste("FCO2 =", calib.FCO2[i]), line=0.5, cex=0.8)
  lines(c(ElapTime[min(x)], ElapTime[max(x)]), c(calib.FCO2[i], calib.FCO2[i]), col="green", lwd=4)
}
#cat('These are your Baselines')
#cat('\n')
readline("Press enter to continue")

for (i in 1:(length(calib.points)-1))
{
  j<-i+1
  x<-unlist(calib.points[i])
  # replace the wvd.na with interpolation
  between.cal<-c(rev(unlist(calib.points[i]))[1], base.starts[j])
  between.WVD<-c(unlist(calib.WVD[i]), unlist(calib.WVD[j]))
  between.FO2<-c(unlist(calib.FO2[i]), unlist(calib.FO2[j]))
  between.FCO2<-c(unlist(calib.FCO2[i]), unlist(calib.FCO2[j]))
  between.cal.points<-seq(between.cal[1]+1, between.cal[2]-1, 1)
  lmbetweencal<-lm(between.WVD~ElapTime[between.cal])
  WVD.na[x]<-unlist(baseline.WVD[i])
  WVD.na[between.cal.points]<-lmbetweencal$coef[1]+lmbetweencal$coef[2]*ElapTime[between.cal.points]
  lmbetweencal<-lm(between.FO2~ElapTime[between.cal])
  FO2.na[x]<-unlist(baseline.FO2[i])
  FO2.na[between.cal.points]<-lmbetweencal$coef[1]+lmbetweencal$coef[2]*ElapTime[between.cal.points]
  lmbetweencal<-lm(between.FCO2~ElapTime[between.cal])
  FCO2.na[x]<-unlist(baseline.FCO2[i])
  FCO2.na[between.cal.points]<-lmbetweencal$coef[1]+lmbetweencal$coef[2]*ElapTime[between.cal.points]
}
i<-length(calib.points)
x<-unlist(calib.points[i])
WVD.na[x]<-unlist(baseline.WVD[i])
baseline.WVD<-WVD.na
FO2.na[x]<-unlist(baseline.FO2[i])
baseline.FO2<-FO2.na
FCO2.na[x]<-unlist(baseline.FCO2[i])
baseline.FCO2<-FCO2.na
# replace baseline with wvd.na, which contains the interpolated 
# and non-interpolated sections

# Show the gas data along with interpolated baselines (red lines)
par(mfrow=c(3,1), mar= c(2, 2, 2, 2))
plot(ElapTime, WVD, type="l")
lines(ElapTime, baseline.WVD, lwd=3, col="red")
plot(ElapTime, FO2, type="l")
lines(ElapTime, baseline.FO2, lwd=3, col="red")
plot(ElapTime, FCO2, type="l")
lines(ElapTime, baseline.FCO2, lwd=3, col="red")

# Calculate Water Vapour Pressures used in FiH2O and FeH2O calculations
WVDi<-baseline.WVD
delta.WVD<-WVD-WVDi

# for those measurements where the DewP is recording a "zero" value since the machine is not 
# connected, generate a DewP variable that is based on the assumption of constant incoming DewP
if(mean(DewP)<0.5) # 0.5 is arbitrary but less than we can reasonably generate
{
  DP<-16.116*log(mean(unlist(calib.WVD)))-25.573
  # equation empirically determined from 
  # http://hyperphysics.phy-astr.gsu.edu/hbase/kinetic/watvap.html#c1
  DewP[1:nrow(align.data)]<-DP
}

WVPi<-(611*10^(7.5*DewP/(237.7+DewP)))/1000 
  # converts dew point into WVP in kPa
  # from Handbook of chemistry and physics
WVPi.t<-WVPi*(26+273)/(DewP+273)
  # convert to water vapour pressure at RH-300 temperature (~26C)
WVPe<-0.1421*WVD-0.1038
  # empirically determined conversion from WVD(ug/mL) to WVP (kPa)
  # from Handbook of chemistry and physics table
DewPte<-15.799*log(WVPe)+6.9848
  # Dewpoint of excurrent air, required for temperature compensation
  # of FeH2O value
WVPe.t<-WVPe*(26+273)/(DewPte+273)
  # temperature corrected WVPe to the temperature of the RH-300 meter

# Finalise FH2O values
FiH2O<-WVPi.t/BP
FeH2O<-WVPe.t/BP
  # define Fi and Fe H2O, required for final respirometry calculations

# Finalise FCO2 values
FiCO2<-baseline.FCO2
delta.FCO2<-FCO2-FiCO2
FeCO2<-FiCO2+delta.FCO2

# Finalise FO2 values
FO2.span<-0.2095 #####  
# Change this if span calibrated to another value
FO2.real<-(1-mean(FiH2O, na.rm=TRUE)-mean(FiCO2, na.rm=TRUE))*0.2095
  # this is the true level the oxygen analyser should be (0.2095 because we're using CO2 free air) 
  # to correct our FO2 measurements, multiply oxygen values by FO2.real/FO2.span
FiO2<-(1-FiH2O-FiCO2)*0.2095
  # multiplied by 0.2095 because the incurrent oxygen is known from our dewpoint measurement
  # and knowledge of incurrent CO2 levels
delta.FO2<-(baseline.FO2-FO2)*FO2.real/FO2.span
FeO2<-FiO2-delta.FO2
detach(align.data)


# 9.  Calcuating EWL, VH2O, VO2, VCO2  ####
# assuming channels are aligned, one can calculate these values
# EWL units: mg H2O/h
attach(align.data)

EWL<-delta.WVD*FR/1000*60
# ug/mL * mL/min = ug/min * 1/1000 mg/ug * 60 min/h = mg H2O/h
EWL[non.chamber.index]<-NA

FiN2<-1-FiO2-FiCO2-FiH2O
FeN2<-1-FeO2-FeCO2-FeH2O

VH2O<-60*FR*((FeH2O*(FiN2)/(FeN2))-FiH2O)
# units: mL/min * 60 min/h = mL H2O / h
VH2O[non.chamber.index]<-NA

VCO2<-60*FR*((FeCO2*(FiN2)/(FeN2))-FiCO2)
# units: mL/min * 60 min/h = mL O2 / h
VCO2[non.chamber.index]<-NA

VO2<-60*FR*(FiO2-(FeO2*(FiN2)/(FeN2)))
# units: mL/min * 60 min/h = mL CO2 / h
VO2[non.chamber.index]<-NA

RQ<-VCO2/VO2
# units: none
RQ[non.chamber.index]<-NA

pdf(file=f.pdf, width=10, height=10, family="Helvetica")
par(mfrow=c(4,1), mar=c(3,5,1,1))
plot(ElapTime[chamber.index], EWL[chamber.index], col="blue", bty="n",
     ylab=expression("EWL (mg " ~ H[2]~O ~ "/h)"))
plot(ElapTime[chamber.index], VO2[chamber.index], col="red",  bty="n",
     ylab=expression(VO[2] ~ "(mL " ~ O[2] ~ "/h)"))
plot(ElapTime[chamber.index], VCO2[chamber.index], col="green", bty="n",
     ylab=expression(VCO[2] ~ "(mL " ~ CO[2] ~ "/h)"))
plot(ElapTime[chamber.index], RQ[chamber.index], col="grey",  bty="n", ylab="RQ",
     ylim=c(0.6,1.1))
detach(align.data)
dev.off()

#10. Output data into new data frame ####
# Note: output file will have NA values for periods of time that correspond to baseline
# or non-chamber time points

output<-data.frame(align.data, index.na, EWL, VO2, VCO2, RQ)
col.names<-c("Time", "ElapTime","WVD","FO2","FCO2","FR","Ti","Tr",
             "Pres", "Tfox", "Baseline", "BP","DewP", "Index","EWL","VO2","VCO2","RQ")
colnames(output)<-col.names
write.csv(output, f.out)
