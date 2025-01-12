rm(list=ls())

# load functions 
source("Functions/WindowFunctions.R")

# specify species as orange tip
SpeciesName<-"Anthocharis_cardamines"
print(SpeciesName)

# load data
AbundanceData<-read.csv(paste0("Data/UKBMS/Site indices/ukbmsSiteIndices2021_",SpeciesName,".csv"),row.names = 1)

# specify windows to search
# this is specifying the period between (i) CountPercentPrev - in the previous year, 30 days before first observation for this species (across years)
# and (ii) CountPercentCurr - in the current year, the day when 99% of counts for this species had been obserevd (across years)
# these values are based on raw ukbms count data
CountPercentCurr<-99 
CountPercentPrev<- -30 
WindowDurations<-c(seq(20,150,10))
WindowStartInterval<-10

# generate folders for outputs if they aren't already there
if(!file.exists("SlidingWindows")) dir.create("SlidingWindows")
if(!file.exists(paste0("SlidingWindows/",SpeciesName))) dir.create(paste0("SlidingWindows/",SpeciesName))
if(!file.exists(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr))) dir.create(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr))
if(!file.exists(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/AR1/Space"))) dir.create(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/AR1/Space"))
if(!file.exists(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/AR1/Time"))) dir.create(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/AR1/Time"))
if(!file.exists(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/Time"))) dir.create(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/Time"))
if(!file.exists(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/Space"))) dir.create(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/Space"))
if(!file.exists(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/AR1/Time"))) dir.create(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/AR1/Time"))
if(!file.exists(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/AR1/Space"))) dir.create(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/AR1/Space"))

# Get climate values for windows with SpAv 1976-90
GetWindowsClimates(SpeciesName, AbundanceData, CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                   WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, SpAvSubset=1976:1990)

#############################

### Run window search with previous year abundance predictor as a yearly deviation
# This is the main window search model

# run temperature window search 
RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE)
# plot temperature window search 
PlotWindows(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
            CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
            WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
            AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE)

# run precipitation window search 
RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="Precip1",LastWindow="MaxTemp1", 
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE)
# plot precipitation window search 
PlotWindows(SpeciesName, Scaled=TRUE, FocWindow="Precip1",
            CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
            WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
            AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE)

# check models
CheckWindowModels(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                  CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                  WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                  AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE)
CheckWindowModels(SpeciesName, Scaled=TRUE, FocWindow="Precip1",
                  CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                  WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                  AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE)

################

### Run window search with previous year abundance predictor (not as a yearly deviation)

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)
PlotWindows(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
            CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
            WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
            AR1=TRUE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="Precip1",LastWindow="MaxTemp1", 
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)
PlotWindows(SpeciesName, Scaled=TRUE, FocWindow="Precip1",
            CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
            WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
            AR1=TRUE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)

CheckWindowModels(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                  CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                  WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                  AR1=TRUE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)
CheckWindowModels(SpeciesName, Scaled=TRUE, FocWindow="Precip1",
                  CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                  WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                  AR1=TRUE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)

################

### Run window search with no previous year abundance predictor 

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=FALSE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)
PlotWindows(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
            CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
            WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
            AR1=FALSE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="Precip1",LastWindow="MaxTemp1", 
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=FALSE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)
PlotWindows(SpeciesName, Scaled=TRUE, FocWindow="Precip1",
            CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
            WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
            AR1=FALSE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)

CheckWindowModels(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                  CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                  WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                  AR1=FALSE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)
CheckWindowModels(SpeciesName, Scaled=TRUE, FocWindow="Precip1",
                  CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                  WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                  AR1=FALSE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)

##############################

### Space only window search

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, SpaceOnly=TRUE)
PlotWindows(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
            CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
            WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
            AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, SpaceOnly=TRUE)

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="Precip1",LastWindow="MaxTemp1", 
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, SpaceOnly=TRUE)
PlotWindows(SpeciesName, Scaled=TRUE, FocWindow="Precip1",
            CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
            WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
            AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, SpaceOnly=TRUE)

CheckWindowModels(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                  CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                  WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                  AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, SpaceOnly=TRUE)
CheckWindowModels(SpeciesName, Scaled=TRUE, FocWindow="Precip1",
                  CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                  WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                  AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, SpaceOnly=TRUE)

##############################

### Time only window search

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, TimeOnly=TRUE)
PlotWindows(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
            CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
            WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
            AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, TimeOnly=TRUE)

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="Precip1",LastWindow="MaxTemp1", 
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, TimeOnly=TRUE)
PlotWindows(SpeciesName, Scaled=TRUE, FocWindow="Precip1",
            CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
            WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
            AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, TimeOnly=TRUE)

CheckWindowModels(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                  CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                  WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                  AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, TimeOnly=TRUE)
CheckWindowModels(SpeciesName, Scaled=TRUE, FocWindow="Precip1",
                  CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                  WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                  AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, TimeOnly=TRUE)

