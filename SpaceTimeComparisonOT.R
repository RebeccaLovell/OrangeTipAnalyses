
rm(list=ls())

# load functions 
source("Functions/SpaceTimeFunction.R")

# specify species as orange tip
SpeciesName<-"Anthocharis_cardamines"
print(SpeciesName)

# specify windows to search
# this is specifying the period between (i) CountPercentPrev - in the previous year, 30 days before first observation for this species (across years)
# and (ii) CountPercentCurr - in the current year, the day when 99% of counts for this species had been obserevd (across years)
# these values are based on raw ukbms count data
CountPercentCurr<-99
CountPercentPrev<--30
WindowDurations<-c(seq(20,150,10))
WindowStartInterval<-10


# generate folders for outputs if there isn't already one
if(!file.exists("SpaceTimeComparison")) dir.create("SpaceTimeComparison")
if(!file.exists(paste0("SpaceTimeComparison/",SpeciesName))) dir.create(paste0("SpaceTimeComparison/",SpeciesName))
if(!file.exists(paste0("SpaceTimeComparison/",SpeciesName,"/AR1/"))) dir.create(paste0("SpaceTimeComparison/",SpeciesName,"/AR1"))
if(!file.exists(paste0("SpaceTimeComparison/",SpeciesName,"/AR1Hurdle/"))) dir.create(paste0("SpaceTimeComparison/",SpeciesName,"/AR1Hurdle"))
if(!file.exists(paste0("SpaceTimeComparison/",SpeciesName,"/Hurdle/"))) dir.create(paste0("SpaceTimeComparison/",SpeciesName,"/Hurdle"))


#############################

### Run space versus time model with the previous year abundance predictor as a yearly deviation (AR1yeardev)
# this is the main model

# run space time comparison model
SpaceTimeModel(SpeciesName, AR1=TRUE, AR1YearDev=TRUE, Hurdle=TRUE,
               CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
               WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
               nitt=5000000, thin=100, SpAvSubset=1976:1990)

# find slope differences and plot
CompareSpaceTime(SpeciesName,  AR1=TRUE, AR1YearDev=TRUE, Hurdle=TRUE,
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval,
                 SpAvSubset=1976:1990, nitt=5000000)

# do model checking 
CheckSpaceTimeModel(SpeciesName,  AR1=TRUE, AR1YearDev=TRUE, Hurdle=TRUE,
                    CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                    WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval,
                    SpAvSubset=1976:1990, nitt=5000000)

#############################

### Run space versus time model with the previous year abundance predictor (not as a yearly deviation, AR1)

SpaceTimeModel(SpeciesName, AR1=TRUE, AR1YearDev=FALSE, Hurdle=TRUE,
               CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
               WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
               nitt=5000000, thin=100, SpAvSubset=1976:1990)

CompareSpaceTime(SpeciesName,  AR1=TRUE, AR1YearDev=FALSE, Hurdle=TRUE,
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval,
                 SpAvSubset=1976:1990, nitt=5000000)

CheckSpaceTimeModel(SpeciesName,  AR1=TRUE, AR1YearDev=FALSE, Hurdle=TRUE,
                    CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                    WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval,
                    SpAvSubset=1976:1990, nitt=5000000)


#############################

### Run space versus time model with no previous year abundance predictor

SpaceTimeModel(SpeciesName, AR1=FALSE, AR1YearDev=FALSE, Hurdle=TRUE,
               CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
               WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
               nitt=5000000, thin=100, SpAvSubset=1976:1990)

CompareSpaceTime(SpeciesName,  AR1=FALSE, AR1YearDev=FALSE, Hurdle=TRUE,
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval,
                 SpAvSubset=1976:1990, nitt=5000000)

CheckSpaceTimeModel(SpeciesName,  AR1=FALSE, AR1YearDev=FALSE, Hurdle=TRUE,
                    CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                    WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval,
                    SpAvSubset=1976:1990, nitt=5000000)

##############################

### Run non-hurdle space versus time model with the previous year abundance predictor as a yearly deviation (AR1yeardev)


SpaceTimeModel(SpeciesName, AR1=TRUE, AR1YearDev=TRUE, Hurdle=FALSE,
               CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
               WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
               nitt=2500000, thin=40, SpAvSubset=1976:1990)

CompareSpaceTime(SpeciesName,  AR1=TRUE, AR1YearDev=TRUE, Hurdle=FALSE,
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval,
                 SpAvSubset=1976:1990, nitt=2500000)

CheckSpaceTimeModel(SpeciesName,  AR1=TRUE, AR1YearDev=TRUE, Hurdle=FALSE,
                    CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                    WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval,
                    SpAvSubset=1976:1990, nitt=2500000)


