
library(ncdf4) # met office data is .nc format
library(lubridate)

# The data used is publicly available from the UKBMS (butterfly data) and the HadUK grid (climate data)
# UKBMS site index and site location data (note that the publicly available site location data excludes sensitive sites) is available here https://catalogue.ceh.ac.uk/documents/571a676f-6c32-489b-b7ec-18dcc617a9f1
# HadUK grid daily temperature and precipitation data is available here https://catalogue.ceda.ac.uk/uuid/4dc8450d889a491ebb20e724debe2dfb 

#############################################
### Match UKBMS sites to met office grid ####
#############################################

# butterfly and climate data are both on OS grid

### UKBMS site data ####

# load UKBMS site data 
SiteData<-read.table("UKBMS/site data.txt", header = TRUE, sep = ",")

# rename columns 
colnames(SiteData)[which(colnames(SiteData)=="Site.Number")]<-"Site" 
colnames(SiteData)[which(colnames(SiteData)=="Easting")]<-"SiteEasting" 
colnames(SiteData)[which(colnames(SiteData)=="Northing")]<-"SiteNorthing" 

# For each transect start point we have an ordnance survey grid reference easting and northing

### Climate data ####

# each year's climate data is a separate file 
# load one weather data file to get the grid coordinates from 
ncdata<-nc_open("tasmax_hadukgrid_uk_1km_day_20211201-20211231.nc")

# attributes(ncdata$var) 

# Get projection coordinates - corners of grid cells
# projection_x_coordinate_bnds and projection_y_coordinate_bnds are the coordinates of grid cell boundaries (transverse mercator projection)
projection_x_coordinate_bnds <- ncvar_get(ncdata, "projection_x_coordinate_bnds")
projection_y_coordinate_bnds <- ncvar_get(ncdata, "projection_y_coordinate_bnds")

### Match UKBMS site locations to climate data cells ####

# add columns to store cell info for each site
SiteData$CellNorthingMax<-SiteData$CellNorthingMin<-NA
SiteData$CellEastingMax<-SiteData$CellEastingMin<-NA

# get vectors of the max and min coordinates for all cells
xmaxs<-projection_x_coordinate_bnds[2,]
xmins<-projection_x_coordinate_bnds[1,]
ymaxs<-projection_y_coordinate_bnds[2,]
ymins<-projection_y_coordinate_bnds[1,]

# find which grid cell, the UKBMS sites fall within 

for(focrow in 1:nrow(SiteData)){ # for each site
  
  # get its easting and northing 
  sitex<-SiteData$SiteEasting[focrow]
  sitey<-SiteData$SiteNorthing[focrow]

  # find the cell which the site falls within  
  cellx_position<-which(sitex<xmaxs & sitex>=xmins)
  celly_position<-which(sitey<ymaxs & sitey>=ymins)
  
  # store the coordinates for the corners of the cell
  SiteData$CellNorthingMax[focrow]<-projection_y_coordinate_bnds[2,celly_position]
  SiteData$CellNorthingMin[focrow]<-projection_y_coordinate_bnds[1,celly_position]
  SiteData$CellEastingMax[focrow]<-projection_x_coordinate_bnds[2,cellx_position]
  SiteData$CellEastingMin[focrow]<-projection_x_coordinate_bnds[1,cellx_position]
  
  # and the position in projection_x/y_coordinate_bnds
  SiteData$celly_index[focrow]<-celly_position
  SiteData$cellx_index[focrow]<-cellx_position
  # this indicates the columns of projection_x_coordinate_bnds and projection_y_coordinate_bnds
  # corresponding to the cell and can be used later to get the cell's temperature data 
  # this is projection_y_coordinate and projection_x_coordinate
}

# only retain columns we want 
SiteData<-SiteData[,match(c("Site","Gridreference","SiteEasting", "SiteNorthing",
                            "Length", "First.year.surveyed","Last.year.surveyed",
                            "Sensitive..1..no.0...yes.", "CellNorthingMin", "CellNorthingMax", 
                            "CellEastingMin", "CellEastingMax", "cellx_index", "celly_index"),colnames(SiteData))] 


### Add a 50 km grid to partially account for spatial autocorrelation in the models ###

# Easting and northing are in m. To make a grid, divide site easting/northing by desired grid cell size (in m).  
# use floor to get the lower x or y value for the cell that each site is in (ie bottom left coordinate)
# paste these together to give a grid cell ID

# a 50 km x 50 km grid is 50,000 m
SiteData$GridCell50km<-as.factor(paste0(floor(SiteData$SiteEasting/50000),"_", floor(SiteData$SiteNorthing/50000)))

write.csv(SiteData, file="UKBMS/SiteData_ClimateGrid.csv")


##############################
### Get site climate data ####
##############################

# SiteData<-read.csv("UKBMS/SiteData_ClimateGrid.csv")

# identify sites and years of interest 
Sites<-SiteData$Site
Years<-1975:2021
SiteYears<-paste0(expand.grid(Sites,Years)[,1],"_",expand.grid(Sites,Years)[,2])
nsiteyears<-length(SiteYears)

### get daily max temperature values for each siteyear ###

# make data frame with a row for each siteyear and a column for each day's climate
MaxTempData<-as.data.frame(matrix(nrow=nsiteyears,ncol=(366+3)))
colnames(MaxTempData)<-c("SiteYear","Site","Year",paste0("Day",1:366))

MaxTempData$SiteYear<-SiteYears # input siteyears

# get site and years for each row
MaxTempData$Site<-sapply(strsplit(MaxTempData$SiteYear,split="_"),function(x) x[1])
MaxTempData$Year<-sapply(strsplit(MaxTempData$SiteYear,split="_"),function(x) x[2])

for(focyear in unique(MaxTempData$Year)){ # for each year
  
  yearrows<-which(MaxTempData$Year==focyear)
  Sites<-MaxTempData$Site[yearrows] # get the sites surveyed in this year 
  
  # the data is stored in monthly files so need to loop over months 
  for(focmonth in 1:12){
    
    print(c(focyear,focmonth))
    
    # find the ordinal day number of the first day of the month
    # this gives 1st Jan as day 0 so added one to get ordinal day 
    MonthStartDay<-as.POSIXlt(paste(1,focmonth,focyear),format="%d %m %Y")$yday+1 
    
    # find the number of days in this month 
    ndays<-as.numeric(days_in_month(as.POSIXlt(paste(1,focmonth,focyear),format="%d %m %Y"))) 
    
    # load this month's data
    ncdata<-nc_open(paste0("tasmax_hadukgrid_uk_1km_day_",focyear,sprintf("%02d",focmonth),"01-",focyear,sprintf("%02d",focmonth),ndays,".nc"))
    
    # extract max temp             
    maxtemp <- ncvar_get(ncdata, "tasmax") 
    # double tasmax[projection_x_coordinate,projection_y_coordinate,time],where time is day
    # celly_index/cellx_index in AbundanceOutput indicates the columns of projection_x_coordinate_bnds 
    # and projection_y_coordinate_bnds corresponding to the cell 
    # this is projection_y_coordinate and projection_x_coordinate                
    
    for(focsite in Sites){ # for each site in this year, get the data for this month
      
      focrow<-intersect(which(MaxTempData$Site==focsite),yearrows)    
      
      for(day in 1:ndays){
        
        # find the day number for this day
        # minus 1 as we're adding day of the month to the number for first day of the month
        focdayno<-MonthStartDay+day-1  
        
        # Find the column corresponding to this day
        foccol<-which(colnames(MaxTempData)==paste0("Day",focdayno)) 
        
        SiteDataRow<-which(SiteData$Site==focsite)
        
        MaxTempData[focrow,foccol]<-maxtemp[SiteData$cellx_index[SiteDataRow],
                                            SiteData$celly_index[SiteDataRow],day]  
        
      }
    }
  }
}

write.csv(MaxTempData,file="Climate/Temperature/SiteYearMaxTemp.csv")

### get daily rainfall values for each siteyear ###

# make data frame with a row for each siteyear and a column for each day's climate
RainfallData<-as.data.frame(matrix(nrow=nsiteyears,ncol=(366+3)))
colnames(RainfallData)<-c("SiteYear","Site","Year",paste0("Day",1:366))

RainfallData$SiteYear<-SiteYears # input siteyears

# get site and years for each row
RainfallData$Site<-sapply(strsplit(RainfallData$SiteYear,split="_"),function(x) x[1])
RainfallData$Year<-sapply(strsplit(RainfallData$SiteYear,split="_"),function(x) x[2])

for(focyear in unique(RainfallData$Year)){ # for each year 
  
  yearrows<-which(RainfallData$Year==focyear)
  Sites<-RainfallData$Site[yearrows] # get the sites surveyed in this year 
  
  # the data is stored in monthly files so need to loop over months 
  for(focmonth in 1:12){
    
    print(c(focyear,focmonth))
    
    # find the ordinal day number of the first day of the month
    # this gives 1st Jan as day 0 so added one to get ordinal day 
    MonthStartDay<-as.POSIXlt(paste(1,focmonth,focyear),format="%d %m %Y")$yday+1 
    
    # find the number of days in this month 
    ndays<-as.numeric(days_in_month(as.POSIXlt(paste(1,focmonth,focyear),format="%d %m %Y"))) 
    
    # load this month's data
    ncdata<-nc_open(paste0("rainfall_hadukgrid_uk_1km_day_",focyear,sprintf("%02d",focmonth),"01-",focyear,sprintf("%02d",focmonth),ndays,".nc"))
    # extract rainfall             
    rainfall <- ncvar_get(ncdata, "rainfall") 
    # double rainfall[projection_x_coordinate,projection_y_coordinate,time],where time is day
    # celly_index/cellx_index in AbundanceOutput indicates the columns of projection_x_coordinate_bnds 
    # and projection_y_coordinate_bnds corresponding to the cell 
    # this is projection_y_coordinate and projection_x_coordinate                
    
    for(focsite in Sites){ # for each site in this year, get the data for this month
      
      focrow<-intersect(which(RainfallData$Site==focsite),yearrows)    
      
      for(day in 1:ndays){
        
        # find the day number for this day
        # minus 1 as we're adding day of the month to the number for first day of the month
        focdayno<-MonthStartDay+day-1  
        
        # Find the column corresponding to this day
        foccol<-which(colnames(RainfallData)==paste0("Day",focdayno)) 
        
        SiteDataRow<-which(SiteData$Site==focsite)
        
        RainfallData[focrow,foccol]<-rainfall[SiteData$cellx_index[SiteDataRow],
                                              SiteData$celly_index[SiteDataRow],day]  
        
      }
    }
  }
}

write.csv(RainfallData,file="Climate/Precipitation/SiteYearPrecip.csv")


###################################
### Sort UKBMS Site Index Data ####
###################################

# load site index data
SiteIndices<-read.csv("UKBMS/Site indices/ukbmsSiteIndices2021.csv")

# replace spaces with _ in species names
SiteIndices$SPECIES<-gsub(" ","_",SiteIndices$SPECIES)

# remove those before 1976 (UKBMS standardised recording period)
SiteIndices<-SiteIndices[-which(SiteIndices$YEAR<1976),] 

# add SiteYear Column 
SiteIndices$SiteYear<-paste0(SiteIndices$SITE.CODE,"_",SiteIndices$YEAR)

# remove duplicates
SiteIndices<-unique(SiteIndices) 

# Change site and year column names for consistency
colnames(SiteIndices)[which(colnames(SiteIndices)=="SITE.CODE")]<-"Site"
colnames(SiteIndices)[which(colnames(SiteIndices)=="YEAR")]<-"Year"

# load the site data/climate grid and merge this with counts
SiteData<-read.csv("UKBMS/SiteData_ClimateGrid.csv",row.names = 1)
SiteIndices<-merge(SiteIndices,SiteData,by="Site")
# this will discard any site indices that don't have site data available
# site data is UKBMS transects only, so non-transect site indices will be removed


### get orange tip (Anthocharis cardamines) data ###

# set up store for data
# this includes all siteyears for which we have any site indices (i.e. for any species) - these will be used to fill in zero counts
SiteYears<-SiteIndices[,-match(c("SPECIES.CODE","SPECIES","COMMON.NAME","SITE.INDEX"),colnames(SiteIndices))]
SiteYears<-unique(SiteYears)

OrangeTipSI<-SiteIndices[which(SiteIndices$COMMON.NAME=="Orange-tip"),] 

# check there's no NAs before we add zero counts
length(which(is.na(OrangeTipSI$SITE.INDEX)))

# merge orange tip site indices with table of all site years surveyed 
OTSiteIndices<-merge(SiteYears,OrangeTipSI,by=colnames(SiteYears), all.x=TRUE)

# assign zeros to those with no site index
OTSiteIndices$SITE.INDEX[which(is.na(OTSiteIndices$SITE.INDEX))]<-0

# remove those where the species was present but no site index could be generated (-2 values)
OTSiteIndices<-OTSiteIndices[-which(OTSiteIndices$SITE.INDEX==-2),] 

# add residual column
OTSiteIndices$resid<-1:nrow(OTSiteIndices) 

write.csv(OTSiteIndices,"UKBMS/Site indices/ukbmsSiteIndices2021_Anthocharis_cardamines.csv")



