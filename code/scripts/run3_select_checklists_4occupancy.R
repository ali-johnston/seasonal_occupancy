###########################################################################
###########################################################################
###### CODE TO RANDOMLY SELECT A SUBSET OF CHECKLISTS FOR OCCUPANCY MODELS
###### REDUCES CHECKLISTS TO MAX N CHECKLISTS PER SPATIAL-TEMPORAL UNIT
###########################################################################
###########################################################################

# created by: ali johnston


 # load in required packages
library(dplyr)


# set parameters for whole script
data.tag <- "BCR23_2016"


# File directories in which to search for models and plot results:
data.root <- "data/"
data.folder <- paste(data.root, data.tag, "_data/", sep="")

result.root <- "results/"
result.folder <- paste(result.root, data.tag, "/", sep="")
dir.create(result.folder, recursive=T)

# set working directory as overall folder
setwd("/Users/ali/Documents/REPOS/seasonal_occupancy")


# Source in relevant functions:
source("code/functions/functions_ebird_occupancy.R")
source("code/functions/functions_misc.R")



###########################################################################
###########################################################################
###### SET RUN PARAMETERS
###########################################################################
###########################################################################


#####################################################################
###  YEARS FOR ANALYSIS

min_yr <- 2016
max_yr <- 2016


#####################################################################
###  OCCUPANCY MODELLING PARAMS

# define minimum and maximum visits used per site. 
# if below minimum, site is not included
# if above maximum, a random selection of visits up to max is used
min_visits_site <- 2
max_visits_site <- 10



###########################################################################
###########################################################################
######   READ AND MANIPULATE EBIRD DATA
###########################################################################
###########################################################################


##################################################################### 
###  READ IN EXPERTISE DATA

# Read in observer expertise scores:
expertise.folder <- paste0(result.folder, "expertise/")

expertise_scores <- read.csv(paste0(expertise.folder, "observer_expertise.txt"), header=TRUE, sep=" ")



#####################################################################
###  READ IN BIRD DATA

# Collect data for these species
dat_loc <- paste0(data.folder, "BCR23_erd_gridcells.RData") # output from Common_dataset_processing.R
load(dat_loc)

covariates <- erd$X
locations <- erd$locs
birds <- erd$y

# combine together checklist info and bird data for selected species
all_bird_lists <- cbind(locations, covariates, birds)



#####################################################################
###  SELECT ONLY RELEVANT YEARS

bird_lists <- filter(all_bird_lists, YEAR<=max_yr, YEAR>=min_yr)

# tidy up workspace and remove large file 
rm(all_bird_lists)



#####################################################################
###  MERGE EXPERTISE SCORES WITH CHECKLIST INFO

# ensure both OBSERVER_IDs are in the same format
bird_lists$OBSERVER_ID <- as.character(bird_lists$OBSERVER_ID)
expertise_scores$OBSERVER_ID <- as.character(expertise_scores$OBSERVER_ID)
colnames(expertise_scores)[2] <- "expertise"

# merge 2 dataframes together by OBSERVER_ID
allspec_mod_data <- left_join(bird_lists, expertise_scores, by="OBSERVER_ID")



#####################################################################
###  ADD QUADRATIC FOR THE TIME OF DAY
allspec_mod_data$TIME2 <- allspec_mod_data$TIME^2


#####################################################################
###  tidy up workspace

rm(bird_lists)




###########################################################################
###########################################################################
######    SUBSET TO KEEP ONLY RANDOM MAX N CHECKLISTS PER SPATIAL UNIT
###########################################################################
###########################################################################


##################################################################### 
###  variables for analysis (this should be exhaustive - include all vars possible for modelling)

modis.numbers <- c(0:10, 12, 13, 16)
modis.umd.names <- paste0("UMD_FS_C", modis.numbers, "_1500_PLAND")
modis.umd.names <- update_modis_names(modis.umd.names)

site_covar_names <- c(modis.umd.names, "LATITUDE", "LONGITUDE", "ELEV")
obs_covar_names <- c("EFFORT_HRS", "TIME", "TIME2", "expertise", "EFFORT_DISTANCE_KM", 
   "NUMBER_OBSERVERS", "I.STATIONARY")


##################################################################### 
###  remove rows with incomplete information for any relevant covariates

cols_model <- allspec_mod_data[,c(site_covar_names, obs_covar_names)]

any_na <- apply(cols_model, 1, function(x){any(is.na(x))})

allspec_mod_data2 <- allspec_mod_data[!any_na,]


#####################################################################
###  tidy up workspace

rm(allspec_mod_data)




######################################################################################
### Filter Site-Observer Combinations with Too Few or Too Many Observation Records ###

allspec_mod_data2$LOC_ID <- allspec_mod_data2$grid.cell.number

#  RetainRangeCounts function randomly selects up to MaxNPerSiteObsrvr visits
# for each site. set.seed ensures the random selection is repeatable. 

# WARNING. 
# This step can be time-consuming for large datasets

set.seed(192)
mod_data <- RetainRangeCountsPerLocObsvrWeek(InputData = allspec_mod_data2, 
                                                      MinNPerSiteObsrvr = min_visits_site, 
                                                      MaxNPerSiteObsrvr = max_visits_site,
                                                      NWeek=1,
                                                      repeat_visit_vars=c("LOC_ID")) # "OBSERVER_ID", "YEAR", ...


##################################################################### 
###  randomly select 10% of the LOCATIONS for test data

all_loc <- table(allspec_mod_data2$LOC_ID)

set.seed(259)
test_loc <- names(all_loc)[sample(1:length(all_loc), round(length(all_loc)*0.1), replace=FALSE)]

mod_data$train_loc <- ifelse(mod_data$LOC_ID %in% test_loc, FALSE, TRUE)



#####################################################################
###  average landcover variables and elevation within each grid cell

landcover_vars <- c(modis.umd.names, "ELEV", "LATITUDE", "LONGITUDE")
landcover_only <- mod_data[,c(landcover_vars, "grid.cell.number")]

by_grid <- group_by(landcover_only, grid.cell.number)
avg_landcover <- as.data.frame(summarise_all(by_grid, mean))




#####################################################################
###  replace the landcover covariates for each checklist with the average landcover covariates per grid cell. 

mod_data_no_landcover <- select(mod_data, -one_of(landcover_vars))

# check that the grid cell columns are the same format
mod_data_no_landcover$grid.cell.number <- as.character(mod_data_no_landcover$grid.cell.number)
avg_landcover$grid.cell.number <- as.character(avg_landcover$grid.cell.number)

# merge together
mod_data2 <- left_join(mod_data_no_landcover, avg_landcover, by="grid.cell.number")

write_loc <- paste0(data.folder, "Checklist_by_grid.cell_week_max10vis.txt")
write.table(mod_data2, write_loc, row.names=FALSE)




