##################################################################### 
###  CODE FROM WES, STRIPPED DOWN TO KEEP ONLY THE FUNCTIONS TO SOURCE
###  IN FROM OTHER SCRIPTS


#The R script contained in this file provides an example of using eBird data to fit
# occupancy models.  The example makes use of the same data from Wood Thrushes that
# Daniel Fink has created and used in his examples of using eBird data in species
# distribution models.  Specifically, the data come from the following two data.frames
# "testing.data" and "training.data" that are being distributed in the two files
# "Wood_Thrush_testing.data.RData" and "Wood_Thrush_training.data.RData"; the data 
# within these two data.frames are combined to form the basis for the occupancy
# modelling example.
#The work done within this script is divided into four parts:
#  (1) Filtering the data to remove rows that do not meet some basic criteria
#       to be included (count protocol, count duration, distance travelled, and
#       number of observers).  This filtering is hard-wired to retain only data
#       from stationary counts (i.e. point counts) and travelling counts (i.e.
#       transect counts), but the distance, duration, and number of observers
#       criteria can be set by the analyst.
#  (2) Filtering the data based on the number of repeated visits by the same 
#       observer to the same location.  The analyst can specify both the minimum
#       and maximum allowable number of visits to be used for data within the
#       final occupancy model analyses.
#  (3) Reformatting of the data so that they conform to the rather persnickety
#       format that is required by the "unmarked" package in R.  In general terms
#       this reformatting involves taking data that are organized so that each
#       visit to a site by an observer is described in a single row of data, so 
#       that the resultant table of data has all of the information from the
#       repeated visits (of a single observer at a single site) arranged within
#       a single row.  There are a series of essential details to this process
#       including naming conventions for the variables in the table that is
#       produced.
#  (4) Fitting occupancy models and examining the results.
# Model filtering and reformatting are done by three R functions that are 
# created within this script.
#NOTE: there are a small number of changes that each user will need to make
# in order to run the example on their own computer.  You can find these locations
# by searching for the text phrase "CHANGE REQUIRED".  At this time (6 August 2016), 
# there is only a single location at which changes will need to be made to two lines
# of the script.
#NOTE: variable names are those used in the 2014 version of the eBird ERD for all
# of the variables in the input data. The data filtering and formating sections of 
# this script depends on the existence 6 variables with specific names, as follows:
#  - COUNT_TYPE
#  - EFFORT_HRS
#  - EFFORT_DISTANCE_KM
#  - NUMBER_OBSERVERS
#  - LOC_ID
#  - OBSERVER_ID
#  - DAY
#  - TIME
#NOTE: The example data used with this script has already been partially filtered
# and does not represent the full range of variation in the eBird Reference Dataset
# (ERD).  Specifically, prior to the start of this script the following has already
# been done to create a subset of the entire ERD:
#  - only data from a limited range of dates within a single year have been retained
#  - only a single checklist from a group of shared checklists has been retained
#  - the response variable has been converted from the counts of bird of the focal
#     species (or "X" values that indicate that at least one bird of the species 
#     was present) into a binary, presence/absence variable: zero counts are kept 
#     as zeros (i.e. absences), but all count greater than one and all "X" values 
#     have been converted to "1" (i.e. presence).  This script is perfectly content
#     with having multiple species' presence/absence response data included in the
#     input data, because the species whose data are to be analysed needs to be
#     specified in order for the data-formatted functions to run (other species
#     are simply ignored)
#This script was created by Wesley Hochachka in early August 2016, based on scripts
# previously created by Wesley Hochachka with some bug fixes provided by Richard
# Schuster.



######################################################################
### Filter Out Unwanted Data by Column Values & Low Row Prevalence ###
######################################################################
#Define a function that filters the input data, retaining only those rows (checklists)
# and columns (species) that the analyst deems suitable for further use.  This
# filtering is done based on the following criteria:
#  (1) Count protocol: hardwired removal of rows not from stationary or travelling 
#       counts.  Note that this is not relevant for the specific example.  However,
#       more generally for use of this function with a less filtered version of the
#       eBird ERD data, I have included this filtering.  This function tests for
#       the existence of the variable that identifies the count protocol, the 
#       ERD variable "COUNT_TYPE", and if this variable exists within the data.frame
#       then removal of all but stationary or travelling count data is done.
#  (2) Count duration: removal of counts longer in time than user-specified maximum
#  (3) Count distance: removal of counts longer than user-specified distance
#  (4) Number observers: remove checklists with more than some maximum number of 
#       observers listed.  These numbers can be in excess of 50 observers listed(!),
#       but typically the numbers range from between 1 and 5 observers, at least in
#       the data (California & New York, June & December 2013) examined in prototyping
#       this function.  There are some records for which the number of observers is
#       not provided, and to deal with these rare cases, we also have to test whether
#       there is a number provided.
# As inputs, this function requires the name of the input data.frame, the
# threshold maximum count duration in hours, the maximum threshold count
# distance travelled in kilometers, and the minimum proportion (i.e. zero to one
# value) of checklists on which a species is reported.
#Additionally, make a slight change to content of the field describing count
# protocols, by turning the code values into human-intelligible labels.
FilterChecklists <- function(InputData, MaxHrs, MaxDistKm, MaxNObservers){
   #The first task for this function is to remove data from protocols other than
   # stationary ("P21") and travelling ("P22"), and counts from longer times or
   # distances than specified in the input variables.  This filtering is done by
   # appending a variable to be treated as boolean to the input data.frame, and
   # then running each row through a set of 3 tests that will flip the boolean 
   # from "Keep" to "remove" states (i.e. from the default value of 1 to zero).
   Keep <- rep(1, nrow(InputData))
   InputData <- cbind(InputData, Keep)
   rm(Keep)
   #
   if ("COUNT_TYPE" %in% names(InputData))
      InputData[InputData$COUNT_TYPE != "P21" & InputData$COUNT_TYPE != "P22",]$Keep <- 0
   #
   InputData[InputData$EFFORT_HRS > MaxHrs,]$Keep <- 0
   #
   InputData[InputData$EFFORT_DISTANCE_KM > MaxDistKm,]$Keep <- 0
   #
   InputData[InputData$NUMBER_OBSERVERS > MaxNObservers & !is.na(InputData$NUMBER_OBSERVERS),]$Keep <- 0
   #
   Filtered <- InputData[InputData$Keep == 1,]
   #
   #We do not need the variable "Keep" anymore, so let's tidy up by removing it
   # with the following command.  This is a bit convoluted looking, and has to
   # be understood by reading the command from the outside in.  Start with the
   # inner bit saying 'which(names(Filtered) %in% "Keep")'.  This part of the
   # command gives us the column number for the variable with the name "Keep"
   # (by returning the place number from the entire list of names of columns
   # that has the name "Keep").  All of this is embedded within a pair of 
   # square brackets, which are used to specify a subset of rows and columns.
   # rows are specified by information that comes before the comma within the 
   # brackets, and columns after the comma; so with nothing before the comma
   # we are telling R that we want to select all of the available rows of data,
   # and so it is just the columns for which we want a subset.  The "which..."
   # command selects a specific column number or column numbers, and by adding
   # a negative sign before this we are telling R to remove just the specified
   # column numbers.  So when all of this is combined, we are telling R to
   # select all of the rows and all of the columns of data except for the 
   # column named "Keep".  We are overwriting the original version of the
   # data.frame named "Filtered" with the version from which this single
   # column has been removed.
   Filtered <- Filtered[, -(which(names(Filtered) %in% "Keep"))]
   #Now change the labels for the count protocols to make sense to mere mortals.
   # Oh, and test to see if this is necessary, before doing the renaming.
   if ("COUNT_TYPE" %in% names(InputData)){
      Filtered[Filtered$COUNT_TYPE == "P21",]$COUNT_TYPE <- "stationary"
      Filtered[Filtered$COUNT_TYPE == "P22",]$COUNT_TYPE <- "travelling"   
   }
   #Send completed data.frame back to main program.
   return(Filtered)
}



######################################################################################
### Filter Site-Observer Combinations with Too Few or Too Many Observation Records ###
######################################################################################
#Now we have to deal with the instances in which an observer has submitted a larger 
# than MaxNPerSiteObsrvr number of checklists.  The approach that I am taking is to 
# select a random subset of these lists to use, discarding the others.  We start by
# splitting the data into two groups: those person-site combinations with up to the 
# the maximum number of checklists, and those with more than the maximum.  Then, for 
# this latter group of person-site combinations we assign a random number to 
# each of the records, and finally delete all ranks above the maximum number that 
# we want to keep.  Unfortunately, R's "rank" function is global in scope, not 
# allowing us to easily assign ranks to each separate observer-location combination.  
# As best as I can tell we can't easily take advantage of an SQL-based solution 
# because RSQLite also does not have a "rank" function, so we need to use 
# "aggregate" to do the job for us.  Using "aggregate" is somewhat painful, 
# because "aggregate" returns a separate list of ranks for each person-location-year, 
# and I need to convert these to a vector using the "unlist" function before 
# adding the new column and removing the rows with rank values over the maximum 
# allowed.  Additional funkiness is that aggregate seems to be ordering its output 
# in alphabetical order of the factor numbers of a factor version of the grouping 
# variable, PersonLocYr, and not of the character-value alphabetical order.  So, 
# I need to insure that the data in DataSubset_ExcessLists are ordered by this
# same way.  Only then can I concatinate the two subsets of rows back together 
# into a single data.frame, and then recalculate the NTimesPerLoc values for 
# this subset of the data.  From this we can finally build the data objects 
# for occupancy modelling.  
#While we are at it, we will also remove the data from location-observer 
# combination with fewer than the minimum number of lists in the set of data.
#All of this is wrapped into function so that it can be applied to multiple, 
# independent subsets of data.  We set default values of the minimum and maximum
# numbers of checklists for each site-observer combination to 1 and 10, respectively
# just because that's what I'm currently doing with my actual work while I am 
# developing this script.  Without defaults in place, if someone forgets to
# specify these values the script doesn't generate an error, but the resultant 
# output file is missing the NLists column (for some reason I have not bothered to
# try to understand).
RetainRangeCountsPerLocObsvr <- function(InputData, MinNPerSiteObsrvr=1, MaxNPerSiteObsrvr=10,
                                          repeat_visits=c("LOC_ID", "OBSERVER_ID")){
   #
   #Create a new variable, "ID", that is a single-variable identifier of unique
   # combinations of location and observer, needed because we will be grouping
   # individual observations into records groups of repeated records based on
   # this combination of location and observer.  We concatinate the identifiers
   # of location (LOC_ID) and observer (OBSERVER_ID) with an underscore separating
   # them, and name this new column of information "ID".  Then we bind this new
   # column at the front of the data.frame.
   ID <- InputData[,which(colnames(InputData)==repeat_visits[1])]
   if(length(repeat_visits)>1){
      for(i in 2:length(repeat_visits)){
         ID <- paste(ID, InputData[,which(colnames(InputData)==repeat_visits[i])], sep="_")
      }
   }
   ID <- as.data.frame(ID, stringsAsFactors=FALSE)
   names(ID) <- "ID"
   InputData <- cbind(ID, InputData)


   time_stamp1 <- Sys.time()
      #First, we need to add a count of the number of checklists in the data set for each
   # combination of location (LOC_ID) and observer (OBSERVER_ID).  The four steps are:
   # (1) use the "aggregate" function to count the times that each LOC_ID - OBSERVER_ID
   # combination was seen, (2) name the columns in the output for merging, (3) merge the
   # new information to append the counts to the ends of the original data lines, and
   # (4) tidy up.
   NPerLoc <- aggregate(InputData$LOC_ID, by=list(InputData$ID), FUN="length")
   names(NPerLoc) <- c("ID", "NLists")
   InputData <- merge(x=InputData, y=NPerLoc, by=c("ID"))
   rm(NPerLoc)
   #
   #Now, divide the data into two subsets based on whether we need to remove excess lists.
   DataSubset_EnoughLists <- InputData[InputData$NLists <= MaxNPerSiteObsrvr,]
   DataSubset_ExcessLists <- InputData[InputData$NLists > MaxNPerSiteObsrvr,]
   #
   #For the subset that needs trimming, add a column of random numbers to be converted 
   # into ranks, and then sort the original data by a factor-coded version of the location
   # ID and observer ID to insure that we can concatinate the column of rank information back 
   # onto the original data.  The data should already be sorted appropriately, so this
   # sorting is a bit of paranoid self defense.
   DataSubset_ExcessLists <- cbind(DataSubset_ExcessLists, runif(n=nrow(DataSubset_ExcessLists)))
   names(DataSubset_ExcessLists) <- c(names(InputData), "RandNum")
   DataSubset_ExcessLists <- DataSubset_ExcessLists[order(as.factor(DataSubset_ExcessLists$LOC_ID), as.factor(DataSubset_ExcessLists$OBSERVER_ID)),]
   #
   #Create the ranks based on location and observer groups, and then extract this information from
   # the list object that is returned by "aggregate" in the form of a vector using "unlist".  Note 
   # that even though the data have been sorted first by LOC_ID and then by OBSERVER_ID, these two
   # sorting variables are given in reverse order in the "by" option of aggregate.  We need to do
   # this, because the list produced as output to "aggregate" is sorted in reverse order to the
   # order in which the variables are stated in the "by" list.
   RanksList <- aggregate(x=DataSubset_ExcessLists$RandNum, 
                       by=list(DataSubset_ExcessLists$ID),
                       FUN="rank", simplify=T)
   Ranks <- unlist(RanksList$x)
   #
   #Join the column of ranks back onto the original data, and remove those records with
   # ranks greater than the maximum desired number of records per location-observer.
   # Remove the two columns (random number and rank) that are no longer needed, and
   # add back the records onto the subset of the data that was not in need of removing
   # excess records.
   DataSubset_ExcessLists <- cbind(DataSubset_ExcessLists, Ranks)
   DataSubset_ExcessLists <- DataSubset_ExcessLists[DataSubset_ExcessLists$Ranks <= MaxNPerSiteObsrvr,]
   DataSubset_ExcessLists <- DataSubset_ExcessLists[, 1:(ncol(DataSubset_ExcessLists) - 2)]
   DataSubset_new <- rbind(DataSubset_EnoughLists, DataSubset_ExcessLists)
   #
   #Tidy up.
   rm(RanksList, Ranks, DataSubset_EnoughLists, DataSubset_ExcessLists)
   #
   #Trim off the data from location-observer combinations with fewer than specified minimum
   # number of lists.
   DataSubset_new <- DataSubset_new[DataSubset_new$NLists >= MinNPerSiteObsrvr,]
   DataSubset_new <- DataSubset_new[DataSubset_new$NLists <= MaxNPerSiteObsrvr,]


   #The counts NLists are now wrong for any of the site-observer combinations
   # for which only a random subset of size MaxNPerSiteObsrvr was retained.
   # So, we need to recount the numbers of checklists per site-observer
   # combination, before returning the updated data set for further use.
   NPerLoc <- aggregate(DataSubset_new$ID, 
                     by=list(DataSubset_new$ID), 
                     FUN="length")
   names(NPerLoc) <- c("ID", "NLists_new")
   DataSubset_new <- merge(x=DataSubset_new, y=NPerLoc, by=c("ID"))
   DataSubset_new <- DataSubset_new[,-(which(names(DataSubset_new) == "NLists"))]
   names(DataSubset_new)[which(names(DataSubset_new) == "NLists_new")] <- "NLists"
   rm(NPerLoc)  
   #Ship out the final data.frame for use outside of the function.


   time_stamp2 <- Sys.time()
   print(time_stamp2 - time_stamp1)
   return(DataSubset_new)
}


###########################################################################
###### ADAPTATION BY TO ALLOW TEMPORAL SPECIFICATION OF 'SITE'

RetainRangeCountsPerLocObsvrWeek <- function(InputData, MinNPerSiteObsrvr=1, MaxNPerSiteObsrvr=10, NWeek=2,
   repeat_visit_vars=c("LOC_ID", "OBSERVER_ID")){
   # specify the week of the observation
   week <- ceiling(InputData$DAY/7)
   week[week>52] <- 52 # cap at week 52

   # group together in number of weeks
   InputData$week_cat <- ceiling(week/NWeek)

   #Create a new variable, "ID", that is a single-variable identifier of unique
   # combinations of location and observer, needed because we will be grouping
   # individual observations into records groups of repeated records based on
   # this combination of location and observer.  We concatinate the identifiers
   # of location (LOC_ID) and observer (OBSERVER_ID) with an underscore separating
   # them, and name this new column of information "ID".  Then we bind this new
   # column at the front of the data.frame.
   ID <- InputData[,which(colnames(InputData)==repeat_visit_vars[1])]
   if(length(repeat_visit_vars)>1){
      for(i in 2:length(repeat_visit_vars)){
         ID <- paste(ID, InputData[,which(colnames(InputData)==repeat_visit_vars[i])], sep="_")
      }
   }

   # add in week category
   ID <- paste(ID, InputData$week_cat, sep="_")

   # convert to a data frame
   ID <- as.data.frame(ID, stringsAsFactors=F)
   names(ID) <- "ID"
   InputData <- cbind(ID, InputData)
#   rm(ID)


   time_stamp1 <- Sys.time()
      #First, we need to add a count of the number of checklists in the data set for each
   # combination of location (LOC_ID) and observer (OBSERVER_ID).  The four steps are:
   # (1) use the "aggregate" function to count the times that each LOC_ID - OBSERVER_ID - week_cat
   # combination was seen, (2) name the columns in the output for merging, (3) merge the
   # new information to append the counts to the ends of the original data lines, and
   # (4) tidy up.
   NPerLoc <- aggregate(InputData$LOC_ID, by=list(InputData$ID), FUN="length")
   names(NPerLoc) <- c("ID", "NLists")
   InputData <- merge(x=InputData, y=NPerLoc, by=c("ID"))
   rm(NPerLoc)
   #
   #Now, divide the data into two subsets based on whether we need to remove excess lists.
   DataSubset_EnoughLists <- InputData[InputData$NLists <= MaxNPerSiteObsrvr,]
   DataSubset_ExcessLists <- InputData[InputData$NLists > MaxNPerSiteObsrvr,]

   #
   #For the subset that needs trimming, add a column of random numbers to be converted 
   # into ranks, and then sort the original data by a factor-coded version of the location
   # ID and observer ID to insure that we can concatinate the column of rank information back 
   # onto the original data.  The data should already be sorted appropriately, so this
   # sorting is a bit of paranoid self defense.
   DataSubset_ExcessLists <- cbind(DataSubset_ExcessLists, runif(n=nrow(DataSubset_ExcessLists)))
   names(DataSubset_ExcessLists) <- c(names(InputData), "RandNum")
   DataSubset_ExcessLists <- DataSubset_ExcessLists[order(as.factor(DataSubset_ExcessLists$LOC_ID), as.factor(DataSubset_ExcessLists$OBSERVER_ID)),]
   #
   #Create the ranks based on location and observer groups, and then extract this information from
   # the list object that is returned by "aggregate" in the form of a vector using "unlist".  Note 
   # that even though the data have been sorted first by LOC_ID and then by OBSERVER_ID, these two
   # sorting variables are given in reverse order in the "by" option of aggregate.  We need to do
   # this, because the list produced as output to "aggregate" is sorted in reverse order to the
   # order in which the variables are stated in the "by" list.
   RanksList <- aggregate(x=DataSubset_ExcessLists$RandNum, 
                       by=list(DataSubset_ExcessLists$ID),
                       FUN="rank", simplify=T)
   Ranks <- unlist(RanksList$x)
   #
   #Join the column of ranks back onto the original data, and remove those records with
   # ranks greater than the maximum desired number of records per location-observer.
   # Remove the two columns (random number and rank) that are no longer needed, and
   # add back the records onto the subset of the data that was not in need of removing
   # excess records.
   DataSubset_ExcessLists <- cbind(DataSubset_ExcessLists, Ranks)
   DataSubset_ExcessLists <- DataSubset_ExcessLists[DataSubset_ExcessLists$Ranks <= MaxNPerSiteObsrvr,]
   DataSubset_ExcessLists <- DataSubset_ExcessLists[, 1:(ncol(DataSubset_ExcessLists) - 2)]
   DataSubset_new <- rbind(DataSubset_EnoughLists, DataSubset_ExcessLists)
   #
   #Tidy up.
   rm(RanksList, Ranks, DataSubset_EnoughLists, DataSubset_ExcessLists)
   #
   #Trim off the data from location-observer combinations with fewer than specified minimum
   # number of lists.
   DataSubset_new <- DataSubset_new[DataSubset_new$NLists >= MinNPerSiteObsrvr,]
   DataSubset_new <- DataSubset_new[DataSubset_new$NLists <= MaxNPerSiteObsrvr,]


   #The counts NLists are now wrong for any of the site-observer combinations
   # for which only a random subset of size MaxNPerSiteObsrvr was retained.
   # So, we need to recount the numbers of checklists per site-observer
   # combination, before returning the updated data set for further use.
   NPerLoc <- aggregate(DataSubset_new$ID, 
                     by=list(DataSubset_new$ID), 
                     FUN="length")
   names(NPerLoc) <- c("ID", "NLists_new")
   DataSubset_new <- merge(x=DataSubset_new, y=NPerLoc, by=c("ID"))
   DataSubset_new <- DataSubset_new[,-(which(names(DataSubset_new) == "NLists"))]
   names(DataSubset_new)[which(names(DataSubset_new) == "NLists_new")] <- "NLists"
   rm(NPerLoc)  
   #Ship out the final data.frame for use outside of the function.


   time_stamp2 <- Sys.time()
   print(time_stamp2 - time_stamp1)
   return(DataSubset_new)
}





#######################################################################################
################# Creating Data.Frame Suitable for Use With unmarked ##################
#######################################################################################
#
#Once the primordial data.frame has been created, we need to take the one-observation-
# per-row data, and turn them into a form in which the repeated observations from a 
# single location form a series on a single row.  That is the work for this next section 
# of the script. Our presumed target will be a data.frame object that can be fed into 
# the unmarked package.  While it is possible for unmarked to use one-observation-per-row 
# data, using unmarked's "formatLong" function, this is only possible of there are
# no site-level covariates (i.e. predictors specific to a site, e.g., elevation, that
# do not change from observation to observation within that site).  Because this
# is not always the case, the script below will create a data.frame that is more
# universally useable, to be processed by the "formatWide" function of unmarked
# for subsequent use by functions in the "unmarked" package.
#
#In our example, we want to fit a model in which either detection or
# occupancy rates are dependent on the following predictors (the various 
# forest types' percentage cover can probably be summed prior to analyses): 
#
#  DAY - continuous variable, Julian date
#  TIME - continuous variable, decimal hours since midnight
#  COUNT_TYPE - factor, whether stationary or travelling count
#  EFFORT_HRS - continuous variable
#  EFFORT_DISTANCE_KM - continuous variable, value of zero for stationary counts
#  NUMBER_OBSERVERS - continuous variable
#  NLists - the number of repeated counts during the period of interest
#  Deciduous_broad_PLAND - an estimate of the propostion of land in deciduous broadleaf forest within 1500m of the site
#  Urban_Built_PLAND - an estimate of the proportion of urbanized landscape within 1500m of the site
#
# Additionally, we will have the identifier ID, created as a combination of
# location and observer ID numbers, as an identifier label at the start of each
# line.
#Aside from the ID variable and the habitat information, each of the other variables 
# (the response and the 6 predictors) will need to appear N times in each row of the 
# data file that is being created, where N is the maximum number of observation periods 
# specified by MaxNPerSiteObsrvr, that was set immediately above.  We will create this
# data file by first creating an empty template record, and after that we will
# go through our one-observation-per-row data line by line and population this
# template, starting a new row for each separate ID number.
#
#In order to create both create a template to be filled with data from each
# site, and to populate the template, we need to know the names of the
# variables in the data to be reformatted.  We will be providing this information
# to R in the form of 4 vectors of character data:
#  - The name of the column containing the site-identifier information.  This
#     is a vector with only a single item.  If there is no site-identifier
#     variable (this is optional for data going into unmarked), simply do not
#     provide any information to the function about SiteIDName.
#  - The name of the column containing the response (detection/non-detection
#     data).  This is a vector with only a single variable.  The name of this
#     variable is required.
#  - The names of the site-level covariates.  One or more items can be in
#     this vector.  This is optional, and if no site-level covariates are
#     needed, simply do not mention SiteCovarSet when calling the function.
#  - The names of the observation-level covariates.  Again, one or more 
#     names can be placed in this vector.  This is also optional, and when
#     no observation-level covariates are needed, omit mention of ObsCovarSet
#     when calling the function.
#NOTE: the names in these lists must be written exactly as the variables
# are named in the data.frame that is to be reformatted.  Failure to 
# have these names exactly correct will result in bad things happening...or 
# at least good things not happening.
#NOTE: This function is *almost* complete agnostic to the data fed into it.  I.e.,
# with just some slight modification, you could use this function to format data
# from any source so that the data of one-observation-per-line data into the
# all-observation-from-a-site format needed by the "unmarked" package.  There are
# only two modifications needed (either to another source of data, or to this
# function) in order to have it work with sources of data other than the eBird
# Reference Dataset (ERD):
#  (1) The data need to have one or more columns on which the rows can be sorted
#       so that all of the data that to be combined into a single record of
#       repeated observations can be sorted so that they form a contiguous block
#       or rows.  This function sorts data into these blocks based on the 
#       variables named: ID, DAY, and TIME.  Search the function for the text
#       "InputData[order(Inputdata" in order to find the location at which
#       this sorting happens.  You really only need to sort on the single grouping
#       identifier "ID"; sorting within each ID group by date and time-within-day
#       is just a nicety that results in the output data.frame having the 
#       individual observations being presented in chronological order.  For
#       non-eBird ERD data, either rename the variables for group, date, and time
#       within this function to match those in your data, or alter the names in 
#       your data to match those used by the function.  Minimally, you will need
#       a variable in your data that corresponds to "ID" in order to tell the
#       function which rows of input data to group together into a single output
#       record.  If you do not have variables that correspond to date and time-
#       within-day, then you can safely delete mention of these two sorting 
#       variables from the script.
#  (2) The function depends on a variable named "NLists" in order to know when it is
#       finished with assembling all of the data within a group into a single
#       row of output.  "NLists" is a variable that tells the script how many
#       observation sessions were conducted at each of the locations by an
#       observer.  This variable is needed in order to allow the function to 
#       work properly when varying numbers of visits were made to each of the
#       sites: if fewer than the maximum number of surveys were conducted at
#       a specific site by an observer, the output record is padded with 
#       missing values.  The occupancy model fitting using "unmarked" has no
#       problem dealing with varying numbers of repeat visits per site.  In 
#       order to use this function with some other source of input data, you 
#       would either need to create the column of information "NLists" and add
#       it to your data, rename the variable within the script to match the
#       name of the variable within your own data (i.e. search and replace
#       all instances of "NLists" within this function), or reword this script
#       in order to create another mechanism for identifying when the end of 
#       a group of input records has been reached and it's time to print out
#       the assembled data for this one group and start assembling a record
#       for the next group of observations.
#
#Here we tell R that we are creating a function to produce an output file in the 
# format of an unmarked "wide" input data.frame.  Note that all non-essential 
# information is by default given the value of NULL (i.e. the information 
# is assumed not to exist unless explicitly mentioned when calling the 
# function).  I.e., you do not need to have a variable that carried into
# the output that identifies your site, nor do you need to have any 
# site-level or observation-level predictor variables.
MakeUnmarkedWideInput <- function(InputData, 
                                  SiteIDName = NULL, 
                                  ResponseName, 
                                  SiteCovarSet = NULL, 
                                  ObsCovarSet = NULL, 
                                  MaxNPerSiteObsrvr){
   #
   #Prepare the input data for formatting for unmarked.  First sort the data
   # in the order required for reformatting.
   InputData <- InputData[order(InputData$ID, 
                                InputData$DAY, 
                                InputData$TIME),]
   #
   #Assemble a blank template row containing appropriately labelled columns.
   # Start by creating a single-item data.frame from which to build units 
   # of repeated counts and their covariates.
   BuildingBlock <- data.frame(x = NA)
   #
   #From BulidingBlock create a template for the site identifier, and another for
   # the vectors of response variable and observation-specific covariations that
   # have the length of the number of maximum number of possible observation 
   # periods.
   SingletonTemplate <- BuildingBlock
   #
   GroupTemplate <- BuildingBlock
   for (i in 2:MaxNPerSiteObsrvr){   
      GroupTemplate <- cbind(GroupTemplate, BuildingBlock)
   }
   rm(i, BuildingBlock)
   #
   #Copy these templates to form the base on which to build a complete
   # record of the response plus observation-level covariates, filled 
   # with NA values.  For use with unmarked, the format for data should
   # be the following:
   #  (i)   The site-identifier variable (optional) named "site"
   #  (ii)  The sets of response-variable values for all observations at
   #         that site.
   #  (iii) Columns for site-level covariates (i.e. covariates whose values 
   #         do not change among observation periods)
   #  (iv)  Sets of columns for the observation-level covariates (i.e.
   #         covariates that can change in value among the observation 
   #         periods within each site).  For each covariate there needs to
   #         be a set equal in length to the maximum number of observation
   #         periods.  The names of covariates within each set have the same
   #         prefix, a dot, and then na integer number from 1 to the maximum
   #         number (e.g., Time.1, Time.2, Time.3 ... Time.n).
   # Our first step in building this structure will be to assemble together
   # parts (i) and (ii).  Because the site label/column is optional for data
   # intended for unmarked, we need to test if there is an identifier variable.
   # We do this by testing if the R object "SiteIDName" exists.
   if (!is.null(SiteIDName)){
      BlankRecord <- cbind(SingletonTemplate, GroupTemplate)
   } else if (is.null(SiteIDName))
      BlankRecord <- GroupTemplate  
   #
   #Assign names of individual items in the GroupTemplate, initially for the
   # response variable, following the naming convention that unmarked
   # will need: variable name, followed by a dot, followed by the integer
   # number of that observation in the sequence.  As before, we test to 
   # make sure that we actually have a variable to be used as a site
   # identifier, and the counter used in the for-loop for creating the
   # column name labels deals with the possible absence of the site 
   # identifier by using the "1+j" construct.
   if (!is.null(SiteIDName))
      names(BlankRecord)[1] <- "site"
   ColNameBase <- "y"
   if (!is.null(SiteIDName)){     
      j <- 1
   } else j <- 0
   for (i in (1):(MaxNPerSiteObsrvr)){
      names(BlankRecord)[i+j] <- paste(ColNameBase, (i), sep=".")
   }
   rm(i, j, ColNameBase)
   #
   #Now we add component (iii), the site-level covariates and name them. 
   # The list of site-level covariates needs to have been specified already 
   # in a vector of names SiteCovarSet.  If there are no site-level 
   # covariates specified, this step is skipped, so the option exists to
   # not have any site-level covariates (this functionality was suggested by
   # Richard Schuster).  Note that we are cutting a corner a wee bit, but 
   # recycling and renaming the single variable in SingletonTemplate for 
   # each of the site-level covariates.
   if (!is.null(SiteCovarSet)){
      for (i in 1:length(SiteCovarSet)){
         ColName <- SiteCovarSet[i]
         names(SingletonTemplate) <- ColName
         BlankRecord <- cbind(BlankRecord, SingletonTemplate)
      }
      rm(i, ColName)     
   }
   #
   #Now we add the parts of component (iv) to the structure. We add each
   # observation-level covariate in turn as its complete set and name the
   # observation-specific covariate variable.  The list of names of
   # observation-level covariates needs to have already been specified in
   # the vector ObsCovarSet.  It is possible to not use any site-level
   # covariates by simply not specifying any covariate names (modification
   # suggested by Richard Schuster). Note that we are cutting corners a bit
   # by recycling and renaming the components of GroupTemplate for each of
   # the observation-level covariates in turn.
   if (!is.null(ObsCovarSet)){
      for (i in 1:length(ObsCovarSet)){
         ColNameBase <- ObsCovarSet[i]
         for (j in 1:MaxNPerSiteObsrvr){
            names(GroupTemplate)[j] <- paste(ColNameBase, j, sep=".")
         }
         BlankRecord <- cbind(BlankRecord, GroupTemplate)  
      }
      rm(i, j, ColNameBase, SingletonTemplate, GroupTemplate) 
   }
   #
   #
   #Now we can assemble our data.frame containing one row of data for the
   # complete set of observations from each of the site-observer combinations.  
   # We do this by using an outer "while" loop that moves through the input data 
   # for each site-observer set of lines in turn.  For each new site-observer 
   # combination the first data to population this new record are the site ID 
   # variable and then information from the site-level covariates.  We then  
   # step through the records of this site with a middle "for" loop that walks 
   # through the lines that contain data from each observation periods for that 
   # specific location-observer combination, making  use of NLists (the variable 
   # that tells the number of observation periods) in the middle loop's counter 
   # of the records for each site-observation combination.  Within each of the 
   # input lines we need a final, inner-most "for" loop that steps through each 
   # observation period's response (detected/not detected) information and the 
   # observation-level covariates in turn within each of the lines specified 
   # in the middle loop.
   #Note that this process only works when the input data are sorted by the 
   # ID column (the concatenated information from LOC_ID and OBSERVER_ID).  
   # While this should already be the case, out of an overabundance of paranoia 
   # we will be sorting the data again, and in doing this sorting not only by 
   # the ID variable but also by date and time of day within ID (mostly just 
   # because it feels more orderly to have the final record presenting 
   # observations in chronological order).
   #
   #Check if there is a column containing "site" identifier information.  We do 
   # this by first assuming that no such information exists, but then if there
   # is a value read into the function for SiteIDName, then we change update
   # this assumed "no" to a "yes".  The variable IsSiteName is numeric, because 
   # we need to use this information in order to correctly place the information
   # into columns in the table, depending on whether there is a first column of
   # site ID data.
   IsSiteName <- 0
   if (!is.null(SiteIDName))
      IsSiteName <- 1
   #
   #We use "i" as the counter in the outer loop, working down through the first 
   # lines of each group of input records, starting with the very first line of the 
   # input file. We only stop when the next possible line number is beyond the 
   # end of the file.  This starting and stopping is set and controlled by the 
   # next two lines.
   i <- 1
   while (i <= nrow(InputData)){        #outer loop start: stepping through unique site/observer combinations
      NewWideRow <- BlankRecord          #create blank record to be populated with data
      if (!is.null(SiteIDName)){         #if there is site ID information, add this to record being built
         NewWideRow[1] <- InputData[i, which(names(InputData) == SiteIDName)]
      }
      if (!is.null(SiteCovarSet)){
         for (s in 1:length(SiteCovarSet)){            #population record with site-level covariate information
            NewWideRow[IsSiteName+MaxNPerSiteObsrvr+s] <- InputData[i, which(names(InputData) == SiteCovarSet[s])]
         }
         rm(s)
      }
      for (j in (1):(InputData[i,]$NLists)){   #middle loop start: each observation in site-observer combination
         NewWideRow[j+IsSiteName] <- InputData[i+j-1, which(names(InputData) == ResponseName)]   #response variable data
         if (!is.null(ObsCovarSet)){
            for (k in 1:length(ObsCovarSet)){       #inner loop: each obs-level covariate within site-observer
               NewWideRow[j+IsSiteName+length(SiteCovarSet)+(MaxNPerSiteObsrvr*k)] <- 
                  InputData[i+j-1, which(names(InputData) == ObsCovarSet[k])]
            }
            rm(k)
         }
      }
      #Create new table if this is the first output record; the "ifelse" clause added by Richard Schuster in
      # order to fix an error that would cause loss some rows of data.  I had previously just written
      # "InputData[i,]$NLists)" as the condition after "==".  Anyway, here's what the "if" statement does.
      # We need to create a new table only for the first fully assembled row that is intended for this table,
      # and that is where the "i == 1" condition comes it: this test is true only for the first record because 
      # the counter variable "i", which is a counter of the row number for the starting row of each site-observer
      # combination, can only have a value of one for the first record.  However, we also need to wait until this
      # first record is completely filled in before we save it to the file that we are building, and this test of
      # completion is the responsibility of the second part of the logical test, the "(j) ==" clause. "j" is 
      # the counter of the number of observation periods whose data have been added to the row that is being built,
      # used in the "for" loop that is run in the part of the script immediately above.  Knowing the number of 
      # observation periods whose data are to be added to this row, we compare this value "j" against the number of
      # observation periods whose data should have been incorporated into the row as a double-checking whether this 
      # row of data is complete.  If the number of observation periods of available data was less than
      # or equal to the maximum number of observation periods allowed for each site-observer, then we just want to
      # compare the value of "j" to the count of the number of lists, the "NLists" variable.  However, we also have
      # the special case in which there are more observations than allowable time periods, and in that case we need
      # to compare "j" to the maximum allowable number of observation periods.  So, "j" needs to be compared with
      # one of two things, depending on how many observation periods' data were available, and the selection of
      # which number to compare to "j" is the job of the "ifelse" (note the lack of a space!) function.
      if (i == 1 & (j) == ifelse(InputData[i,]$NLists > MaxNPerSiteObsrvr, MaxNPerSiteObsrvr, InputData[i,]$NLists)){
         NewWideTable <- NewWideRow
         NewWideRow <- BlankRecord    
      #Append latest record to output table, if it exists; the "ifelse" clause added by Richard Schuster to fix 
      # an error in my origional version that did not correctly account for site-observer combinations that had
      # greater than MaxNPerSiteObsrvr counts available.  
      } else if (i > 1 & (j) == ifelse(InputData[i,]$NLists > MaxNPerSiteObsrvr, MaxNPerSiteObsrvr, InputData[i,]$NLists)){ 
         NewWideTable <- rbind(NewWideTable, NewWideRow)
         NewWideRow <- BlankRecord
      }   
      i <- i + InputData[i,]$NLists     #increment the counter to the start of the next location-observer block
   }
   rm(i, j, IsSiteName, NewWideRow)  #tidy up
   #
   #Any predictor variables that are represented by character data need to be converted into 
   # factors for analysis.  While it might seem logical (to me at least) to convert the 
   # character variables prior to assembling the unmarked wide-format data.frame, for some 
   # reason that I do not understand R seems to be back-converting factor data to character 
   # in some cases (I think it has to do with a column's data type being set in the first row, 
   # the presence of missing values not being factors, and causing newly appended rows to have 
   # their data coerced back to character from factors).  Anyway, as a result, we need to do a 
   # post-assembly coersion of all character data to R's factor data type.  Here we do this 
   # by identifying the first column of predictor data (the one immediately after the last "y.*" 
   # column), and then from there to the last column testing if each column has the data type 
   # "character" and if so converting it to the factor data type. Identify the column number 
   # for the first column of predictors, but only if there are any predictor variables.
   if (!is.null(SiteCovarSet) | !is.null(ObsCovarSet)){
      PredictCol1 <- (which(names(NewWideTable) == paste("y.", MaxNPerSiteObsrvr, sep="")) + 1)
      #Loop through the columns in turn, and change character data types to factors.
      for (i in PredictCol1:length(NewWideTable)){
         if (is.character(NewWideTable[,i]))
            NewWideTable[,i] <- as.factor(NewWideTable[,i])
      }
      rm(i)
   }
   #   
   #Finally, we can return the table.
   return(NewWideTable)
}







###########################################################################
###########################################################################
######   FUNCTION TO CONDUCT STEPWISE FORWARD SELECTION ON OCCU MODEL
###########################################################################
###########################################################################


occu_forward_aic <- function(data, all_state_covs, all_det_covs, 
                              start_state_covs="1", start_det_covs="1"){

   # create a data frame with one row per potential covariate
   cov_df <- data.frame(covs=c(all_state_covs, all_det_covs), param=c(rep("state", length(all_state_covs)), rep("det", length(all_det_covs))), in_model=0, aic=NA, delta_aic=NA)

   # specify which covariates are in the first model
   if(start_state_covs[1]!="1") cov_df$in_model[which(cov_df$covs %in% start_state_covs)] <- 1
   if(start_det_covs[1]!="1") cov_df$in_model[which(cov_df$covs %in% start_det_covs)] <- 1


   ##################################################################### 
   ###  create a model only the starting covariates

   model_so_far <- run_occu_mod(det_covs=start_det_covs, state_covs=start_state_covs, data=data)
   det_covariates <- start_det_covs
   state_covariates <- start_state_covs
   any_improve <- TRUE
   new_cov_count <- 0
   # create variables to fill with selected variables
   cov_name <- aic_improvement <- order <- vector()
   h <- 0

   while(any_improve){

      # keep track of how many covariates added (for assistance only)
      new_cov_count <- new_cov_count + 1
      print("#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #")
      print(paste("testing addition of covariate", new_cov_count))
      print(paste("current state covariates:", paste(state_covariates, collapse=", ")))
      print(paste("current det covariates:", paste(det_covariates, collapse=", ")))
      print("#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #")

      ##################################################################### 
      ###  test sequentially introducing more covariates to reduce AIC

      # find a df with all the covariates that are available to be added
      new_cov_df <- cov_df[cov_df$in_model==0,]

      # loop through the covariates that can be added
      for(iii in 1:nrow(new_cov_df)){

         print(paste("testing", new_cov_df$param[iii], "covariate", new_cov_df$covs[iii]))


         # define the new covariates as the covariates so far
         det_covs_new <- det_covariates
         state_covs_new <- state_covariates

         # add the new covariate to the sets
         if(new_cov_df$param[iii]=="state") state_covs_new <- c(state_covs_new, as.character(new_cov_df$covs[iii]))
         if(new_cov_df$param[iii]=="det") det_covs_new <- c(det_covs_new, as.character(new_cov_df$covs[iii]))

         # run an occu model with the new covariate
         add_1_cov <- run_occu_mod(det_covs=det_covs_new, state_covs=state_covs_new, data=data)

         # if any SEs are NA or error from model run
         if(inherits(add_1_cov, "try-error")) new_cov_df$aic[iii] <- NA
         if(!inherits(add_1_cov, "try-error")) {

            # check whether any of the SEs are NA
            all_se <- c(summary(add_1_cov)$state[,"SE"], summary(add_1_cov)$det[,"SE"])

            if(any(is.na(all_se))) {
               new_cov_df$aic[iii] <- NA
               } else {
               new_cov_df$aic[iii] <- add_1_cov@AIC
            }
         } # close if

      } # close iii

      # if all AIC are NA than set any_improve to FALSE (treat this as final model)
      if(all(is.na(new_cov_df$aic))) any_improve <- FALSE

      # find whether any of the new covariates reduced AIC
      # insert a catch so that if the original model didn't work, you can pick the lowest from any
      if(!inherits(model_so_far, "try-error")) if(any_improve) any_improve <- any((new_cov_df$aic[!is.na(new_cov_df$aic)]) < model_so_far@AIC)

      # find which covariate gave the largest improvement in AIC
      if(any_improve | inherits(model_so_far, "try-error")){

      # find which covariate gave the largest improvement in AIC
         if(!inherits(model_so_far, "try-error")) aic_diff <- model_so_far@AIC - new_cov_df$aic

         # if the model_so_far doesn't work, then pick the covariate that provides the lowest AIC
         if(inherits(model_so_far, "try-error")) aic_diff <- -new_cov_df$aic

         w <- which(aic_diff==max(aic_diff, na.rm=TRUE))
         new_cov <- as.character(new_cov_df$covs[w])

         w2 <- which(as.character(cov_df$covs)==new_cov)
         cov_df$in_model[w2] <- 1

         if(cov_df$param[w2]=="det") det_covariates <- c(det_covariates, new_cov)
         if(cov_df$param[w2]=="state") state_covariates <- c(state_covariates, new_cov)

         model_so_far <- run_occu_mod(det_covs=det_covariates, state_covs=state_covariates, data=data)

         print(paste("added covariate", new_cov))

         # save info about the new selected covariates
         cov_name <- c(cov_name, new_cov)
         aic_improvement <- c(aic_improvement, aic_diff[w])
         h <- h + 1
         order <- c(order, h)

      } # close if(any_improve)

   } # close while(any_improve)

   model_sel_record <- data.frame(order, cov_name, aic_improvement)

   # collect objects together to return
   ret <- list(model=model_so_far, state_covariates=state_covariates, det_covariates=det_covariates, add1_aic=cov_df, model_sel_record=model_sel_record)
   return(ret)

}



#######################################################################################
################ Wrapper function to run occu model by specifying covs ################
#######################################################################################

# This wrapper function takes vectors of the detectability covariates and the state
# covariates and constructs the formulae and runs the occupancy model. 
# It makes it easier to change the covariates and run the models, without having to 
# write out the formulae each time

run_occu_mod <- function(det_covs="1", state_covs="1", data){

   # create a catch if there are no det_covs or state_covs or they are NULL
   if(length(det_covs)==0 | is.null(det_covs)) det_covs <- "1"
   if(length(state_covs)==0 | is.null(state_covs)) state_covs <- "1"

   # construct detectability and state formulae
   det_form <- paste("~", paste(det_covs, collapse=" + "))
   state_form <- paste("~", paste(state_covs, collapse=" + "))

   # combine formulae together 
   form <- as.formula(paste(det_form, state_form, sep="  "))

   # calculate number of covariates
   no_covs_det <- length(strsplit(det_form, "+", fixed=TRUE)[[1]]) + ifelse(det_covs[1]=="1", 0, 1)
   no_covs_state <- length(strsplit(state_form, "+", fixed=TRUE)[[1]]) + ifelse(state_covs[1]=="1", 0, 1)
   no_covs <- no_covs_state + no_covs_det

   # run occu model
   occu_mod <- try(occu(form, 
                        starts=rep(0, no_covs),
                        data=data))

   return(occu_mod)

}



###########################################################################
###########################################################################
######  CREATE PREDICTION DATAFRAME FOR A FIXED DAY, BUT VARIABLE EXPERTISE
###########################################################################
###########################################################################

create_expertise_predict_df <- function(expertise, day=1, data){
   w <- which(data$DAY==day)

   df <- data.frame(log_pred=expertise)
   column_names <- vector()
   for(i in 1:ncol(data)){
      if(colnames(data)[i]=="log_pred") next
      df <- cbind(df, new1=as.vector(data[w,i]))
      column_names <- c(column_names, colnames(data)[i])
   }

   colnames(df)[2:ncol(df)] <- column_names

   return(df)
}


