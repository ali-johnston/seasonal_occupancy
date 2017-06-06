###########################################################################
###########################################################################
###### CODE TO RUN SEASONAL OCCUPANCY MODELS ON EBIRD DATA 
###########################################################################
###########################################################################

# created by: ali johnston


 # load in required packages
library(unmarked)
library(scales)
library(dplyr)


# set parameters for whole script
data.tag <- "BCR23"


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
source("code/functions/functions_create_prediction_dataframe.R")



###########################################################################
###########################################################################
###### SET RUN PARAMETERS
###########################################################################
###########################################################################

#####################################################################
###  SET SPECIES TO RUN 

specs_vec <- c("Hylocichla_mustelina")


#####################################################################
###  YEARS FOR ANALYSIS

min_yr <- 2004
max_yr <- 2012



#####################################################################
###  SET PLOT PARAMETERS
# plot_wd <- 9
# plot_ht <- 9
# plot_res <- 600


#####################################################################
###  OCCUPANCY MODELLING PARAMS

# define minimum and maximum visits used per site. 
# if below minimum, site is not included
# if above maximum, a random selection of visits up to max is used
min_visits_site <- 2
max_visits_site <- 10



#####################################################################
###  FLAGS FOR WHICH PARTS OF THE CODE TO RUN

run_mods_new <- TRUE



###########################################################################
###########################################################################
######   READ AND MANIPULATE EBIRD DATA
###########################################################################
###########################################################################


##################################################################### 
###  READ IN EXPERTISE DATA

# Read in observer expertise scores:
expertise.folder <- paste0(data.folder, "expertise/")

# find most recent version of scores
expertise_recent <- find_most_recent(folder_name=expertise.folder, file_name="observer_expertise")
expertise_scores <- read.csv(expertise_recent[1], h=T, sep=" ")



#####################################################################
###  READ IN BIRD DATA

# Collect data for these species
checklist_path <- paste0(data.folder, data.tag, "_all_checklist_data.csv")
bird_path <- paste0(data.folder, data.tag, "_all_bird_data.csv", sep="")

# combine together checklist info and bird data for selected species
all_bird_lists <- collate_bird_data(checklist_loc=checklist_path, bird_loc=bird_path, species=specs_vec, abund=FALSE)



#####################################################################
###  SELECT ONLY RELEVANT YEARS

bird_lists <- filter(all_bird_lists, YEAR<=max_yr, YEAR>=min_yr)

# tidy up workspace and remove large file 
rm(all_bird_lists)



#####################################################################
###  MERGE EXPERTISE SCORES WITH CHECKLIST INFO

# ensure both OBSERVER_IDs are in the same format
bird_lists$OBSERVER_ID <- as.character(bird_lists$OBSERVER_ID)
if(substr(expertise_scores$OBSERVER_ID[1], 1, 4)=="1/1/") expertise_scores$OBSERVER_ID <- substr(as.character(expertise_scores$OBSERVER_ID), 5, nchar(as.character(expertise_scores$OBSERVER_ID)))

# tidy up expertise scores
expertise_scores_merge <- expertise_scores[,c("OBSERVER_ID", "log_pred")]
colnames(expertise_scores_merge)[which(colnames(expertise_scores_merge)=="log_pred")] <- "expertise"

# merge 2 dataframes together by OBSERVER_ID
allspec_mod_data <- left_join(bird_lists, expertise_scores_merge, by="OBSERVER_ID")



#####################################################################
###  ADD QUADRATIC FOR THE TIME OF DAY
allspec_mod_data$TIME2 <- allspec_mod_data$TIME^2


#####################################################################
###  tidy up workspace

rm(bird_lists)


###########################################################################
###########################################################################
######   CREATE PREDICTION DATAFRAMES
###########################################################################
###########################################################################

srd_path <- paste(data.folder, "BCR23.srd.3km.data.RData", sep="")

# define low, high and median expertise
low_expertise <- quantile(expertise_scores$log_pred, 0.025)
median_expertise <- quantile(expertise_scores$log_pred, 0.5)
high_expertise <- quantile(expertise_scores$log_pred, 0.975)


# create seasonal prediction dataframe
seas_pred_df <- create_prediction_dataframe(srd_path=srd_path, 
                                             var_alter="DAY", min_var=1, max_var=366, by=1, 
                                             expertise=median_expertise)

# # create expertise dataframe
# expertise_pred_df <- create_prediction_dataframe(srd_path=srd_path, 
#                                                 var_alter="expertise", min_var=low_expertise, max_var=high_expertise, by=(high_expertise-low_expertise)/100,
#                                                 DAY=90)




###########################################################################
###########################################################################
###### ADD SEASONAL SMOOTH BASES TO MODEL AND PREDICTION DATAFRAMES
###########################################################################
###########################################################################

#####################################################################
###  create daily smooth across year (DAY=1:366)

day_smooth <- create_cyclic_bases(df=5, min_var=1, max_var=366, col_prefix="smd", var_name="DAY")


#####################################################################
###  merge with the checklist data and prediction dataframes

allspec_mod_data2 <- left_join(allspec_mod_data, day_smooth, by="DAY")

seas_pred_df2 <- left_join(seas_pred_df, day_smooth, by="DAY")


#####################################################################
###  Add fortnight variable to checklist data and prediction dataframes

allspec_mod_data2$week_cat <- add_week_cat(allspec_mod_data2$DAY)

seas_pred_df2$week_cat <- add_week_cat(seas_pred_df2$DAY)



#####################################################################
###  create weekly smooth across fortnights (week_cat=1:26)

fortnight_smooth <- create_cyclic_bases(df=5, min_var=1, max_var=26, col_prefix="smw", var_name="week_cat")



#####################################################################
###  merge with the checklist data and prediction dataframes

allspec_mod_data3 <- left_join(allspec_mod_data2, fortnight_smooth, by="week_cat")

seas_pred_df3 <- left_join(seas_pred_df2, fortnight_smooth, by="week_cat")



#####################################################################
###  tidy up workspace

rm(allspec_mod_data2, seas_pred_df2, allspec_mod_data, seas_pred_df)



###########################################################################
###########################################################################
######    PREPARE DATA FOR MODELLING
###########################################################################
###########################################################################


##################################################################### 
###  variables for analysis (this should be exhaustive - include all vars possible for modelling)

modis.umd.names <- c("Water", "Evergreen_Needleleaf", "Evergreen_Broadleaf",
   "Deciduous_Needleleaf", "Deciduous_Broadleaf", "Mixed_Forest",
   "Woodland", "Wooden_Grassland", "Closed_Shrubland",
   "Open_Shrubland", "Grassland", "Cropland",
   "Urban_Built", "Bare")

site_covar_names <- c(modis.umd.names, "LATITUDE", "LONGITUDE", paste0("smw_", 1:5))
obs_covar_names <- c("TIME", "TIME2", "I.STATIONARY", "EFFORT_HRS", "EFFORT_DISTANCE_KM", 
   "NUMBER_OBSERVERS", "expertise", paste0("smd_", 1:5))


##################################################################### 
###  remove rows with incomplete information for any relevant covariates

cols_model <- allspec_mod_data3[,c(site_covar_names, obs_covar_names)]

any_na <- apply(cols_model, 1, function(x){any(is.na(x))})

allspec_mod_data4 <- allspec_mod_data3[!any_na,]


#####################################################################
###  tidy up workspace

rm(allspec_mod_data3)



######################################################################################
### Filter Site-Observer Combinations with Too Few or Too Many Observation Records ###

allspec_mod_data4$LOC_ID <- paste(allspec_mod_data4$LONGITUDE, allspec_mod_data4$LATITUDE, sep="_")

#  RetainRangeCounts function randomly selects up to MaxNPerSiteObsrvr visits
# for each site. set.seed ensures the random selection is repeatable. 
set.seed(192)
mod_data <- RetainRangeCountsPerLocObsvrWeek(InputData = allspec_mod_data4, 
                                                      MinNPerSiteObsrvr = min_visits_site, 
                                                      MaxNPerSiteObsrvr = max_visits_site,
                                                      NWeek=2,
                                                      repeat_visit_vars=c("LOC_ID")) # "OBSERVER_ID", "YEAR", ...



##################################################################### 
###  randomly select 10% of the LOCATIONS for test data

all_loc <- table(allspec_mod_data4$LOC_ID)
set.seed(259)
test_loc <- names(all_loc)[sample(1:length(all_loc), round(length(all_loc)*0.1), replace=FALSE)]

mod_data$train_loc <- ifelse(mod_data$LOC_ID %in% test_loc, FALSE, TRUE)

train_data <- mod_data[mod_data$train_loc,]

write.table(train_data, paste0(data.folder, "model_data_train_", Sys.Date(), ".txt"), row.names=FALSE)
write.table(mod_data, paste0(data.folder, "model_data_all_", Sys.Date(), ".txt"), row.names=FALSE)


##################################################################### 
###  validation data

validation_data <- mod_data[!mod_data$train_loc,]

validation_data2 <- filter(mod_data, YEAR==2012)

write.table(validation_data, paste0(data.folder, "model_data_test1_10pct_", Sys.Date(), ".txt"), row.names=FALSE)
write.table(validation_data, paste0(data.folder, "model_data_test2_2012_", Sys.Date(), ".txt"), row.names=FALSE)



##################################################################### 
###  variables to use for modelling

SiteCovarNamesMod <- c(modis.umd.names, "LATITUDE", "LONGITUDE", paste0("year", 2007:2011), paste(paste0("smw_", 1:5), collapse=" + "))
ObsCovarNamesMod <- c("TIME", "TIME2", "I.STATIONARY", "EFFORT_HRS", "EFFORT_DISTANCE_KM", 
   "NUMBER_OBSERVERS", "log_pred")



AIC_exp <- AIC_noexp <- AIC_exp_smdet <- vector()

##################################################################### 
###  loop through species

# for(species in specs_vec){
   run_spec <- 1

   species <- specs_vec[run_spec]
   print(paste("Running for species:", species, "  species number run =", which(specs_vec==species)))
   spec.folder <- paste(result.folder, "Spec_", species, "/" , sep="")


   if(run_mods_new){
      train_data$species_response <- train_data[,c(species)]
      dir.create(spec.folder, recursive=T)

      ######################################################################################
      ### Turn into unmarked-friendly wide format

      print("Formatting data")

      spec_data1 <- MakeUnmarkedWideInput(InputData = train_data3,
                                                SiteIDName = "ID", 
                                                ResponseName = "species_response", 
                                                SiteCovarSet = SiteCovarNames, 
                                                ObsCovarSet = ObsCovarNames,
                                                MaxNPerSiteObsrvr = MaxNPerSiteObsrvr)
      # for debugging
      # InputData = train_data2; SiteIDName = "ID"; ResponseName = "species_response"; SiteCovarSet = SiteCovarNames; ObsCovarSet = ObsCovarNames; MaxNPerSiteObsrvr = MaxNPerSiteObsrvr

      spec_data2 <- formatWide(dfin=spec_data1, type="unmarkedFrameOccu")

      ###################################################################
      ### Fit Occupancy Models to eBird Data and Examine Model Output ###
      print("Fitting models")

      expertise_covariate_selection <- occu_forward_aic(data=spec_data2, 
         all_state_covs=SiteCovarNamesMod, all_det_covs=ObsCovarNamesMod,
         start_state_covs=SiteCovarNamesStart, 
         start_det_covs=ObsCovarNamesStart)
      mod_expertise <- expertise_covariate_selection$model

      print("Finished expertise model selection")
      print("Final model is:")

      summary(mod_expertise)

      det_cov_noexp <- expertise_covariate_selection$det_covariates
      if(any(det_cov_noexp=="log_pred")) det_cov_noexp <- det_cov_noexp[-which(det_cov_noexp=="log_pred")]
      smooth_int <- paste(paste0("smi_", 1:5), collapse=" + ")      
      if(any(det_cov_noexp==smooth_int)) det_cov_noexp <- det_cov_noexp[-which(det_cov_noexp==smooth_int)]

      # add expertise covariate to model
      mod_no_expertise <- run_occu_mod(
         det_covs=det_cov_noexp, 
         state_covs=expertise_covariate_selection$state_covariates, 
         data=spec_data2)

      # summarise the three models
      summary(mod_no_expertise)
      summary(mod_expertise)

      AIC_noexp <- c(AIC_noexp, AIC(mod_no_expertise))
      AIC_exp <- c(AIC_exp, AIC(mod_expertise))

#      gof_exp <- parboot(mod_expertise)
#      gof_noexp <- parboot(mod_no_expertise)


      # the predicted values are automatically calculated.
      saveRDS(mod_no_expertise, paste0(spec.folder, "mod_no_expertise.rds"))
      saveRDS(mod_expertise, paste0(spec.folder, "mod_expertise.rds"))
 
      write_modsel_noexp <- paste0(spec.folder, "modsel_noexpertise.txt")
      write.table(noexpertise_covariate_selection$model_sel_record, write_modsel_noexp, row.names=FALSE)

   } # close if(run_mods_new)

   # load up the results of old models (if using)
   if(!run_mods_new){
      mod_no_expertise <- readRDS(paste0(spec.folder, "mod_no_expertise.rds"))
      mod_expertise <- readRDS(paste0(spec.folder, "mod_expertise.rds"))
   }


   print("Predicting from models")


   # predict season occupancy
   occ_no_expertise <- predict(object=mod_no_expertise,
      newdata=srd_year_med, type="state")
   occ_expertise <- predict(object=mod_expertise,
      newdata=srd_year_med, type="state")

   # predict season detectability
   det_no_expertise <- predict(object=mod_no_expertise,
      newdata=srd_year_med, type="det")
   det_expertise <- predict(object=mod_expertise,
      newdata=srd_year_med, type="det")

   # check for estimates of SE that are NA


   ##################################################################### 
   ###  VALIDATION DATASET 1 RANDOM 10%

   val1_occ_expertise <- predict(object=mod_expertise, newdata=validation_data1, type="state")
   val1_occ_no_expertise <- predict(object=mod_no_expertise, newdata=validation_data1, type="state")

   val1_det_expertise <- predict(object=mod_expertise, newdata=validation_data1, type="det")
   val1_det_no_expertise <- predict(object=mod_no_expertise, newdata=validation_data1, type="det")

   # estimated recording species on list

   val1_obs_expertise <- val1_occ_expertise$Predicted * val1_det_expertise$Predicted
   val1_obs_no_expertise <- val1_occ_no_expertise$Predicted * val1_det_no_expertise$Predicted

   val1_obs <- validation_data1[,c(species)]

   validation1 <- data.frame(ID=1:length(val1_obs), obs=val1_obs, pred_no_exp=val1_obs_no_expertise, pred_exp=val1_obs_expertise)

   ##################################################################### 
   ### VALIDATION DATASET 2 HIGH DETECTABILITY

   # find location x week combinations with high detectability across all visits (across all observers)
   # only 2012 data here, so no issues with multiple years
   validation_data2$val2_det_expertise <- predict(object=mod_expertise, newdata=validation_data2, type="det")$Predicted
   validation_data2$val2_occ_no_expertise <- predict(object=mod_no_expertise, newdata=validation_data2, type="state")$Predicted
   validation_data2$val2_occ_expertise <- predict(object=mod_expertise, newdata=validation_data2, type="state")$Predicted
   validation_data2$species_obs <- validation_data2[,c(species)]

   group_loc_fortnight <- group_by(validation_data2, LONGITUDE, LATITUDE, week_cat)
   group_det <- function(x) { 1 - prod(1-x) }

   ID_det <- as.data.frame(summarise(group_loc_fortnight, 
        total_det_exp = group_det(val2_det_expertise), # total detectability across all visits
        pred_no_exp = mean(val2_occ_no_expertise), # average estimated occupancy (should be constant)
        pred_exp = mean(val2_occ_expertise),
        obs = max(species_obs)))

   # reduce to only locations with at least 90% detectability
   validation_data2.1 <- filter(ID_det, total_det_exp>0.9)

   # add ID column
   validation_data2.1$ID <- 1:nrow(validation_data2.1)
   validation2 <- select(validation_data2.1, ID, obs, pred_no_exp, pred_exp)


   ##################################################################### 
   ###  AIC
   model <- c("noexp", "exp")
   aic <- c(AIC(mod_no_expertise), AIC(mod_expertise))


   ##################################################################### 
   ###  AUC

   auc_ran10 <- presence.absence.accuracy(validation1)$AUC
   auc_det0.9 <- presence.absence.accuracy(validation2)$AUC



   ##################################################################### 
   ###  calculate brier scores on both validation sets
   brier <- function(pred, obs){ sum((pred-obs)^2)/length(pred) }

   brier_noexp_ran10 <- brier(pred=validation1$pred_no_exp, obs=validation1$obs)
   brier_exp_ran10 <- brier(pred=validation1$pred_exp, obs=validation1$obs)

   brier_noexp_det0.9 <- brier(pred=validation2$pred_no_exp, obs=validation2$obs)
   brier_exp_det0.9 <- brier(pred=validation2$pred_exp, obs=validation2$obs)


   brier_ran10 <- c(brier_noexp_ran10, brier_exp_ran10)
   brier_det0.9 <- c(brier_noexp_det0.9, brier_exp_det0.9)


   key_ppm <- data.frame(species, model, aic, auc_ran10, auc_det0.9, brier_ran10, brier_det0.9, season="may", run_date=Sys.Date(), run_time=Sys.time())

   # append to a file with all info combined
   key_ppm_loc <- paste0(res.folder, "all_species_key_ppm_season.txt")
   if(file.exists(key_ppm_loc)){
      write.table(key_ppm, key_ppm_loc, row.names=FALSE, append=TRUE)
    } else {
      write.table(key_ppm, key_ppm_loc, row.names=FALSE)
    }


   if(plot_graphs){

      par(mfrow=c(2,2))
      plot(srd_year_med$DAY, occ_expertise$Predicted, ylim=c(0,1), type="l")
      plot(srd_year_med$DAY, det_expertise$Predicted, type="l", ylim=c(0,1))
      plot(srd_year_med$DAY, occ_no_expertise$Predicted, ylim=c(0,1), type="l")
      plot(srd_year_med$DAY, det_no_expertise$Predicted, type="l", ylim=c(0,1))


      # predict season-expertise interactions
      det_expertise_0.5 <- predict(object=mod_expertise,
         newdata=srd_year_med, type="det")
      det_no_expertise_0.5 <- predict(object=mod_no_expertise,
         newdata=srd_year_med, type="det")
      det_expertise_0.025 <- predict(object=mod_expertise,
         newdata=srd_year_novice, type="det")
      det_expertise_0.975 <- predict(object=mod_expertise,
         newdata=srd_year_expert, type="det")

      det_exp_day1 <- predict(object=mod_expertise, 
         newdata=expertise_day1, type="det")
      det_exp_day90 <- predict(object=mod_expertise, 
         newdata=expertise_day90, type="det")
      det_exp_day180 <- predict(object=mod_expertise, 
         newdata=expertise_day180, type="det")


      add_year_axis <- function(label.size=0.7){
         month1 <- as.POSIXlt(paste0("2016-",1:12, "-02"))
         day1 <- c(month1$yday, 366)

         month15 <- as.POSIXlt(paste0("2015-", 1:12, "-16"))
         day15 <- month15$yday

         axis(side=1, at=day1, labels=rep("", 13))
         axis(side=1, tick=FALSE, at=day15, padj=-1, cex.axis=label.size,
            labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
      }


      print("Plotting graphs and maps")
      # PLOT THE COMPARISON IN PREDICTIONS BETWEEN THE TWO MODELS
      plot.loc <- paste(spec.folder, "Seasonal_occ_se_srd_", Sys.Date(), ".png", sep="")
      png(plot.loc, width=14, height=8, units="cm", res=300, pointsize=9)
         par(mfrow=c(1,2), mar=c(5,1,1,1), oma=c(0, 4, 0, 0))
         plot(srd_year_med$DAY, occ_no_expertise$Predicted, xlim=c(0.5, 366.5), ylim=c(0, 1),
            xlab="Date", ylab="Estimated occupancy", type="l", lwd=2, xpd=NA, yaxt="n", xaxt="n")
         add_year_axis()
         axis(side=2, at=c(0, 0.5, 1))
         lines(srd_year_med$DAY, occ_no_expertise$lower, lty=2)
         lines(srd_year_med$DAY, occ_no_expertise$upper, lty=2)
         text(x=0.02, y=0.95, labels="a", font=2, pos=4, cex=0.8)
         text(x=20, y=0.95, labels="No expertise", pos=4, cex=0.8)

         plot(srd_year_med$DAY, occ_expertise$Predicted, xlim=c(0.5, 366.5), ylim=c(0, 1),
            xlab="Date", ylab="", type="l", lwd=2, yaxt="n", xaxt="n")
         add_year_axis()
         lines(srd_year_med$DAY, occ_expertise$lower, lty=2)
         lines(srd_year_med$DAY, occ_expertise$upper, lty=2)
         text(x=0.02, y=0.95, labels="b", font=2, pos=4, cex=0.8)
         text(x=20, y=0.95, labels="Expertise interaction", pos=4, cex=0.8)

         text(x=-10, y=1.3, labels=gsub("_", " ", species), xpd=NA, font=3, pos=4)
      dev.off()



      ylimits <- c(0, 1)
      ylimits[1] <- floor(min(c(occ_no_expertise$Predicted, occ_expertise$Predicted))*10)/10
      ylimits[2] <- ceiling(max(c(occ_no_expertise$Predicted, occ_expertise$Predicted))*10)/10
      y_mid <- sum(ylimits)/2
      plot.loc <- paste(spec.folder, "Fortnight_occ_compare_", Sys.Date(), ".png", sep="")
      png(plot.loc, width=7, height=8, units="cm", res=300, pointsize=9)
         plot(occ_no_expertise$Predicted, occ_expertise$Predicted, xlim=ylimits, ylim=ylimits,
            xlab="Estimated occupancy - no expertise", ylab="Estimated occupancy - expertise",
            pch=16, cex=0.8, xpd=NA, yaxt="n", xaxt="n")
         axis(side=2, at=c(ylimits[1], y_mid, ylimits[2]))
         axis(side=1, at=c(ylimits[1], y_mid, ylimits[2]))
         abline(0, 1, col="grey")
         segments(x0=occ_no_expertise$lower, x1=occ_no_expertise$upper, y0=occ_expertise$Predicted, y1=occ_expertise$Predicted)
         segments(x0=occ_no_expertise$Predicted, x1=occ_no_expertise$Predicted, y0=occ_expertise$lower, y1=occ_expertise$upper)
      dev.off()


      # PLOT THE COMPARISON IN DETECTABILITY BETWEEN THE TWO MODELS
      plot.loc <- paste(spec.folder, "Seasonal_det_se_srd_", Sys.Date(), ".png", sep="")
      png(plot.loc, width=14, height=8, units="cm", res=300, pointsize=9)
         par(mfrow=c(1,2), mar=c(5,1,1,1), oma=c(0, 4, 0, 0))
         plot(srd_year_med$DAY, det_no_expertise_0.5$Predicted, xlim=c(0.5, 366.5), ylim=c(0, 1),
            xlab="Date", ylab="Estimated detectability", type="l", lwd=2, xpd=NA, yaxt="n", xaxt="n")
         axis(side=2, at=c(0, 0.5, 1))
         add_year_axis()
         lines(srd_year_med$DAY, det_no_expertise_0.5$lower, lty=2)
         lines(srd_year_med$DAY, det_no_expertise_0.5$upper, lty=2)
         text(x=0.02, y=0.95, labels="a", font=2, pos=4, cex=0.8)
         text(x=20, y=0.95, labels="No expertise", pos=4, cex=0.8)

         plot(srd_year_med$DAY, det_expertise_0.5$Predicted, xlim=c(0.5, 366.5), ylim=c(0, 1),
            xlab="Date", ylab="Estimated detectability", type="l", lwd=2, yaxt="n", xaxt="n")
         add_year_axis()
         lines(srd_year_med$DAY, det_expertise_0.5$lower, lty=2)
         lines(srd_year_med$DAY, det_expertise_0.5$upper, lty=2)
         text(x=0.02, y=0.95, labels="b", font=2, pos=4, cex=0.8)
         text(x=20, y=0.95, labels="Expertise interaction", pos=4, cex=0.8)

         text(x=-2, y=1.3, labels=gsub("_", " ", species), xpd=NA, font=3, pos=4)
      dev.off()



      dup <- duplicated(srd_year_med$week_cat)
      w <- which(!dup)
      cols <- c("grey75", "grey50", "black")

      limits <- range(c(det_no_expertise_0.5$Predicted, det_expertise_0.5$Predicted, det_expertise_0.025$Predicted, det_expertise_0.975$Predicted))
      limits[1] <- floor(limits[1]*5)/5
      limits[2] <- ceiling(limits[2]*5)/5
      limits_mid <- mean(limits)
      plot.loc <- paste(spec.folder, "Fortnight_det_compare_", Sys.Date(), ".png", sep="")
      png(plot.loc, width=7, height=8, units="cm", res=300, pointsize=9)
         plot(det_no_expertise_0.5$Predicted[w], det_expertise_0.5$Predicted[w], xlim=limits, ylim=limits,
            xlab="Estimated detectability - no expertise", ylab="Estimated detectability - expertise",
            pch=16, cex=0.8, xpd=NA, col=cols[2], xaxt="n", yaxt="n")
         axis(side=1, at=c(limits[1], limits_mid, limits[2]))
         axis(side=2, at=c(limits[1], limits_mid, limits[2]))
         abline(0, 1, col="grey")
         segments(x0=det_no_expertise_0.5$Predicted[w], x1=det_no_expertise_0.5$Predicted[w],
            y0=det_expertise_0.025$Predicted[w], y1=det_expertise_0.975$Predicted[w], col=cols[2])
         points(det_no_expertise_0.5$Predicted[w], det_expertise_0.025$Predicted[w], pch=16, cex=0.8, col=cols[1])
         points(det_no_expertise_0.5$Predicted[w], det_expertise_0.975$Predicted[w], pch=16, cex=0.8, col=cols[3])

         x_pos <- limits[1]
         if(species=="Setophaga_americana") x_pos <- 0.4
         legend(x=x_pos, y=limits[2], pch=16, col=rev(cols), legend=c("97.5th quantile expertise", "50th quantile expertise", "2.5th quantile expertise"), cex=0.6)
      dev.off()



      # PLOT THE COMPARISON IN PREDICTIONS BETWEEN THE TWO MODELS
      plot.loc <- paste(spec.folder, "Seasonal_det_exp_nov_srd_", Sys.Date(), ".png", sep="")
      png(plot.loc, width=14, height=8, units="cm", res=300, pointsize=9)
         par(mfrow=c(1,2), mar=c(5,1,1,1), oma=c(0, 4, 0, 0))
         plot(srd_year_med$DAY, det_no_expertise_0.5$Predicted, xlim=c(0.5, 366.5), ylim=c(0, 1),
            xlab="Date", ylab="Estimated detectability", type="l", lwd=2, xpd=NA, xaxt="n", yaxt="n")
         add_year_axis()
         axis(side=2, at=c(0, 0.5, 1))
         lines(srd_year_med$DAY, det_no_expertise_0.5$lower, lty=2)
         lines(srd_year_med$DAY, det_no_expertise_0.5$upper, lty=2)
         text(x=0.02, y=0.95, labels="a", font=2, pos=4, cex=0.8)
         text(x=20, y=0.95, labels="No expertise", pos=4, cex=0.8)

         plot(srd_year_med$DAY, det_expertise_0.5$Predicted, xlim=c(0.5, 366.5), ylim=c(0, 1),
            xlab="Date", ylab="Estimated detectability", type="l", lwd=2, xaxt="n", yaxt="n", col=cols[2])
         add_year_axis()
         axis(side=2, at=c(0, 0.5, 1))
         lines(srd_year_med$DAY, det_expertise_0.5$lower, lty=2, col=cols[2])
         lines(srd_year_med$DAY, det_expertise_0.5$upper, lty=2, col=cols[2])

         lines(srd_year_med$DAY, det_expertise_0.025$Predicted, lwd=2, col=cols[1])
         lines(srd_year_med$DAY, det_expertise_0.025$lower, lty=2, col=cols[1])
         lines(srd_year_med$DAY, det_expertise_0.025$upper, lty=2, col=cols[1])

         lines(srd_year_med$DAY, det_expertise_0.975$Predicted, lwd=2, col=cols[3])
         lines(srd_year_med$DAY, det_expertise_0.975$lower, lty=2, col=cols[3])
         lines(srd_year_med$DAY, det_expertise_0.975$upper, lty=2, col=cols[3])

         text(x=0.02, y=0.95, labels="b", font=2, pos=4, cex=0.8)
         text(x=20, y=0.95, labels="Expertise interaction", pos=4, cex=0.8)

         text(x=-10, y=1.3, labels=gsub("_", " ", species), xpd=NA, font=3, pos=4)
      dev.off()


      # PLOT THE COMPARISON IN PREDICTIONS BETWEEN THE TWO MODELS
      ymax <- ceiling(max(c(det_exp_day1$upper, det_exp_day90$upper, det_exp_day180$upper))/5)*5
      cols <- c("black", "grey50", "grey75")
      plot.loc <- paste(spec.folder, "Seasonal_det_expertise_", Sys.Date(), ".png", sep="")
      png(plot.loc, width=7, height=8, units="cm", res=300, pointsize=9)
         par(mar=c(5,4,4,2))
         plot(expertise_day1$log_pred, det_exp_day180$Predicted, ylim=c(0, 1),
            xlab="Observer expertise score", ylab="Estimated detectability", type="l", lwd=2,
            main=gsub("_", " ", species), col=cols[3], xpd=NA, yaxt="n")
         axis(side=2, at=c(0, 0.5, 1))
         lines(expertise_day1$log_pred, det_exp_day180$lower, lty=2, col=cols[3])
         lines(expertise_day1$log_pred, det_exp_day180$upper, lty=2, col=cols[3])

         lines(expertise_day90$log_pred, det_exp_day90$Predicted, lwd=2, col=cols[2])
         lines(expertise_day90$log_pred, det_exp_day90$lower, lty=2, col=cols[2])
         lines(expertise_day90$log_pred, det_exp_day90$upper, lty=2, col=cols[2])

         lines(expertise_day180$log_pred, det_exp_day1$Predicted, lwd=2, col=cols[1])
         lines(expertise_day180$log_pred, det_exp_day1$lower, lty=2, col=cols[1])
         lines(expertise_day180$log_pred, det_exp_day1$upper, lty=2, col=cols[1])

      dev.off()

   } # close if(plot_graphs)





   # ##################################################################### 
   # ###   PREDICT EFFECT OF NUMBER OF HOURS
   # hrs <- seq(0, 2, by=1/60)
   # n <- length(hrs)
   # nd_hrs <- data.frame(EFFORT_HRS=hrs, TIME=rep(srd_data$TIME[1], n), DAY=rep(srd_data$DAY[1], n), I.STATIONARY=rep(srd_data$I.STATIONARY[1], n),
   #                      NUMBER_OBSERVERS=rep(srd_data$NUMBER_OBSERVERS[1], n), EFFORT_DISTANCE_KM=rep(srd_data$EFFORT_DISTANCE_KM[1], n),
   #                      NLists=rep(srd_data$NLists[1], n), log_pred=rep(srd_data$log_pred[1], n))
   # predict_hrs <- predict(object = mod_expertise, newdata = nd_hrs, type = "det")


   # ##################################################################### 
   # ###   PREDICT EFFECT OF EXPERTISE
   # q1 <- quantile(expertise_scores$log_pred, 0.025)
   # q2 <- quantile(expertise_scores$log_pred, 0.975)
   # nd_expertise <- data.frame(TIME=srd_data$TIME[1], DAY=srd_data$DAY[1], I.STATIONARY=srd_data$I.STATIONARY[1],
   #                      NUMBER_OBSERVERS=srd_data$NUMBER_OBSERVERS[1], EFFORT_DISTANCE_KM=srd_data$EFFORT_DISTANCE_KM[1],
   #                      EFFORT_HRS=srd_data$EFFORT_HRS[1], NLists=srd_data$NLists[1], log_pred=seq(q1, q2, length.out=100))
   # predict_expertise <- predict(object = mod_expertise, newdata = nd_expertise, type = "det")

   # max_y <- ceiling(max(c(predict_hrs$upper, predict_expertise$upper))*5)/5

   # ##################################################################### 
   # ###   PLOT EFFECT OF NUMBER OF HOURS
   # for(k in 1:2){
   #    ylimits <- c(0, 1)
   #    if(k==2) ylimits[2] <- max_y
   #    plot.loc <- paste(spec.folder, "Det_effort_hrs_", model.month, "_scale", k,  Sys.Date(), ".png", sep="")
   #    png(plot.loc, width=8, height=8, units="cm", res=300, pointsize=9)
   #       plot(x = range(nd_hrs$EFFORT_HRS), 
   #            y = ylimits, 
   #            xlab = "Length of checklist (hrs)",
   #            ylab = "Probability of Detection",
   #            type = "n",
   #            xaxt="n", yaxt="n")
   #       axis(side=1, at=c(0, 1, 2))
   #       axis(side=2, at=c(0, ylimits[2]/2, ylimits[2]))
   #       lines(nd_hrs$EFFORT_HRS, predict_hrs$Predicted, lwd = 2, col = "black")
   #       lines(nd_hrs$EFFORT_HRS, predict_hrs$lower, lwd = 1, col = "black", lty = 2)
   #       lines(nd_hrs$EFFORT_HRS, predict_hrs$upper, lwd = 1, col = "black", lty = 2)
   #    dev.off()
   # } # close k


   # ##################################################################### 
   # ###   PLOT EFFECT OF EXPERTISE
   # for(k in 1:2){
   #    ylimits <- c(0, 1)
   #    if(k==2) ylimits[2] <- max_y
   #    plot.loc <- paste(spec.folder, "Det_expertise_", model.month, "_scale",k, "_", Sys.Date(), ".png", sep="")
   #    png(plot.loc, width=8, height=8, units="cm", res=300, pointsize=9)
   #       plot(x = range(nd_expertise$log_pred), 
   #            y = ylimits, 
   #            xlab = "Expertise score",
   #            ylab = "Probability of Detection",
   #            type = "n", yaxt="n", xaxt="n")
   #       axis(side=2, at=c(0, ylimits[2]/2, ylimits[2]))
   #       axis(side=1, at=c(2.5, 3, 3.5))
   #       lines(nd_expertise$log_pred, predict_expertise$Predicted, lwd = 2, col = "black")
   #       lines(nd_expertise$log_pred, predict_expertise$lower, lwd = 1, col = "black", lty = 2)
   #       lines(nd_expertise$log_pred, predict_expertise$upper, lwd = 1, col = "black", lty = 2)
   #    dev.off()
   # }

# } # close sss species

