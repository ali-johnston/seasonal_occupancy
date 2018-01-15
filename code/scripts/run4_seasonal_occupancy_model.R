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
library(PresenceAbsence)


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
source("code/functions/functions_create_prediction_dataframe.R")



###########################################################################
###########################################################################
###### SET RUN PARAMETERS
###########################################################################
###########################################################################

#####################################################################
###  SET SPECIES TO RUN 

# wood thrush c("Hylocichla_mustelina")
specs_vec <- c("Hylocichla_mustelina", "Cardinalis_cardinalis", "Catharus_guttatus", "Columba_livia", 
               "Passerina_cyanea", "Pheucticus_ludovicianus", "Picoides_pubescens", "Poecile_atricapillus",
               "Seiurus_aurocapilla", "Sitta_carolinensis")


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
plot_graphs <- TRUE



#####################################################################
###  NUMBER OF WEEKS TO COMBINE IN A 'SITE'

no_weeks <- 1


#####################################################################
###  DEGREES OF FREEDOM FOR SMOOTHS

det_df <- 3
occ_df <- 3



###########################################################################
###########################################################################
######   READ AND MANIPULATE EBIRD DATA
###########################################################################
###########################################################################


#####################################################################
###  READ IN BIRD DATA

# Collect data for these species
dat_loc <- paste0(data.folder, "Checklist_by_grid.cell_week_max10vis.txt")
dat_mod <- read.table(dat_loc, header=TRUE)
dat_mod$LOC_ID <- dat_mod$grid.cell.number
dat_mod <- dat_mod[,-c(which(colnames(dat_mod)=="ID"))]



###########################################################################
###########################################################################
######   CREATE PREDICTION DATAFRAMES
###########################################################################
###########################################################################

srd_path <- paste(data.folder, "BCR23_SRD2016x1_srd.RData", sep="")


# Read in observer expertise scores:
expertise.folder <- paste0(result.folder, "expertise/")

# find most recent version of scores
expertise_scores <- read.table(paste0(expertise.folder, "observer_expertise.txt"), header=TRUE)
colnames(expertise_scores)[2] <- "expertise"

# define low, high and median expertise
low_expertise <- quantile(expertise_scores$expertise, 0.025)
median_expertise <- quantile(expertise_scores$expertise, 0.5)
high_expertise <- quantile(expertise_scores$expertise, 0.975)

# create seasonal prediction dataframe
seas_pred_df <- create_prediction_dataframe(srd_path=srd_path, 
                                             var_alter="DAY", min_var=1, max_var=366, by=1, 
                                             expertise=median_expertise)

colnames(seas_pred_df) <- update_modis_names(colnames(seas_pred_df))

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

day_smooth <- create_cyclic_bases(df=det_df, min_var=1, max_var=366, col_prefix="smd", var_name="DAY")

#####################################################################
###  merge with the checklist data and prediction dataframes

mod_data2 <- left_join(dat_mod, day_smooth, by="DAY")

seas_pred_df2 <- left_join(seas_pred_df, day_smooth, by="DAY")


#####################################################################
###  Add fortnight variable to checklist data and prediction dataframes

mod_data2$week_cat <- add_week_cat(mod_data2$DAY, nweek=no_weeks)

seas_pred_df2$week_cat <- add_week_cat(seas_pred_df2$DAY, nweek=no_weeks)



#####################################################################
###  create weekly smooth across fortnights (week_cat=1:26)

week_smooth <- create_cyclic_bases(df=det_df, min_var=1, max_var=52/no_weeks, col_prefix="smw", var_name="week_cat")



#####################################################################
###  merge with the checklist data and prediction dataframes

mod_data3 <- left_join(mod_data2, week_smooth, by="week_cat")

seas_pred_df3 <- left_join(seas_pred_df2, week_smooth, by="week_cat")



#####################################################################
###  tidy up workspace

rm(mod_data2, seas_pred_df2, dat_mod, seas_pred_df)


######################################################################
###   save prediction dataset to file to share with dave harris

seas_pred_df4 <- select(seas_pred_df3, -starts_with("smd_"), -starts_with("smw_")) 
write_loc <- paste0(data.folder, "prediction_dataframe.txt")
write.table(seas_pred_df4, write_loc, row.names=FALSE)




###########################################################################
###########################################################################
######    PREPARE DATA FOR MODELLING
###########################################################################
###########################################################################


##################################################################### 
###  variables for analysis (this should be exhaustive - include all vars possible for modelling)

modis.numbers <- c(0:10, 12, 13, 16)
modis.umd.names <- paste0("UMD_FS_C", modis.numbers, "_1500_PLAND")
modis.umd.names <- update_modis_names(modis.umd.names)

site_covar_names <- c(modis.umd.names, "LATITUDE", "LONGITUDE", "ELEV", paste0("smw_", 1:occ_df))
obs_covar_names <- c("TIME", "TIME2", "I.STATIONARY", "EFFORT_HRS", "EFFORT_DISTANCE_KM", 
   "NUMBER_OBSERVERS", "expertise", paste0("smd_", 1:det_df))


##################################################################### 
###  remove rows with incomplete information for any relevant covariates

cols_model <- mod_data3[,c(site_covar_names, obs_covar_names)]

any_na <- apply(cols_model, 1, function(x){any(is.na(x))})

mod_data4 <- mod_data3[!any_na,]


#####################################################################
###  tidy up workspace

rm(mod_data3)


##################################################################### 
###  split train and validation data

validation_data1 <- mod_data4[!mod_data4$train_loc,]

validation_data2 <- filter(mod_data4, YEAR==2016)

write.table(validation_data1, paste0(data.folder, "model_data_test1_10pct_", Sys.Date(), ".txt"), row.names=FALSE)
write.table(validation_data2, paste0(data.folder, "model_data_test2_2016", Sys.Date(), ".txt"), row.names=FALSE)


train_data <- mod_data4[mod_data4$train_loc,]

write.table(train_data, paste0(data.folder, "model_data_train_", Sys.Date(), ".txt"), row.names=FALSE)




###########################################################################
###########################################################################
###### READ IN LOOKUP WITH VARIABLES TO RUN FOR EACH SPECIES
###########################################################################
###########################################################################

spec_var_lookup <- read.csv("data/STEM_results/spec_top3_vars.csv")


###########################################################################
###########################################################################
###### RUN OCCUPANCY MODELS
###########################################################################
###########################################################################


##################################################################### 
###  variables to use for modelling

site_covar_names_mod <- c(modis.umd.names, "LATITUDE", "LONGITUDE", paste(paste0("smw_", 1:occ_df), collapse=" + "))

# sig variables from full model
site_covar_names_mod <- c(modis.umd.names,  
   "LATITUDE", "LONGITUDE", paste(paste0("smw_", 1:occ_df), collapse=" + "))

# most important variables from STEM models
site_covar_names_mod <- c("Deciduous_broad_PLAND", "Evergreen_broad_PLAND", 
   "ELEV", "LATITUDE", "LONGITUDE", paste(paste0("smw_", 1:occ_df), collapse=" + "))

site_covar_base <- c("ELEV", "LATITUDE", "LONGITUDE", paste(paste0("smw_", 1:occ_df), collapse=" + "))


site_smooth <- c(paste(paste0("smw_", 1:occ_df), collapse=" + "))


obs_covar_names_mod <- c("TIME", "TIME2", "I.STATIONARY", "EFFORT_HRS", "EFFORT_DISTANCE_KM",
   "NUMBER_OBSERVERS", "expertise", paste(paste0("smd_", 1:det_df), collapse=" + "))

obs_smooth <- c(paste(paste0("smd_", 1:det_df), collapse=" + "))

# obs_covar_names_mod <- c("EFFORT_HRS", "expertise", "EFFORT_DISTANCE_KM", paste(paste0("smd_", 1:3)))

# paste0("year", 2007:2011), 

# AIC_exp <- AIC_noexp <- AIC_exp_smdet <- vector()

final_ndf <- vector(length=length(specs_vec))

##################################################################### 
###  loop through species

for(sss in c(1, 5, 6, 9, 10)){ # 1:length(specs_vec)){

   species <- specs_vec[sss]

   print(paste("Running for species:", species, "  species number run =", which(specs_vec==species)))
   spec.folder <- paste(result.folder, "Spec_", species, "/" , sep="")

   if(run_mods_new){
      train_data$species_response <- train_data[,c(species)]
      train_data$species_response <- ifelse(train_data$species_response==0, 0, 1)
      dir.create(spec.folder, recursive=T)

      ######################################################################################
      ### Turn into unmarked-friendly wide format

      print("Formatting data")

      train_data$ID <- paste(train_data$LOC_ID, train_data$week_cat, sep="_")
      spec_data1 <- MakeUnmarkedWideInput(InputData = train_data,
                                                SiteIDName = "ID", 
                                                ResponseName = "species_response", 
                                                SiteCovarSet = site_covar_names, 
                                                ObsCovarSet = obs_covar_names,
                                                MaxNPerSiteObsrvr = max_visits_site)
      # for debugging
      # InputData = train_data; SiteIDName = "ID"; ResponseName = "species_response"; SiteCovarSet = SiteCovarNames; ObsCovarSet = ObsCovarNames; MaxNPerSiteObsrvr = MaxNPerSiteObsrvr

      spec_data2 <- formatWide(dfin=spec_data1, type="unmarkedFrameOccu")


      write_loc1 <- paste0(spec.folder, "final_model_data_", Sys.Date(), ".txt")
      write.table(train_data, write_loc1, row.names=FALSE)

      write_loc2 <- paste0(spec.folder, "final_model_data_unmarked_", Sys.Date(), ".rds")
      saveRDS(spec_data2, write_loc2)


      #####################################################################
      ###  find variables to use in model

      site_smooth <- c(paste(paste0("smw_", 1:occ_df), collapse=" + "))
      obs_smooth <- c(paste(paste0("smd_", 1:det_df), collapse=" + "))

      spec_var <- filter(spec_var_lookup, species_scientific==species, model==1)

      det_var <- filter(spec_var, var_type=="det")
      occ_var <- filter(spec_var, var_type=="occ")

      site_vars <- as.vector(as.matrix(occ_var$var))
      
      site_covar_names_mod <- c(site_vars, site_smooth)

      det_vars <- as.vector(as.matrix(det_var$var))
      if(any(det_vars=="TIME")) {
         w <- which(det_vars=="TIME")
         det_vars[w] <- "TIME + TIME2"
      }

      obs_covar_names_mod <- c(det_vars, obs_smooth)



      ###################################################################
      ### Fit Occupancy Models to eBird Data and Examine Model Output ###
      print("Fitting models")

      occu_mod <- run_occu_mod(det_covs=obs_covar_names_mod, state_covs=site_covar_names_mod, data=spec_data2)

      coefs <- summary(occu_mod)


      occ_se <- any(is.na(coefs$state[,c("SE")]))
      det_se <- any(is.na(coefs$det[,c("SE")]))

      # the predicted values are automatically calculated.
      saveRDS(occu_mod, paste0(spec.folder, "occu_mod.rds"))

   } # close if(run_mods_new)

   # load up the results of old models (if using)
   if(!run_mods_new){
      occu_mod <- readRDS(paste0(spec.folder, "occu_mod.rds"))
   }


   print("Predicting from models")


   # predict season occupancy and detectability
   est_occ <- predict(object=occu_mod,
      newdata=seas_pred_df3, type="state")
   est_det <- predict(object=occu_mod,
      newdata=seas_pred_df3, type="det")


   ##################################################################### 
   ###  VALIDATION DATASET 1; RANDOM 10% OF LOCATIONS NOT IN MODEL

   val1_occ <- predict(object=occu_mod, newdata=validation_data1, type="state")
   val1_det <- predict(object=occu_mod, newdata=validation_data1, type="det")

   # estimated recording species on list

   val1_obs <- val1_occ$Predicted * val1_det$Predicted

   val1_obs <- validation_data1[,c(species)]
   val1_obs <- ifelse(val1_obs==0, 0, 1)

   validation1 <- data.frame(ID=1:length(val1_obs), obs=val1_obs, pred=val1_occ$Predicted)


   ##################################################################### 
   ### VALIDATION DATASET 2; LOCATIONS WITH HIGH DETECTABILITY (IN BOTH TRAIN AND TEST)

   # find location x week combinations with high detectability across all visits (across all observers)
   # only 2012 data here, so no issues with multiple years
   validation_data2$val2_det <- predict(object=occu_mod, newdata=validation_data2, type="det")$Predicted
   validation_data2$val2_occ <- predict(object=occu_mod, newdata=validation_data2, type="state")$Predicted
   validation_data2$species_obs <- validation_data2[,c(species)]
   validation_data2$species_obs <- ifelse(validation_data2$species_obs==0, 0, 1)

   group_loc_fortnight <- group_by(validation_data2, grid.cell.number, week_cat)
   group_det <- function(x) { 1 - prod(1-x) }

   ID_det <- as.data.frame(summarise(group_loc_fortnight, 
        total_det = group_det(val2_det), # total detectability across all visits
        pred = mean(val2_occ), # average estimated occupancy (should be constant)
        obs = max(species_obs)))

   # reduce to only locations with at least 90% detectability
   validation_data2.1 <- filter(ID_det, total_det>0.9)

   # add ID column
   validation_data2.1$ID <- 1:nrow(validation_data2.1)
   validation2 <- select(validation_data2.1, ID, obs, pred)


   ##################################################################### 
   ###  AUC

   auc_ran10 <- presence.absence.accuracy(validation1)$AUC
   auc_det0.9 <- presence.absence.accuracy(validation2)$AUC



   ##################################################################### 
   ###  calculate brier scores on both validation sets
   brier <- function(pred, obs){ sum((pred-obs)^2)/length(pred) }

   brier_ran10 <- brier(pred=validation1$pred, obs=validation1$obs)

   brier_det0.9 <- brier(pred=validation2$pred, obs=validation2$obs)


   key_ppm <- data.frame(species, auc_ran10, auc_det0.9, brier_ran10, brier_det0.9, run_date=Sys.Date(), run_time=Sys.time(), week="smooth")

   # append to a file with all info combined
   key_ppm_loc <- paste0(result.folder, "all_species_key_ppm_season.txt")
   if(file.exists(key_ppm_loc)){
      write.table(key_ppm, key_ppm_loc, row.names=FALSE, append=TRUE)
    } else {
      write.table(key_ppm, key_ppm_loc, row.names=FALSE)
    }


    #####################################################################
    ###  combine results for saving

   res <- data.frame(DAY=seas_pred_df3$DAY, occ_est=est_occ$Predicted, occ_lcl=est_occ$lower, occ_ucl=est_occ$upper,
      det_est=est_det$Predicted, det_lcl=est_det$lower, det_ucl=est_det$upper)
   
   write_loc <- paste0(spec.folder, species, "_estimated_det_occ_by_DAY.txt")
   write.table(res, write_loc, row.names=FALSE) 


   if(plot_graphs){


      # find weeks when species is not present
      species_obs <- train_data[,c(species)]
      species_obs <- ifelse(species_obs==0, 0, 1)
      no_obs_week_cat <- aggregate(species_obs, by=list(train_data$week_cat), sum)
      colnames(no_obs_week_cat) <- c("week_cat", "no_obs")

      # find the first week_cat there is an observation
      max_week_cat <- max(no_obs_week_cat$week_cat)
      sub <- filter(no_obs_week_cat, no_obs>0)
      first_week_cat1 <- min(sub$week_cat)
      last_week_cat1 <- max(sub$week_cat)

      # convert week categories to days of year
      first_day <- (first_week_cat1-1)*(7*no_weeks)+2
      last_day <- last_week_cat1*(7*no_weeks)-1
      if(last_week_cat1==max_week_cat) last_day <- 366


      print("Plotting graphs and maps")
      # PLOT THE ESTIMATED OCCUPANCY

      plot.loc <- paste(spec.folder, "Seasonal_occ_det_restricted_time_", species, "_", Sys.Date(), ".png", sep="")
      png(plot.loc, width=14, height=8, units="cm", res=300, pointsize=9)
         par(mfrow=c(1,2), mar=c(5,1,1,1), oma=c(0, 4, 2, 0))
         plot(seas_pred_df3$DAY, est_occ$Predicted, xlim=c(0.5, 366.5), ylim=c(0, 1),
            xlab="Date", ylab="", type="l", lwd=2, xpd=NA, yaxt="n", xaxt="n")
         add_year_axis()
         axis(side=2, at=c(0, 0.5, 1))
         polygon(x=c(seas_pred_df3$DAY, rev(seas_pred_df3$DAY), seas_pred_df3$DAY[1]), 
            y=c(est_occ$lower, rev(est_occ$upper), est_occ$lower[1]), col="grey", border="white")
         lines(seas_pred_df3$DAY, est_occ$Predicted, lty=1)

         # cover up area with no observations (cannot separate two processes)
         polygon(x=c(-2, first_day-1, first_day-1, -2, -2), y=c(-0.01, -0.01, 1.01, 1.01, -0.01), col="white", border="white")
         polygon(x=c(368, last_day+1, last_day+1, 368, 368), y=c(-0.01, -0.01, 1.01, 1.01, -0.01), col="white", border="white")

         # add label
         text(x=0.02, y=0.95, labels="a", font=2, pos=4, cex=0.8)
         text(x=20, y=0.95, labels="Occupancy", pos=4, cex=0.8)

         plot(seas_pred_df3$DAY, est_det$Predicted, xlim=c(0.5, 366.5), ylim=c(0, 1),
            xlab="Date", ylab="", type="l", lwd=2, xpd=NA, yaxt="n", xaxt="n")
         add_year_axis()
         polygon(x=c(seas_pred_df3$DAY, rev(seas_pred_df3$DAY), seas_pred_df3$DAY[1]), 
            y=c(est_det$lower, rev(est_det$upper), est_det$lower[1]), col="grey", border="white")
         lines(seas_pred_df3$DAY, est_det$Predicted, lty=1)

         # cover up area with no observations (cannot separate two processes)
         polygon(x=c(-2, first_day-1, first_day-1, -2, -2), y=c(-0.01, -0.01, 1.01, 1.01, -0.01), col="white", border="white")
         polygon(x=c(368, last_day+1, last_day+1, 368, 368), y=c(-0.01, -0.01, 1.01, 1.01, -0.01), col="white", border="white")

         # add label
         text(x=0.02, y=0.95, labels="b", font=2, pos=4, cex=0.8)
         text(x=20, y=0.95, labels="Detectability", pos=4, cex=0.8)

         text(x=-20, y=1.15, labels=gsub("_", " ", species), xpd=NA, font=3)
      dev.off()

   } # close if(plot_graphs)

} # close species

