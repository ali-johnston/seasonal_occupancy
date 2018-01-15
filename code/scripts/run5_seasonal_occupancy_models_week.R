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

specs_vec <- "Cardinalis_cardinalis"


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

day_smooth <- create_cyclic_bases(df=5, min_var=1, max_var=366, col_prefix="smd", var_name="DAY")


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

week_smooth <- create_cyclic_bases(df=5, min_var=1, max_var=52/no_weeks, col_prefix="smw", var_name="week_cat")



#####################################################################
###  merge with the checklist data and prediction dataframes

mod_data3 <- left_join(mod_data2, week_smooth, by="week_cat")

seas_pred_df3 <- left_join(seas_pred_df2, week_smooth, by="week_cat")



#####################################################################
###  tidy up workspace

rm(mod_data2, seas_pred_df2, dat_mod, seas_pred_df)



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

site_covar_names <- c(modis.umd.names, "LATITUDE", "LONGITUDE", "ELEV", paste0("week_", 1:52))
obs_covar_names <- c("TIME", "TIME2", "I.STATIONARY", "EFFORT_HRS", "EFFORT_DISTANCE_KM", 
   "NUMBER_OBSERVERS", "expertise", paste0("week_det_", 1:52))


# ##################################################################### 
# ###  remove rows with incomplete information for any relevant covariates

# cols_model <- mod_data3[,c(site_covar_names, obs_covar_names)]

# any_na <- apply(cols_model, 1, function(x){any(is.na(x))})

# mod_data4 <- mod_data3[!any_na,]
mod_data4 <- mod_data3

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
###### ADD DUMMY WEEK VARIABLES
###########################################################################
###########################################################################

add_dummy_week <- function(data, suffix=""){
   week_cat <- data$week_cat
   max_week <- max(week_cat)
   for(i in 1:max_week){
      week_dummy <- ifelse(week_cat==i, 1, 0)
      data <- cbind(data, week_dummy)
      colnames(data)[ncol(data)] <- paste0("week_", suffix, i)
   }
   return(data)
}


train_data2 <- add_dummy_week(train_data)
train_data3 <- add_dummy_week(train_data2, suffix="det_")

seas_pred_df4 <- add_dummy_week(seas_pred_df3)
seas_pred_df5 <- add_dummy_week(seas_pred_df4, suffix="det_")



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

site_covar_names_mod <- c(modis.umd.names, "LATITUDE", "LONGITUDE", paste(paste0("week_", 1:52), collapse=" + "))

# sig variables from full model
site_covar_names_mod <- c(modis.umd.names,  
   "LATITUDE", "LONGITUDE", paste(paste0("smw_", 1:5), collapse=" + "))

# most important variables from STEM models
site_covar_names_mod <- c("Deciduous_broad_PLAND", "Evergreen_broad_PLAND", 
   "ELEV", "LATITUDE", "LONGITUDE", paste(paste0("smw_", 1:5), collapse=" + "))

site_covar_names_mod <- c("Deciduous_broad_PLAND", "Evergreen_broad_PLAND", "ELEV",
   paste(paste0("week_", 1:52), collapse=" + "))


obs_covar_names_mod <- c("TIME", "TIME2", "I.STATIONARY", "EFFORT_HRS", "EFFORT_DISTANCE_KM",
   "NUMBER_OBSERVERS", "expertise", paste(paste0("week_det_", 1:52), collapse=" + "))



##################################################################### 
###  loop through species

for(sss in 1:length(specs_vec)){

   species <- specs_vec[sss]

   print(paste("Running for species:", species, "  species number run =", which(specs_vec==species)))
   spec.folder <- paste(result.folder, "Spec_", species, "_week_fac/" , sep="")
   dir.create(spec.folder)



   if(run_mods_new){
      train_data3$species_response <- train_data3[,c(species)]
      train_data3$species_response <- ifelse(train_data3$species_response==0, 0, 1)
      dir.create(spec.folder, recursive=T)


      ###### SUBSET BY WEEK

      no_obs_week_cat <- aggregate(train_data3$species_response, by=list(train_data3$week_cat), sum)
      colnames(no_obs_week_cat) <- c("week_cat", "no_obs")

      # find the first week_cat there are at least 20 observations of the species
      max_week_cat <- max(no_obs_week_cat$week_cat)
      sub <- filter(no_obs_week_cat, no_obs>29)
      first_week_cat1 <- min(sub$week_cat)
      last_week_cat1 <- max(sub$week_cat)

      # convert week categories to days of year
      first_day <- (first_week_cat1-1)*(7*no_weeks)+2
      last_day <- last_week_cat1*(7*no_weeks)-1
      if(last_week_cat1==max_week_cat) last_day <- 366

      train_data4 <- filter(train_data3, week_cat>=first_week_cat1, week_cat<=last_week_cat1)

      site_week_fac <- c(paste(paste0("week_", (first_week_cat1+1):last_week_cat1), collapse=" + "))

      # find only weeks that have at least 20 observations
      week_gt20 <- no_obs_week_cat$week_cat[no_obs_week_cat$no_obs>29]

      train_data4 <- filter(train_data3, week_cat %in% week_gt20)

      site_week_fac <- c(paste(paste0("week_", week_gt20[2:length(week_gt20)]), collapse=" + "))
      obs_week_fac <- c(paste(paste0("week_det_", week_gt20[2:length(week_gt20)]), collapse=" + "))

      ######################################################################################
      ### Turn into unmarked-friendly wide format

      print("Formatting data")

      train_data4$ID <- paste(train_data4$LOC_ID, train_data4$week_cat, sep="_")
      spec_data1 <- MakeUnmarkedWideInput(InputData = train_data4,
                                                SiteIDName = "ID", 
                                                ResponseName = "species_response", 
                                                SiteCovarSet = site_covar_names, 
                                                ObsCovarSet = obs_covar_names,
                                                MaxNPerSiteObsrvr = max_visits_site)
      # for debugging
      # InputData = train_data; SiteIDName = "ID"; ResponseName = "species_response"; SiteCovarSet = SiteCovarNames; ObsCovarSet = ObsCovarNames; MaxNPerSiteObsrvr = MaxNPerSiteObsrvr

      spec_data2 <- formatWide(dfin=spec_data1, type="unmarkedFrameOccu")

      write_loc1 <- paste0(spec.folder, "final_model_data_", Sys.Date(), ".txt")
      write.table(train_data4, write_loc1, row.names=FALSE)

      write_loc2 <- paste0(spec.folder, "final_model_data_unmarked_", Sys.Date(), ".rds")
      saveRDS(spec_data2, write_loc2)


      #####################################################################
      ###  find variables to use in model

      spec_var <- filter(spec_var_lookup, species_scientific==species, model==1)

      det_var <- filter(spec_var, var_type=="det")
      occ_var <- filter(spec_var, var_type=="occ")

      site_vars <- as.vector(as.matrix(occ_var$var))
      
      site_covar_names_mod <- c(site_vars, site_week_fac)

      det_vars <- as.vector(as.matrix(det_var$var))
      if(any(det_vars=="TIME")) {
         w <- which(det_vars=="TIME")
         det_vars[w] <- "TIME + TIME2"
      }

      obs_covar_names_mod <- c(det_vars, obs_week_fac)



      ###################################################################
      ### Fit Occupancy Models to eBird Data and Examine Model Output ###
      print("Fitting models")

      occu_mod <- run_occu_mod(det_covs=obs_covar_names_mod, state_covs=site_covar_names_mod, data=spec_data2)

      summary(occu_mod)

      # the predicted values are automatically calculated.
      saveRDS(occu_mod, paste0(spec.folder, "occu_mod.rds"))

   } # close if(run_mods_new)

   # load up the results of old models (if using)
   if(!run_mods_new){
      occu_mod <- readRDS(paste0(spec.folder, "occu_mod.rds"))
   }


   if(!inherits(occu_mod, "try-error")){

      print("Predicting from models")


      # predict season occupancy and detectability
      est_occ <- predict(object=occu_mod,
         newdata=seas_pred_df5, type="state")
      est_det <- predict(object=occu_mod,
         newdata=seas_pred_df5, type="det")

       #####################################################################
       ###  combine results for saving

       w <- which(seas_pred_df5$week_cat %in% week_gt20)

       res <- data.frame(DAY=seas_pred_df3$DAY[w], occ_est=est_occ$Predicted[w], occ_lcl=est_occ$lower[w], occ_ucl=est_occ$upper[w],
         det_est=est_det$Predicted[w], det_lcl=est_det$lower[w], det_ucl=est_det$upper[w])
      
      write_loc <- paste0(spec.folder, species, "_estimated_det_occ_by_DAY_weekfac.txt")
      write.table(res, write_loc, row.names=FALSE) 


      if(plot_graphs){


         print("Plotting graphs and maps")
         # PLOT THE ESTIMATED OCCUPANCY
         plot.loc <- paste(spec.folder, "Seasonal_occ_det_weekfac_restricted_time_", species, "_", Sys.Date(), ".png", sep="")
         png(plot.loc, width=14, height=8, units="cm", res=300, pointsize=9)
            par(mfrow=c(1,2), mar=c(5,1,1,1), oma=c(0, 4, 2, 0))
            plot(0, 0, xlim=c(0.5, 366.5), ylim=c(0, 1), col="white",
               xlab="Date", ylab="", type="l", lwd=2, xpd=NA, yaxt="n", xaxt="n")
            add_year_axis()
            axis(side=2, at=c(0, 0.5, 1))

            for(i in 1:length(week_gt20)){
               w <- which(seas_pred_df5$week_cat==week_gt20[i])
               sub_occ <- est_occ[w,]
               sub_pred <- seas_pred_df5[w,]
               n <- nrow(sub_occ)
               polygon(x=c(sub_pred$DAY[1]-0.5, sub_pred$DAY[n]+0.5, sub_pred$DAY[n]+0.5, sub_pred$DAY[1]-0.5, sub_pred$DAY[1]-0.5), 
                  y=c(sub_occ$lower[1], sub_occ$lower[1], sub_occ$upper[1], sub_occ$upper[1], sub_occ$lower[1]), col="grey", border="white")
               lines(c(sub_pred$DAY[1]-0.3, sub_pred$DAY[n]+0.3), c(sub_occ$Predicted[1], sub_occ$Predicted[1]), lty=1)
            }

            # add label
            text(x=0.02, y=0.95, labels="a", font=2, pos=4, cex=0.8)
            text(x=20, y=0.95, labels="Occupancy", pos=4, cex=0.8)

            plot(0, 0, xlim=c(0.5, 366.5), ylim=c(0, 1), col="white",
               xlab="Date", ylab="", type="l", lwd=2, xpd=NA, yaxt="n", xaxt="n")
            add_year_axis()

            for(i in 1:length(week_gt20)){
               w <- which(seas_pred_df5$week_cat==week_gt20[i])
               sub_det <- est_det[w,]
               sub_pred <- seas_pred_df5[w,]
               n <- nrow(sub_det)
               polygon(x=c(sub_pred$DAY[1]-0.5, sub_pred$DAY[n]+0.5, sub_pred$DAY[n]+0.5, sub_pred$DAY[1]-0.5, sub_pred$DAY[1]-0.5), 
                  y=c(sub_det$lower[1], sub_det$lower[1], sub_det$upper[1], sub_det$upper[1], sub_det$lower[1]), col="grey", border="white")
               lines(c(sub_pred$DAY[1]-0.3, sub_pred$DAY[n]+0.3), c(sub_det$Predicted[1], sub_det$Predicted[1]), lty=1)
            }

            # add label
            text(x=0.02, y=0.95, labels="b", font=2, pos=4, cex=0.8)
            text(x=20, y=0.95, labels="Detectability", pos=4, cex=0.8)

            text(x=-20, y=1.15, labels=gsub("_", " ", species), xpd=NA, font=3)
         dev.off()

      } # close if(plot_graphs)

   } # close if(inherits(try-error))

} # close species

