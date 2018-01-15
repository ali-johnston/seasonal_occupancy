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
library(rjags)
library(jagsUI)
library(coda)


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

det_df <- 5
occ_df <- 5


#####################################################################
###  USE SMALL EXTENT

small_extent <- TRUE



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
###  reduce spatial extent

if(small_extent){
   reduce_extent <- filter(mod_data4, LONGITUDE>-85) #-87.35)
   mod_data4 <- reduce_extent
}



#####################################################################
###  re-scale some variables

if(max(mod_data4$Deciduous_broad_PLAND, na.rm=TRUE)>1) mod_data4$Deciduous_broad_PLAND <- mod_data4$Deciduous_broad_PLAND/100
if(max(mod_data4$TIME2, na.rm=TRUE)>100) mod_data4$TIME2 <- mod_data4$TIME2/100





##################################################################### 
###  split train and validation data

validation_data1 <- mod_data4[!mod_data4$train_loc,]

validation_data2 <- filter(mod_data4, YEAR==2016)

write.table(validation_data1, paste0(data.folder, "model_data_test1_10pct_", ifelse(small_extent, "small_", ""), Sys.Date(), ".txt"), row.names=FALSE)
write.table(validation_data2, paste0(data.folder, "model_data_test2_2016", ifelse(small_extent, "small_", ""), Sys.Date(), ".txt"), row.names=FALSE)


train_data <- mod_data4[mod_data4$train_loc,]

write.table(train_data, paste0(data.folder, "model_data_train_", ifelse(small_extent, "small_", ""), Sys.Date(), ".txt"), row.names=FALSE)




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

# for(sss in 1:length(specs_vec)){
   sss <- 1

   species <- specs_vec[sss]

   print(paste("Running for species:", species, "  species number run =", which(specs_vec==species)))
   spec.folder <- paste(result.folder, "Spec_", species, "_jags/" , sep="")

   if(run_mods_new){
      train_data$species_response <- train_data[,c(species)]
      train_data$species_response <- ifelse(train_data$species_response==0, 0, 1)
      dir.create(spec.folder, recursive=T)

      ######################################################################################
      ### Turn into jags-friendly wide format

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

      #####################################################################
      ###  find variables to use in model


      site_smooth <- c(paste0("smw_", 1:occ_df))
      obs_smooth <- c(paste0("smd_", 1:det_df))

      spec_var <- filter(spec_var_lookup, species_scientific==species, model==1)

      det_var <- filter(spec_var, var_type=="det")
      occ_var <- filter(spec_var, var_type=="occ")

      site_vars <- as.vector(as.matrix(occ_var$var))
      
      site_covar_names_mod <- c(site_vars, site_smooth)

      det_vars <- as.vector(as.matrix(det_var$var))
      if(any(det_vars=="TIME") & !any(det_vars=="TIME2")) {
         det_vars <- c(det_vars, "TIME2")
      }

      obs_covar_names_mod <- c(det_vars, obs_smooth)


      #####################################################################
      ###  format data

      # species response matrix
      as.mat <- function(x){
         x <- as.matrix(x)
         rownames(x) <- NULL
         colnames(x) <- NULL
         return(x)
      }

      y <- as.mat(select(spec_data1, y.1:y.10))

      jags_data1 <- list(y=y)

      for(i in 1:length(site_covar_names_mod)){
         occ_var <- as.numeric(spec_data1[,site_covar_names_mod[i]])

         # add to jags_data list
         jags_data1[[length(jags_data1)+1]] <- occ_var
         names(jags_data1)[[length(jags_data1)]] <- site_covar_names_mod[i]
      }

      for(i in 1:length(obs_covar_names_mod)){
         covar <- obs_covar_names_mod[i]
         w1 <- which(colnames(spec_data1)==paste0(covar, ".1"))
         w2 <- which(colnames(spec_data1)==paste0(covar, ".10"))
         det_var <- as.mat(spec_data1[,c(w1:w2)])

         # add to jags_data list
         jags_data1[[length(jags_data1)+1]] <- det_var
         names(jags_data1)[[length(jags_data1)]] <- obs_covar_names_mod[i]
      }

      jags_data1[[length(jags_data1)+1]] <- nrow(y)
      names(jags_data1)[[length(jags_data1)]] <- "M"

      # find how many visits for each 'site'
      no_vis <- function(x) max(which(!is.na(x)))
      no_visits <- apply(y, 1, no_vis)

      jags_data1[[length(jags_data1)+1]] <- no_visits
      names(jags_data1)[[length(jags_data1)]] <- "J"
      

      write_loc <- paste0(spec.folder, "final_model_data_jags_", ifelse(small_extent, "small_", ""), Sys.Date(), ".rds")
      saveRDS(jags_data1, write_loc)


      #####################################################################
      ###  model structure
      
      # Specify model in BUGS language
      file_loc <- paste0(spec.folder, "jags_model.txt")
      sink(file_loc)
      cat("
      model {

      # Priors
      mean.p ~ dunif(0, 1)         # Detection intercept on prob. scale
      alpha0 <- logit(mean.p)      # Detection intercept
      alpha1 ~ dunif(-20, 20)      # Detection slope on EFFORT_HRS
      alpha2 ~ dunif(-20, 20)      # Detection slope on expertise
      alpha3 ~ dunif(-20, 20)      # Detection slope on TIME
      alpha4 ~ dunif(-20, 20)      # Detection slope on TIME2
      alpha5 ~ dunif(-20, 20)      # Detection slope on EFFORT_DISTANCE_KM
      alpha6 ~ dunif(-20, 20)      # Detection slope on smd_1
      alpha7 ~ dunif(-20, 20)      # Detection slope on smd_2
      alpha8 ~ dunif(-20, 20)      # Detection slope on smd_3
      alpha9 ~ dunif(-20, 20)      # Detection slope on smd_4
      alpha10 ~ dunif(-20, 20)     # Detection slope on smd_5

      mean.psi ~ dunif(0, 1)       # Occupancy intercept on prob. scale
      beta0 <- logit(mean.psi)     # Occupancy intercept
      beta1 ~ dunif(-20, 20)       # Occupancy slope on Deciduous_broad_PLAND
      beta2 ~ dunif(-20, 20)       # Occupancy slope on smw_1
      beta3 ~ dunif(-20, 20)       # Occupancy slope on smw_2
      beta4 ~ dunif(-20, 20)       # Occupancy slope on smw_3
      beta5 ~ dunif(-20, 20)       # Occupancy slope on smw_4
      beta6 ~ dunif(-20, 20)       # Occupancy slope on smw_5

      # Likelihood
      for (i in 1:M) {
         # True state model for the partially observed true state
         z[i] ~ dbern(psi[i])      # True occupancy z at site i
         logit(psi[i]) <- beta0 + (beta1 * Deciduous_broad_PLAND[i]) + (beta2 * smw_1[i]) + 
                         (beta3 * smw_2[i]) + (beta4 * smw_3[i]) + (beta5 * smw_4[i]) + 
                          (beta6 * smw_5[i])

         for (j in 1:J[i]) {
            # Observation model for the actual observations
            y[i,j] ~ dbern(p.eff[i,j])    # Detection-nondetection at i and j
            p.eff[i,j] <- z[i] * p[i,j]   # 'straw man' for WinBUGS
#             p[i,j] <- max(0.05, p[i,j]) # constrain detectability to be at least 0.05 to aid identifiability
            logit(p[i,j]) <- alpha0 + (alpha1 * EFFORT_HRS[i,j]) + 
                           (alpha2 * expertise[i,j]) + (alpha3 * TIME[i,j]) + 
                           (alpha4 * TIME2[i,j]) + (alpha5 * EFFORT_DISTANCE_KM[i,j]) + 
                           (alpha6 * smd_1[i,j]) + (alpha7 * smd_2[i,j]) + 
                           (alpha8 * smd_3[i,j]) + (alpha9 * smd_4[i,j]) + 
                           (alpha10 * smd_5[i,j])
         }
      }

      # Derived quantities
      # N.occ <- sum(z[])       # Number of occupied sites among sample of M
      # psi.fs <- N.occ/M       # Proportion of occupied sites among sample of M
      # for(k in 1:100){
      #    logit(psi.pred[k]) <- beta0 + beta1 * XvegHt[k] # psi predictions
      #    logit(p.pred[k]) <- alpha0 + alpha1 * Xwind[k]  # p predictions
      # }
      }
      ",fill = TRUE)
      sink()




      ###################################################################
      ### Fit Occupancy Models to eBird Data and Examine Model Output ###
      print("Fitting models")

      # Initial values: must give for same quantities as priors given !
      zst <- apply(y, 1, max, na.rm=TRUE)        # Avoid data/model/inits conflict
      inits <- function(){list(z = zst, mean.p = runif(1), alpha1 = runif(1), alpha2 = runif(1), 
         alpha3 = runif(1), alpha4 = runif(1), alpha5 = runif(1), alpha6 = runif(1), alpha7 = runif(1), 
         alpha8 = runif(1), alpha9 = runif(1), alpha10 = runif(1),  
         mean.psi = runif(1), beta1 = runif(1), beta2 = runif(1), beta3 = runif(1), beta4 = runif(1), 
         beta5 = runif(1), beta6 = runif(1))}

      # Parameters monitored
      params <- c(paste0("alpha", 0:10), paste0("beta", 0:6)) # Also estimate z = "conditional occ. prob."

      # MCMC settings
      na <- 10;   nc <- 3

      # run the occupancy model in jags
      set.seed(506)
      out <- jags.model(file_loc, data=jags_data1, inits=inits,
                      n.chains = nc, n.adapt=na)

      t3 <- Sys.time()      
      start_chains <- jags(data=jags_data1, inits=inits, parameters.to.save=params, 
         model.file=file_loc,
         n.chains=nc, n.adapt=na, n.iter=100, 
         n.burnin=0, n.thin=1,
         modules=c('glm'), parallel=FALSE, n.cores=NULL, DIC=TRUE, store.data=FALSE,
         seed=527, bugs.format=FALSE, verbose=TRUE)

      t4 <- Sys.time()

      start_chains2 <- jags(data=jags_data1, inits=inits, parameters.to.save=params, 
         model.file=file_loc,
         n.chains=nc, n.adapt=na, n.iter=100, 
         n.burnin=0, n.thin=1,
         modules=c('glm'), parallel=TRUE, n.cores=NULL, DIC=TRUE, store.data=FALSE,
         seed=527, bugs.format=FALSE, verbose=TRUE)

      t5 <- Sys.time()


      coefs <- coef(out)

      samps <- coda.samples(out, params, n.iter=1000)

      plot.mcmc(samps, trace=TRUE, density=TRUE, smooth=FALSE, ask=TRUE)

      samps3 <- coda.samples(out, params, n.iter=2000)

      plot(samps3, ask=TRUE)

      # mixing well
      # alpha 1, 5, 6, 
      # beta 1, 2 5, 6

      # save to file
      save_mcmc_loc <- paste0(spec.folder, "mcmc_n2000.rds")
      saveRDS(samps3, save_mcmc_loc)


      samps3a <- readRDS(save_mcmc_loc)
      samps4 <- coda.samples(samps3a, params, n.iter=1000)


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

