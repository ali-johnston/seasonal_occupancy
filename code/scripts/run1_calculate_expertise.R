###########################################################################
###########################################################################
###### CODE TO CALCULATE EXPERTISE SCORES OF EBIRD OBSERVERS
###########################################################################
###########################################################################

# created by: ali johnston


# load in required packages
library(dplyr)
library(mgcv)


# set parameters for whole script
data.tag <- "BCR23_2016"


# File directories in which to search for models and plot results:
data.root <- "data/"
data.folder <- paste(data.root, data.tag, "_data/", sep="")

result.root <- "results/"
result.folder <- paste(result.root, data.tag, "/expertise/", sep="")
dir.create(result.folder, recursive=T)

# set working directory as overall folder
setwd("/Users/ali/Documents/REPOS/seasonal_occupancy")


# Source in relevant functions:
source("code/functions/functions_create_prediction_dataframe.R")




###########################################################################
###########################################################################
######   READ IN ALL SPECIES AND CALCULATE NO SPECIES PER CHECKLIST
###########################################################################
###########################################################################


dat_loc <- paste0(data.folder, "BCR23_All.Spp_erd.RData")
load(dat_loc)

# Convert abundances to presences
species_matrix <- erd$y

species_presence <- ifelse(species_matrix==0, 0, 1)
rm(species_matrix)

# count up total species on each list
total_species <- apply(species_presence, 1, sum)
rm(species_presence)

# combine total species with other covariates
covariates <- cbind(erd$X, total_species)
covariates$OBSERVER_ID <- as.character(covariates$OBSERVER_ID)
locations <- erd$locs




###########################################################################
###########################################################################
###### FILTER CHECKLISTS TO KEEP ONLY THOSE THAT MEET GIVEN CRITERIA
###########################################################################
###########################################################################


# thresholds for max and min values of covariates
max_effort_hrs <- 5
max_distance_km <- 8.05
min_time <- 5
max_time <- 20 # (5.00-20.00 covers >98% of the dataset)
min_yr <- 2016
max_yr <- 2016
max_no_observers <- 10


dat0 <- filter(covariates, EFFORT_HRS<=max_effort_hrs, EFFORT_DISTANCE_KM<=max_distance_km, 
                           TIME<=max_time, TIME>=min_time, YEAR<=max_yr, YEAR>=min_yr,
                           NUMBER_OBSERVERS<=max_no_observers)





###########################################################################
###########################################################################
###### CALCULATE DERIVED VARIABLES
###########################################################################
###########################################################################


#####################################################################
###  add sqrt of effort_hrs
dat0$EFFORT_HRS_sqrt <- sqrt(dat0$EFFORT_HRS)


#####################################################################
###  calculate habitat diversity

# define function for GS
gini_simpson <- function(p) {1 - sum(p^2)}

# ensure that all proportions add up to 100 
lc_nos <- c(0:10, 12, 13, 16)
lc_names <- paste0("UMD_FS_C", lc_nos, "_1500_PLAND")
landcover_covs <- dat0[,c(lc_names)]

total_landcover <- apply(landcover_covs, 1, sum)
min_landcover <- min(total_landcover)

# minimum landcover is 100, so they are all ok. 

dat0$hab_div <- apply(landcover_covs/100, 1, gini_simpson)

# convert OBSERVER_ID to factor
dat0$OBSERVER_ID <- factor(dat0$OBSERVER_ID)




###########################################################################
###########################################################################
###### CREATE PREDICTION DATAFRAME OF STANDARDISED CHECKLISTS
###########################################################################
###########################################################################

srd_path <- paste0(data.folder, "BCR23_SRD2016x1_srd.RData")
pred_df <- create_prediction_dataframe(srd_path=srd_path, var_alter="EFFORT_HRS", 
   min_var=0.05, max_var=5, by=0.05, DAY=259, modis.names=FALSE)

pred_df$EFFORT_HRS_sqrt <- sqrt(pred_df$EFFORT_HRS)

pred_df$hab_div <- apply(pred_df[,lc_names]/100, 1, gini_simpson)

pred_1hr <- filter(pred_df, EFFORT_HRS==1)
pred_1hr$join <- 1





###########################################################################
###########################################################################
###### SUBSET THE OBSERVERS TO RUN IN CHUNKS
###########################################################################
###########################################################################



nrow(dat0)
# [1] 530368 # 2012-2016
# [1] 153404 # 2016 only


length(table(dat0$OBSERVER_ID))
# [1] 13216 # 2012-2016
# [1] 6206 # 2016 only



# Pick some random observers to include in all runs (for comparison)
all_obs <- names(table(dat0$OBSERVER_ID))
set.seed(181)
test_obs <- sample(all_obs, 20, replace=FALSE)

remaining_obs <- all_obs[-which(all_obs %in% test_obs)]

# define no of individuals in each group/batch and create an index for each observer
# defining their group/batch

# no_per_group <- 300 # 2 mins for all groups to run (with discrete=TRUE)
# no_per_group <- 1000 # 2 mins for all groups to run (with discrete=TRUE)
no_per_group <- 300
no_groups <- ceiling(length(remaining_obs) / (no_per_group-length(test_obs)))

# split remaining observers into n groups/batches
set.seed(187)
obs_cat <- sample(rep(1:no_groups, no_per_group), size=length(remaining_obs), replace=FALSE)

obs_remaining_df <- data.frame(OBSERVER_ID=remaining_obs, group=obs_cat)
obs_test_df <- data.frame(OBSERVER_ID=test_obs, group=0) # set group==0 for set of N individuals that are included in each batch

obs_df <- rbind(obs_remaining_df, obs_test_df) # combine the set of 20 observers and those that only appear in one group. 





###########################################################################
###########################################################################
###### RUN THE ESTIMATE OF EXPERTISE IN CHUNKS
###########################################################################
###########################################################################


t3 <- Sys.time()

for(i in 1:no_groups){

   # select only the observers in this group and those in the set that are repeated in every batch
   obs_group <- obs_df$OBSERVER_ID[obs_df$group==i | obs_df$group==0]
   sub_dat <- filter(dat0, OBSERVER_ID %in% obs_group)
   sub_dat$OBSERVER_ID <- sub_dat$OBSERVER_ID[drop=TRUE]


   # run the model to estimate observer expertise
   t1 <- Sys.time()
   m_group <- bam(total_species ~ EFFORT_HRS + EFFORT_HRS_sqrt 
                        + UMD_FS_C0_1500_PLAND + UMD_FS_C1_1500_PLAND + UMD_FS_C2_1500_PLAND 
                        + UMD_FS_C3_1500_PLAND + UMD_FS_C4_1500_PLAND + UMD_FS_C5_1500_PLAND
                        + UMD_FS_C6_1500_PLAND + UMD_FS_C7_1500_PLAND + UMD_FS_C8_1500_PLAND 
                        + UMD_FS_C9_1500_PLAND + UMD_FS_C10_1500_PLAND + UMD_FS_C12_1500_PLAND
                        + UMD_FS_C13_1500_PLAND + UMD_FS_C16_1500_PLAND + hab_div
                        + I.STATIONARY + EFFORT_DISTANCE_KM + NUMBER_OBSERVERS
                        + s(DAY, bs="cc", k=10) + s(TIME, k=5) 
                        + s(OBSERVER_ID, bs="re") 
                        + s(OBSERVER_ID, EFFORT_HRS, bs="re"),
                        data=sub_dat, family="poisson", discrete=TRUE)
   t2 <- Sys.time()

   print(t2 - t1)


   # create a prediction dataframe with one row for each OBSERVER_ID in sub_dat
   pred_observers <- data.frame(OBSERVER_ID=obs_group)   
   pred_observers$join <- 1
   pred_observers$OBSERVER_ID <- factor(as.character(pred_observers$OBSERVER_ID), levels=levels(sub_dat$OBSERVER_ID))

   # join with dataframe of standardised checklists (1hr, 1km, 7am, 1 observer, travelling, mid-september, etc.)
   pred_group <- left_join(pred_observers, pred_1hr, by="join")


   # predict from the m_group model
   p_group <- predict(m_group, newdata=pred_group, type="link", se.fit=TRUE)

   pred_obs <- data.frame(OBSERVER_ID=pred_group$OBSERVER_ID, fit=p_group$fit, se=p_group$se.fit)
   pred_obs$est <- exp(pred_obs$fit)
   pred_obs$lcl <- exp(pred_obs$fit - 1.96*pred_obs$se)
   pred_obs$ucl <- exp(pred_obs$fit + 1.96*pred_obs$se)

   pred_obs$batch <- i
   pred_obs$OBSERVER_ID <- as.character(pred_obs$OBSERVER_ID)

   if(i==1) all_pred <- pred_obs
   if(i>1) all_pred <- rbind(all_pred, pred_obs)

}


t4 <- Sys.time()




###########################################################################
###########################################################################
###### COMPARE THE ESTIMATES FOR 20 TEST OBSERVERS INCLUDED IN EACH BATCH
###########################################################################
###########################################################################



test_results <- filter(all_pred, as.character(OBSERVER_ID) %in% as.character(test_obs))
test_results$OBSERVER_ID <- test_results$OBSERVER_ID[drop=TRUE]


plot.loc <- paste0(result.folder, "compare_20_expertise_measures_gp", no_per_group, "_", Sys.Date(), ".png")
png(plot.loc, width=14, height=8, units="cm", pointsize=9, res=300)
   boxplot(log(est) ~ as.factor(OBSERVER_ID), data=test_results, main="group 300")
dev.off()


plot.loc <- paste0(result.folder, "compare_20_expertise_measures_gp", no_per_group, "_2runs_", Sys.Date(), ".png")
png(plot.loc, width=8, height=9, units="cm", pointsize=9, res=300)
   plot(test_results$fit[test_results$batch==1], test_results$fit[test_results$batch==2])
dev.off()

cor.test(test_results$fit[test_results$batch==1], test_results$fit[test_results$batch==2])


# calculate mean for each of the test variables

expertise_score <- select(all_pred, OBSERVER_ID, fit)
mean_expertise_score <- aggregate(expertise_score$fit, by=list(expertise_score$OBSERVER_ID), mean)
colnames(mean_expertise_score) <- c("OBSERVER_ID", "expertise_score")


write_loc <- paste0(result.folder, "observer_expertise.txt")
write.table(mean_expertise_score, write_loc, row.names=FALSE)







###########################################################################
###########################################################################
###### SENSE CHECK NOVICE AND EXPERT OBSERVER
###########################################################################
###########################################################################


# low_obs <- subset(dat0, as.character(OBSERVER_ID)=="251789")
# high_obs <- subset(dat0, as.character(OBSERVER_ID)=="340131")


# par(mfrow=c(1,2))
# plot(low_obs$EFFORT_HRS, low_obs$total_species, xlim=c(0, 5), ylim=c(0, 30))
# points(1, exp(mean_expertise_score$expertise_score[mean_expertise_score$OBSERVER_ID=="251789"]), pch=19, col="red")

# plot(high_obs$EFFORT_HRS, high_obs$total_species, xlim=c(0, 5), ylim=c(0, 30))
# points(1, exp(mean_expertise_score$expertise_score[mean_expertise_score$OBSERVER_ID=="340131"]), pch=19, col="red")

# test <- mean_expertise_score[with(mean_expertise_score, order(expertise_score)),]

# low_obs <- as.character(tail(test$OBSERVER_ID))
# high_obs <- as.character(head(test$OBSERVER_ID))

# par(mfrow=c(2,3))

# for(i in 1:3){
#    sub <- filter(dat0, as.character(OBSERVER_ID) %in% low_obs[i])

#    plot(sub$EFFORT_HRS, sub$total_species, xlim=c(0, 5), ylim=c(0, 100))
#    points(1, exp(tail(test$expertise_score)[i]), pch=19, col="red")

# }

# for(i in 1:3){
#    sub <- filter(dat0, as.character(OBSERVER_ID) %in% high_obs[i])

#    plot(sub$EFFORT_HRS, sub$total_species, xlim=c(0, 5), ylim=c(0, 100))
#    points(1, exp(head(test$expertise_score)[i]), pch=19, col="red")

# }



