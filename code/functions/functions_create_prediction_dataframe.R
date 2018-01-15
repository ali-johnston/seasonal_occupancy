

###########################################################################
###########################################################################
###### FUNCTIONS TO CREATE AND MANIPULATE A VARIETY OF PREDICTION DATAFRAMES
###########################################################################
###########################################################################

# created by: ali johnston
# created on: 2017.06.06


#####################################################################
###  READ IN BASIC SPATIAL PREDICTION DATASET

spatial_srd_modis <- function(srd_path){
	load(srd_path)

	spatial_covs <- srd$X
	spatial_covs$LONGITUDE <- srd$locs$x
	spatial_covs$LATITUDE <- srd$locs$y

	rm(srd)

	# add in helpful column names for landcover
	modis.umd.names <- c("Water", "Evergreen_Needleleaf", "Evergreen_Broadleaf",
		"Deciduous_Needleleaf", "Deciduous_Broadleaf", "Mixed_Forest",
		"Woodland", "Wooden_Grassland", "Closed_Shrubland",
		"Open_Shrubland", "Grassland", "Cropland",
		"Urban_Built", "Bare")

	colnames(spatial_covs)[3:16] <- modis.umd.names

	return(spatial_covs)
}



#####################################################################
###  READ IN SPATIAL PREDICTION DATASET, BUT WITHOUT MODIS NAMES

spatial_srd <- function(srd_path){
	load(srd_path)

	spatial_covs <- srd$X
	spatial_covs$LONGITUDE <- srd$locs$x
	spatial_covs$LATITUDE <- srd$locs$y

	rm(srd)

	return(spatial_covs)
}


#####################################################################
###  ADD PREDICTION COVARIATES

add_standardised_covs <- function(srd_data, 
	YEAR=2012, TIME=7, EFFORT_HRS=1, EFFORT_DISTANCE_KM=1,
	NUMBER_OBSERVERS=1, I.STATIONARY=0, NLists=3, expertise=3){

	srd_data$YEAR <- YEAR
	srd_data$TIME <- TIME
	srd_data$TIME2 <- TIME^2
	srd_data$EFFORT_HRS <- EFFORT_HRS
	srd_data$EFFORT_DISTANCE_KM <- EFFORT_DISTANCE_KM
	srd_data$NUMBER_OBSERVERS <- NUMBER_OBSERVERS
	srd_data$I.STATIONARY <- I.STATIONARY
	srd_data$NLists <- NLists
	srd_data$expertise <- expertise

	return(srd_data)

}


#####################################################################
###  ADD DAY VARIABLE

	# srd_data$DAY <- DAY

	# ## Add day of year:
	# if(is.null(DAY)) srd_data$DAY <- as.numeric(format(strptime(paste(date, month, sep="/"), "%d/%m"), "%j"))
	# if(!is.null(DAY)) srd_data$DAY <- DAY





#####################################################################
###  CREATE AVERAGE OF SPATIAL VARIABLES

create_spatial_avg <- function(spatial_covs){
	spatial_covs_avg <- apply(spatial_covs, 2, mean)
	return(spatial_covs_avg)
}



#####################################################################
###  CREATE DATAFRAME VARYING IN A GIVEN VARIABLE

create_prediction_dataframe <- function(srd_path, var_alter="DAY", min_var, max_var, by, DAY=120, modis.names=FALSE, ...){


	r <- require(dplyr)

	# read in spatial data frame
	if(modis.names) spatial_covs <- spatial_srd_modis(srd_path=srd_path)
	if(!modis.names) spatial_covs <- spatial_srd(srd_path=srd_path)

	# calculate spatial averages
	spatial_avg <- create_spatial_avg(spatial_covs)
	spatial_avg_df <- as.data.frame(t(spatial_avg))

	# add in standardised checklist variables
	day_spatial_std <- add_standardised_covs(srd_data=spatial_avg_df, ...)

	# add in DAY variable
	day_spatial_std$DAY <- DAY

	# remove the varying variable
	w <- which(colnames(day_spatial_std)==var_alter)
	std_var_slim <- day_spatial_std[,-c(w)]

	# create a dataframe with a varying var_alter
	var <- seq(min_var, max_var, by=by)
	var_df <- data.frame(var=var)
	colnames(var_df) <- var_alter

	# join varying covariate with other standardised covariates
	var_df$join <- 1
	std_var_slim$join <- 1

	if(r) pred_df <- left_join(var_df, std_var_slim, by="join")
	if(!r) pred_df <- merge(var_df, std_var_slim, by="join", all.x=TRUE)

	# remove join
	w <- which(colnames(pred_df)=="join")
	pred_df_final <- pred_df[,-c(w)]

	return(pred_df_final)
}




