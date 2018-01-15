

#####################################################################
###  FUNCTION TO COLLECT EBIRD CHECKLIST DATA AND EBIRD SPECIES DATA AND COMBINE
###  FOR ONLY DATA FOR SPECIFIED SPECIES

collate_bird_data <- function(checklist_loc, bird_loc, species, abund=FALSE){
	# INPUT
	# checklist_loc; path of folder and file name for checklist data
	# bird_loc; path of folder and file name for bird data
	# species; scientific name of species to include in output. Genus and species separated by "_". Genus has capital letter

	# OUTPUT
	# Checklist data combined with detection/non-detection or abundance data for each specified species

	# WARNINGS
	# assumes the bird data are organised with two leading columns: "index" and "type"
	# uses "index" column to match the checklist data and the bird data
	# index column is not regular eBird variable. It has been previously created. 

	r <- require(dplyr)

	all_checklist_covs <- read.csv(checklist_loc, h=T, sep=" ")

	bird_data <- read.csv(bird_loc, h=T, sep=" ")

	bird_data_sub <- bird_data[,c("index", "type", species)]
	rm(bird_data) # large file, so remove immediately

	# if detection/non-detection data, then convert abundance to 1/0 data
	if(!abund) for(i in 3:ncol(bird_data_sub)) bird_data_sub[,i] <- as.numeric(bird_data_sub[,i]>0)

	# if dplyr has loaded, the function inner_join is more efficient, so use that. 
	if(r) combine_data <- inner_join(all_checklist_covs, bird_data_sub, by=c("index", "type"))
	if(!r) combine_data <- merge(all_checklist_covs, bird_data_sub, by=c("index", "type"))

	# tidy up workspace
	rm(all_checklist_covs, bird_data_sub)

	# return the combined dataframe
	return(combine_data)
}



#####################################################################
#### FUNCTION TO ADD A COLOUR RAMP SCALE TO A MAP:

addscale <- function(x, y, zlim, at, cols, cex=0.5, zbreaks, zlabels=NULL){
	# INPUT
	# x; vector length 2 with the limits of the coloured box (order by x[2]>x[1])
	# y; vector length 2 with the limits of the coloured box (order by y[2]>y[1])
	# zlim; vector length 2 with the limits of the scale of the response (order by z[2]>z[1])
	# at; locations of labels (on the z scale)
	# cols; vector with all colours for the ramp. More colours gives smoother transition
	# cex; size of the text for the labels
	# zbreaks; where to change colours (on the z scale). particularly useful if you don't want evenly spaced colour bands
	# 			should be one shorter than cols vector
	# zlabels; optional character vector. If true, provides labels for the break points. Defaults to the locations of breaks (at)
	# 			particularly useful to provide if you want lower SF in the labels than the actual break locations
	# 			should be same length as 'at' vector

	# OUTPUT
	# Adds a colour ramp to the plot device


	y_col <- y[1] + zbreaks/max(zbreaks)*(y[2]-y[1])
	for(i in 1:length(cols)) polygon(x=c(x[1], x[1], x[2], x[2], x[1]), y=c(y_col[i], y_col[i+1], y_col[i+1], y_col[i], y_col[i]), col=cols[i], border=alpha("black", 0))

	# Add tick marks
	y_ticks <- y[1] + (y[2]-y[1])*(1-(zlim[2]-at)/(zlim[2]-zlim[1]))
	x_loc <- x[2]+0.3*(x[2]-x[1])
	segments(x0=x[2], x1=x_loc, y0=y_ticks, y1=y_ticks)

	# Add labels
	if(is.null(zlabels)) zlabels <- at
	text(x=x_loc, y=y_ticks, labels=zlabels, cex=cex, pos=4)
}



#####################################################################
###  FUNCTION TO FIND THE MOST RECENT VERSION OF A FILE

###  Useful when processed files are stamped in their file name with Sys.Date()

find_most_recent <- function(folder_name, file_name){
	# INPUT
	# folder_name; name of the folder to look in
	# file_name; generic file name to search for (with no date)

	# OUTPUT
	# specific file name of the most recent version\

	# WARNING
	# assumes there are no hyphens in the file_name

	# check for presence of Hmisc package
	check_package <- require(Hmisc)
	if(!check_package) print("Install package Hmisc to run this function")

	# list all the files in the folder
	all_files <- list.files(folder_name)

	# find file names that match the string
	match_files <- all_files[grep(file_name, all_files)]

	# set up vectors to keep the names and dates
	name_no_date <- date <- vector(length=length(match_files))
	date <- as.character(date)

	# find the date and separate from file name	
	for(i in 1:length(match_files)){

		s <- substring.location(match_files[i], "-")$first
		# add in code here to deal with hyphens in file name

		if(s[1]==0) {
			# if no date, set to 1st jan 1990
			date[i] <- "1990-01-01"
			name_no_date[i] <- substr(match_files[i], 1, nchar(match_files[i])-4)
		}

		if(length(s)==2) {
			date[i] <- as.character(as.POSIXct(substr(match_files[i], s[1]-4, s[2]+2)))
			name_no_date[i] <- substr(match_files[i], 1, s[1]-5)
		}


	}

	# find the latest date
	date <- as.POSIXct(date)
	w <- which(date==max(date))

	# return the file name with the latest date
	use_file_name <- match_files[w]

	separate <- ""
	if(substr(folder_name, nchar(folder_name), nchar(folder_name))!="/") separate="/"
	use_file_path <- paste0(folder_name, use_file_name, sep=separate)

	return(c(use_file_path, use_file_name))
}




#####################################################################
###  FUNCTION TO ADD DUMMY VARIABLES FOR YEAR

add_year_dummy <- function(data, yr_seq){
	# INPUT
	# data; object name of dataframe
	# yr_seq; sequence of years (numeric vector) for which to create dummy variables

	# OUTPUT
	# dataframe with column added for each year in yr_seq
	# 			variable names are "year2007" for 2007

	# WARNING
	# assumes data has a numeric column called YEAR with numeric 4-figure years

   for(y in 1:length(yr_seq)){
      data$yeardummy <- ifelse(data$YEAR==yr_seq[y], 1, 0)
      colnames(data)[ncol(data)] <- paste0("year", yr_seq[y])
   }
   return(data)
}


#####################################################################
###  HELPER FUNCTIONS TO CREATE SMOOTH BASIS FUNCTIONS


# function to create smooth 

create_cyclic_bases <- function(df=10, min_var, max_var, col_prefix="sm", var_name="DAY"){
	# INPUT
	# df; number of degrees of freedom for the smooth
	# min_var; the minimum value for the smoothed variable (for DAY usually 1)
	# max_var; the maximum value for the smoothed variable (for DAY usually 366)
	# col_prefix; the prefix for the column names for the bases

	# OUTPUT
	# dataframe with the basis functions for the smooth

	# WARNING
	# the cyclic construction is sensitive to the maximum and minimum values

	require(mgcv)

	sm_var <- seq(min_var, max_var, by=1)
	data <- data.frame(sm_var=sm_var)

	temp <- smooth.construct(
		object = s(sm_var, k=df+1, bs="cc", fx=T), 
		data=data,
		knots = NULL)
	new_vars <- temp$X
   	colnames(new_vars) <- paste(col_prefix, 1:df, sep="_")

   	data_return <- cbind(data, new_vars)
   	colnames(data_return)[1] <- var_name

   	return(data_return)
}



#####################################################################
###  FUNCTION TO ADD FORTNIGHT CATEGORY

add_week_cat <- function(day, nweek=2){

	ndays <- nweek*7

	week_cat <- ceiling(day/ndays)

	max_weeks <- 52/nweek
	week_cat <- ifelse(week_cat>max_weeks, max_weeks, week_cat)

	return(week_cat)
}



#####################################################################
###  FUNCTION TO ADD A PRETTY YEAR AXIS TO AN ANNUAL PLOT

add_year_axis <- function(label.size=0.7, week=FALSE){
	month1 <- as.POSIXlt(paste0("2016-",1:12, "-02"))
	day1 <- c(month1$yday, 366)

	if(week) day1 <- day1/7

	month15 <- as.POSIXlt(paste0("2015-", 1:12, "-16"))
	day15 <- month15$yday

	if(week) day15 <- day15/7

	axis(side=1, at=day1, labels=rep("", 13))
	axis(side=1, tick=FALSE, at=day15, padj=-1, cex.axis=label.size,
	labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
}



#####################################################################
###  CHANGE NAMES OF MODIS LANDCOVER VARIABLES

update_modis_names <- function(dnames){

	dnames <- sub( "UMD_FS_C0_1500", "Water", dnames)
	dnames <- sub( "UMD_FS_C1_1500", "Evergreen_needle", dnames)
	dnames <- sub( "UMD_FS_C2_1500", "Evergreen_broad", dnames)
	dnames <- sub( "UMD_FS_C3_1500", "Deciduous_needle", dnames)		
	dnames <- sub( "UMD_FS_C4_1500", "Deciduous_broad", dnames)
	dnames <- sub( "UMD_FS_C5_1500", "Mixed_forest", dnames)
	dnames <- sub( "UMD_FS_C6_1500", "Closed_shrubland", dnames)
	dnames <- sub( "UMD_FS_C7_1500", "Open_shrubland", dnames)
	dnames <- sub( "UMD_FS_C8_1500", "Woody_savannas", dnames)		
	dnames <- sub( "UMD_FS_C9_1500", "Savannas", dnames)	
	dnames <- sub( "UMD_FS_C10_1500", "Grasslands", dnames)
	dnames <- sub( "UMD_FS_C12_1500", "Croplands", dnames)
	dnames <- sub( "UMD_FS_C13_1500", "Urban_Built", dnames)
	dnames <- sub( "UMD_FS_C16_1500", "Barren", dnames)	

	return(dnames)
}


