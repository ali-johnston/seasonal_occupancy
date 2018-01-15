

# ----------------------------------------------------------------------
# The ERD data object is a list and contains the following objects:
#  $ locs : (n x 3) contains longitude & latitude & date of search
#  		where date of serach recorded as "fractional year" 
#  $ X    : (n x k) design matrix with k predictors 
#  $ y    : (n x p) responses as maximum count of p species on search. 
#  $ species : (p x 2) data frame with common & scientific names of 
# 		all p species
# 	NOTE y = -1 ==> Checklist count entered as "X"
#
#  $ spatial.extent.list = list defining spatial extent
#  $ spatial.extent.index = logical vector indicating all rows 
# 		in ERD2016 bigmem file that are within spatial.extent.list
#  $ temporal.extent.list = list defining temporal extent
#  $ temporal.extent.index = logical vector indicating all rows 
# 		in ERD2016 bigmem file that are within temporal.extent.list
# 
# Each row of each object contains inforamtion on 
# the n checklists, all in the the same row order. 
#
# The SRD data object is also a list with analogous structure. 
# The SRD object differs from the ERD object in that 
# there are no observational variables and there is no temporal dimension.   
#
# ----------------------------------------------------------------------
# Directories 
# ----------------------------------------------------------------------	
	parent.dir <- "/Users/df36/projects/2017_CONF_EURING/intra/"
	source.dir <- "/Users/df36/projects/2017_CONF_EURING/intra/source/"
	data.dir <- paste(parent.dir, "data.erd2016/", sep="") 
# ----------------------------------------------------------------------
# Libraries
# ----------------------------------------------------------------------	
	source(paste(source.dir, "grid_sample_functions.R", sep=""))
	library(data.table)
	library(plyr)
# ----------------------------------------------------------------------
# LOAD ERD, SRD, & Expertise 
# ----------------------------------------------------------------------	
	#output.tag <- "TCBB_CA"
	output.tag <- "BCR23"
	load(file = paste(data.dir, output.tag,"_erd.RData", sep=""))
	load(file = paste(data.dir, output.tag,"_SRD2016x1_srd.RData", sep = ""))
	# -------------------------------------------------------------
	# UMD LANDCOVER 500 meter resolution landcover classification 
	# using the UMD classi- fication scheme 
	# https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mcd12q1
	# -------------------------------------------------------------	
	dnames <- names(erd$X)
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
	names(erd$X) <- dnames
	names(erd$X)
	# ---------------
	dnames <- names(srd$X)
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
	names(srd$X) <- dnames
	names(srd$X)	
	# -----------
	str(erd)
	str(srd)	
# ----------------------------------------------------------------------
# Map Locations
# ----------------------------------------------------------------------	
	library(maps)
	par(mfrow=c(1,2), cex=0.5, mar=c(2,2,2,2))
	plot(erd$locs$lon, erd$locs$lat,  
		cex=0.1, 
		col="grey")
	map("world", add=T, col="blue")
	map("state", add=T, col="blue")
	par(cex=0.5, mar=c(2,2,2,2))
	plot(srd$locs$lon, srd$locs$lat,  
		cex=0.1, 
		col="grey")
	map("world", add=T, col="blue")
	map("state", add=T, col="blue")

# ----------------------------------------------------------------------
# Peek at Species Counts
# ----------------------------------------------------------------------
	str(erd$y)
	# Note that these two are in different orders! 
	erd$species
	# Hist of counts
	par(mfrow=c(3,4), cex=0.25)
	for (iii in 1:ncol(erd$y)){
		hist(erd$y[,iii], 
			breaks=100, 
			main = colnames(erd$y)[iii])
	}
# ----------------------------------------------------------------------
# Data Pre-processing 
#	Filter Effort, YEAR, etc
# 	Grid Checklists to 3km x 3km & record grid cell numbers to define 
# 		repeat checklists
# ----------------------------------------------------------------------
	# ----------------------------------------------------------------------
	# Effort Filters
	# ----------------------------------------------------------------------	
		xxx <- erd$X
		ttt.index <- 
			xxx$YEAR == 2016 & 
			xxx$EFFORT_HRS > 0 & 
			xxx$EFFORT_HRS < 5 & 
			xxx$EFFORT_DISTANCE_KM > 0 &
			xxx$EFFORT_DISTANCE_KM < 5 &
			xxx$TIME >= 5 & 
			xxx$TIME <= 12 & 
			xxx$NUMBER_OBSERVERS <= 10
		mean(ttt.index)
		sum(ttt.index)
		erd$X <- erd$X[ttt.index, ]
		erd$y <- erd$y[ttt.index, ]
		erd$locs <- erd$locs[ttt.index, ]
	# ----------------------------------------------------------------------
	# Load & merge 2016 Expertise scores
	# ----------------------------------------------------------------------	
		obs.exp <- fread(paste(data.dir, "observer_expertise.txt", sep = "") )
		obs.exp <- as.data.frame(obs.exp)
		obs.exp$OBSERVER_ID <- as.factor(obs.exp$OBSERVER_ID)
		erd$X <- merge(
				x = erd$X, 
				y = obs.exp, 
				by = "OBSERVER_ID")
		str(erd$X)		
	# ----------------------------------------------------------------------
	# Grid Locations @ 3km x 3km  
	# ----------------------------------------------------------------------	
		# Corresponds to 5km @ midpoint of study extent
		grid.cell.size <- -89.7 + 89.63826
		# Corresponds to 3km @ midpoint of study extent
		#grid.cell.size <- -89.7 + 89.73705
		# -------
		xxx <- erd$locs$lon
		yyy <- erd$locs$lat
		# Assuming a flat earth
		nx <- round(( max(xxx)-min(xxx))/grid.cell.size )
		ny <- round(( max(yyy)-min(yyy))/grid.cell.size )
		lgc <- lookup.grid.cell(
			xxx, 
			yyy, 
			#xlim Default is range of xxx   
			#ylim Default is range of yyy
			nx = nx, 
			ny = ny, 
			jitter = F )
		# Check out results
		str(lgc)
		# -----------------------------
		par(mfrow=c(1,2), cex = 0.5)
		plot(xxx, yyy, 
			pch=20, 
			col="red", 
			cex=0.5)		
		# Reconstruct the Grid that was used
		for (iii in 1:(lgc$nx+1))
			lines( rep(lgc$bb[1,1] + (iii-1)*lgc$xwidth, 2), range(lgc$bb[,"yyy"]), col="grey") 	
		for (iii in 1:(lgc$ny+1))
			lines( range(lgc$bb[,"xxx"]), rep(lgc$bb[1,2] + (iii-1)*lgc$ywidth, 2), , col="grey") 	
		# Zoom in to see labels
		ttt.index <- yyy >  43 &
				yyy < 43.3 & 
				xxx > -89.7 & 
				xxx < -89.4 
		plot(
			xxx[ttt.index], 
			yyy[ttt.index],  
			pch=20, 
			col="red", 
			cex=0.5)
		# Reconstruct the Grid that was used
		for (iii in 1:(lgc$nx+1))
			lines( rep(lgc$bb[1,1] + (iii-1)*lgc$xwidth, 2), range(lgc$bb[,"yyy"]), col="grey") 	
		for (iii in 1:(lgc$ny+1))
			lines( range(lgc$bb[,"xxx"]), rep(lgc$bb[1,2] + (iii-1)*lgc$ywidth, 2), , col="grey") 	
		text(
			xxx[ttt.index][!is.na(lgc$cell.number)], 
			yyy[ttt.index][!is.na(lgc$cell.number)],
			labels = lgc$cell.number[ttt.index][!is.na(lgc$cell.number)], 
			col = "blue") 
	# ----------------------------------------------------------------------
	# Assign grid cell numbers to ERD$X   
	# ----------------------------------------------------------------------	
	erd$X$grid.cell.number <- lgc$cell.number
	str(erd)




