
# -----------------------------------------------------------------------------
# FUNCTION - LOOKUP GRID CELL NUMBER
# INPUT: 
# 	xxx = vector of longitude or E-W coordinates
# 	yyy = vector of latitide or N-S coordinates
# 	xlim & ylim = Together these define a bounding box
# 		within which lookup occurs. I.e. all (xxx, yyy). 
# 		pairs outside of this box are ignored. 
# 	nx & ny = number of grid cells in each direction
# 	jitter = T/F randomize grid location 
#
# OUTPUT: 
# 	The output data structure is a list with the following components. 
# 
#	vector of grid cell numbers (keys) in row-major order
# 
# CORNER CASES: 
# 	partial or no points - pass pack NA if data point is out of the bounding box
#
# -----------------------------------------------------------------------------
# TEST: lookup.grid.cell
# -----------------------------------------------------------------------------
# 	nnn <- 2500
# 	xxx <- runif(nnn, 0, 10)
# 	yyy <- runif(nnn, 0, 10)
# 	par(cex = 0.5)
# 	plot(xxx, yyy, 
# 		xlim = c(-1,11), 
# 		ylim = c(-1,11), 
# 		pch=20, 
# 		col="red", 
# 		cex=0.5)
# 	lgc <- lookup.grid.cell(
# 		xxx, 
# 		yyy, 
# 		xlim = c(4,6), 
# 		ylim = c(4,6),
# 		nx = 2, 
# 		ny = 2, 
# 		jitter = T )
# 	points(
# 		xxx[!is.na(lgc$cell.number)], 
# 		yyy[!is.na(lgc$cell.number)], 
# 		col = "blue") 
# Reconstruct the Grid that was used
# 	for (iii in 1:(lgc$nx+1))
# 		lines( rep(lgc$bb[1,1] + (iii-1)*lgc$xwidth, 2), range(lgc$bb[,"yyy"]), col="grey") 	
# 	for (iii in 1:(lgc$ny+1))
# 		lines( range(lgc$bb[,"xxx"]), rep(lgc$bb[1,2] + (iii-1)*lgc$ywidth, 2), , col="grey") 	
# Timing Tests
# 	nnn <- 1000000
# 	xxx <- runif(nnn, 0, 10)
# 	yyy <- runif(nnn, 0, 10)
# 	system.time(
# 		grid.cell.number <- lookup.grid.cell(
# 			xxx, 
# 			yyy, 
# 			xlim = c(0,10), 
# 			ylim = c(0,10),
# 			nx = 100, 
# 			ny = 100, 
# 			jitter = T )  )
#
# -----------------------------------------------------------------------------
require(plyr)
lookup.grid.cell <- function(
	xxx, 
	yyy, 
	xlim = c(NA,NA), 
	ylim = c(NA,NA),
	nx = 64, 
	ny = 64, 
	jitter = F ){
	# -----------------------------------
	cell.number <- rep(NA, length(xxx))
	if (any(is.na(xlim))) xlim <- range(xxx, na.rm=T)
	if (any(is.na(ylim))) ylim <- range(yyy, na.rm=T)
	xxx.width <- abs(xlim[2]-xlim[1])/nx
	yyy.width <- abs(ylim[2]-ylim[1])/ny
	# Lower Left Corner
	x.ll <- min(xlim)
	y.ll <- min(ylim)
	# If Jittered, domain is bigger & number of grid cells increases	
	if (jitter){
		x.ll <- x.ll - runif(1)*xxx.width
		y.ll <- y.ll - runif(1)*yyy.width
		nx <- 1 + (max(xlim) - min( c(x.ll, xlim) )) %/% xxx.width
		ny <- 1 + (max(ylim) - min( c(y.ll, ylim) )) %/% yyy.width
	}
	# ID data within grid
	ingrid.index <- 
		xxx >= x.ll & xxx <= x.ll + nx*xxx.width & 
		yyy >= y.ll & yyy <= y.ll + ny*yyy.width
	if (sum(ingrid.index)>0){
		# Assign Row, Column, Grid Cell Number
		col.number <- 1 + ( xxx[ingrid.index] - x.ll ) %/% xxx.width
		row.number <- 1 + ( yyy[ingrid.index] - y.ll ) %/% yyy.width
		cell.number[ingrid.index] <- col.number + (row.number-1)*nx 
			# for (iii in 1:length(grid.cell.number)){
			# 	text(
			# 		x = xxx[iii] + xxx.width/2, 
			# 		y = yyy[iii] + yyy.width/2,
			# 		labels = paste( grid.cell.number[iii]) )
			# 	text(
			# 		x = xxx[iii]  + xxx.width/2, 
			# 		y = yyy[iii]  + 0,
			# 		labels = paste( iii ), 
			# 		col = "red" )			
			# 	}
		} 	
	return( list(
		cell.number = cell.number,
		# jittered bounding box 
		bb = matrix( 
			c(x.ll, y.ll, x.ll + nx*xxx.width, y.ll + ny*yyy.width), 2, 2, 
			byrow=T, dimnames=list(c("ll", "ur"), c("xxx", "yyy"))), 
		nx = nx, 
		ny = ny, 
		xwidth = xxx.width, 
		ywidth = yyy.width  ))
}

# -----------------------------------------------------------------------------
# This function generates a stratified sample over a regular grid of strata, 
# i.e. grid cells. It is relatively efficient computationally speaking. 
#
# INPUT:
# 	xxx = vector of longitude or E-W coordinates
# 	yyy = vector of latitide or N-S coordinates
# 	xlim & ylim = Together these define a bounding box
# 		within which lookup occurs. I.e. all (xxx, yyy). 
# 		pairs outside of this box are ignored. 
# 	nx & ny = number of grid cells in each direction
# 	jitter = T/F randomize grid location 
# 	size = maximum number of points per cell to sample
# 	replace = T/F whether to sample with replacement
#
# OUTPUT
# 	index vector of selected locations
#
# NOTE - The sample() function cannot take a sample 
# larger than the population (in a cell) when 'replace = FALSE'
# If 'replace = FALSE' and the size parameter is larger than the cell populatuion size
# this function passes back a vector of length size, 
# but it will contain only as many unique points in the cell
# and the rest of the entries will be NA's.   
#
# -----------------------------------------------------------------------------
# TEST: sample.grid.cell
# -----------------------------------------------------------------------------
# Generate Random Points over 2D field 
# nnn <- 1000
# xxx <- runif(nnn, 0, 10)
# yyy <- runif(nnn, 0, 10)
# par(cex = 0.5)
# plot(xxx, yyy, 
# 	xlim = c(-1,11), 
# 	ylim = c(-1,11), 
# 	pch=20, 
# 	col="red", 
# 	cex=0.5)
# sgc <- sample.grid.cell(
# 		xxx, 
# 		yyy, 
# 		xlim = c(-1,3), 
# 		ylim = c(3,6),
# 		nx = 5, 
# 		ny = 2, 
# 		jitter = T, 
# 		size = 1, 
# 		replace = F )	
# length(sgc$sample.index)
# sum(!is.na(sgc$sample.index))
# points(
# 	xxx[!is.na(sgc$cell.number)], 
# 	yyy[!is.na(sgc$cell.number)], 
# 	col = "blue", 
# 	pch = 20, 
# 	cex = 0.5) 	
# points(xxx[sgc$sample.index], yyy[sgc$sample.index], pch=20, cex=1, col="green")
# # Reconstruct the Grid that was used
# for (iii in 1:(sgc$nx+1))
# 	lines( rep(sgc$bb[1,1] + (iii-1)*sgc$xwidth, 2), range(sgc$bb[,"yyy"]), col="grey") 	
# for (iii in 1:(sgc$ny+1))
# 	lines( range(sgc$bb[,"xxx"]), rep(sgc$bb[1,2] + (iii-1)*sgc$ywidth, 2), , col="grey") 	
# # --------
# table(sgc$cell.number)
# length(table(sgc$cell.number))
# sort(unique(sgc$sample.index))
# sort(unique(c(1:length(xxx))[!is.na(sgc$cell.number)]))
# -----------------------------------------------------------------------------
sample.grid.cell <- function(
	xxx, 
	yyy, 
	xlim = c(NA,NA), 
	ylim = c(NA,NA),
	nx = 64, 
	ny = 64, 
	jitter = F,
	size = 1, 
	replace = F ){	
	require(plyr)
	# Stratified sample over Grid Cell Number
	sample_fun <- function(x, size, replace){
		# Cells without samples are excluded in the tapply call - if (length(x)==0) return(NA) 
		# Cells with a single sample cause problems, see help(sample)
		# So, I am going to handle this situation "by hand"
		result <- rep(NA, size)
		if (length(x)==1 & replace==F) {
			#cat("sf: length(x)==1 & replace==F",x,"\n")
			result <- rep(NA, size)
			result[1] <- x 
		}
		if (length(x)==1 & replace==T) {
			#cat("sf: length(x)==1 & replace==T",x,"\n")
			result <- rep(x, size)
		}
		if (length(x)>1 & replace == F & size > length(x) ){
			result <- rep(NA, size)
			result[1:length(x)] <- x 
		}
		if (length(x)>1 & replace == F & size <= length(x) ){
			result <- sample(x=x, size=size, replace=replace)	
		}		
		if (length(x)>1 & replace == T ){
			result <- sample(x=x, size=size, replace=replace)	
		}
		return(result) 
	}
	lgc <- lookup.grid.cell(
				xxx, yyy, xlim, ylim, nx, ny, jitter)  
	n.index <- tapply( 
		c(1:length(xxx))[!is.na(lgc$cell.number)], 
		as.factor(lgc$cell.number[!is.na(lgc$cell.number)]), 
		sample_fun, size, replace) 
	n.index <- rbind.fill.matrix(n.index)	
	return(list(
		cell.number = lgc$cell.number, 
		sample.index = n.index, 
		bb = lgc$bb, 
		nx = lgc$nx, 
		ny = lgc$ny, 
		xwidth = lgc$xwidth, 
		ywidth = lgc$ywidth  ))
} # END FUNCTION 
# -----------------------------------------------------------------------------
	



