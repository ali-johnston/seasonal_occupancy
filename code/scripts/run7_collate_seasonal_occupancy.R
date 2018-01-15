###########################################################################
###########################################################################
###### CODE TO COLLATE SEASONAL OCCUPANCY MODELS 
###########################################################################
###########################################################################

# created by: ali johnston
# on: 22.11.2017

 # load in required packages
library(unmarked)
library(dplyr)


# set parameters for whole script
data.tag <- "BCR23_2016"


# File directories in which to search for models and plot results:
data.root <- "data/"
data.folder <- paste(data.root, data.tag, "_data/", sep="")

result.root <- "results/"
result.folder <- paste(result.root, data.tag, "/", sep="")

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



###########################################################################
###########################################################################
######   READ IN PREVIOUS RESULTS
###########################################################################
###########################################################################



##################################################################### 
###  loop through species

for(sss in 1:length(specs_vec)){
 
   species <- specs_vec[sss]

   print(paste("Running for species:", species, "  species number run =", which(specs_vec==species)))
   spec.folder <- paste(result.folder, "Spec_", species, "/" , sep="")

   results_loc <- paste0(spec.folder, species, "_estimated_det_occ_by_DAY.txt")
   res <- read.table(results_loc, header=TRUE) 

   colnames(res)[2:7] <- paste0(species, "_", colnames(res)[2:7])

   if(sss==1) all_species <- res
   if(sss>1) all_species <- cbind(all_species, res[,2:7])

} # close species

write_loc <- paste0(result.folder, "all_species_DAY_occ_det.csv")
write.table(all_species, write_loc, row.names=FALSE)



