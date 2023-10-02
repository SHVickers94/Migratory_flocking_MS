#### Vickers et al. (in review) - Flocking behaviour facilitates range shifts in migratory species ####

#Stephen H. Vickers1*, Timothy D. Meehan2, Nicole L. Michel2, Aldina M.A. Franco1 and James J. Gilroy1

#1 School of Environmental Sciences, University of East Anglia, NR4 7TJ Norwich, UK.
#2 National Audubon Society, 225 Varick Street, New York, NY, 10014 USA.
#* Corresponding author: svickers@rvc.ac.uk/sjedwards94@hotmail.co.uk

# Code compiled March 2023 on R version 4.2.2 (2022-10-31) -- "Innocent and Trusting" Copyright (C) 2022 The R Foundation for Statistical Computing. Platform: x86_64-apple-darwin17.0 (64-bit). 

# 0. Setup ----
# Loading packages:
# devtools::install_github("thomasp85/patchwork")
# remotes::install_github("RS-eco/traitdata")
if(!require(pacman))install.packages('pacman');library(pacman)
pacman::p_load(tidyr,
               ggplot2,
               dplyr,
               plyr,
               formatR,
               pander,
               ggpubr,
               kableExtra,
               scico,
               cowplot,
               rnaturalearth,
               rnaturalearthdata,
               sf,
               ggspatial,
               gridExtra,
               ggnewscale,
               leaflet,
               ggeffects,
               errors,
               data.table,
               geosphere,
               rotl,
               ape,
               ggtree,
               nlme,
               phytools,
               ape,
               caper,
               rgdal,
               traitdata,
               tictoc,
               mgcv,
               truncnorm,
               readr,
               ebirdst,
               MuMIn)
options(scipen=99999999)

# 1. Production of annual centre of abundance and slope coefficients ----
## 1.1 Analysis strata ----

centroids <- st_read("~/Desktop/MIG STRAT MS sub/Data/analytical_strata[78]/bcr_by_state_clean_dissolved_centroid_points_102010.shp", quiet=T)

centroids <- centroids %>% mutate(country2 = case_when(country == 'USA' ~ 'US',
                                                       country == 'CANADA' ~ 'CA'),
                                  state_prov2 = case_when(state_prov == 'NORTH CAROLINA' ~ 'NC',
                                                          state_prov == 'SOUTH DAKOTA' ~ 'SD',
                                                          state_prov == 'OHIO' ~ 'OH',
                                                          state_prov == 'MANITOBA' ~ 'MB',
                                                          state_prov == 'MISSISSIPPI' ~ 'MS',
                                                          state_prov == 'NEVADA' ~ 'NV',
                                                          state_prov == 'ARIZONA' ~ 'AZ',
                                                          state_prov == 'WYOMING' ~ 'WY',
                                                          state_prov == 'OREGON' ~ 'OR',
                                                          state_prov == 'TEXAS' ~ 'TX',
                                                          state_prov == 'IOWA' ~ 'IA',
                                                          state_prov == 'NEBRASKA' ~ 'NE',
                                                          state_prov == 'OKLAHOMA' ~ 'OK',
                                                          state_prov == 'QUEBEC' ~ 'QC',
                                                          state_prov == 'NEW YORK' ~ 'NY',
                                                          state_prov == 'RHODE ISLAND' ~ 'RI',
                                                          state_prov == 'ALABAMA' ~ 'AL',
                                                          state_prov == 'BRITISH COLUMBIA' ~ 'BC',
                                                          state_prov == 'TENNESSE' ~ 'TN',
                                                          state_prov == 'INDIANA' ~ 'IN',
                                                          state_prov == 'YUKON' ~ 'YT',
                                                          state_prov == 'MASSACHUSETTS' ~ 'MA',
                                                          state_prov == 'CALIFORNIA' ~ 'CA',
                                                          state_prov == 'NORTHWEST TERRITORIES' ~ 'NT',
                                                          state_prov == 'GEORGIA' ~ 'GA',
                                                          state_prov == 'MINNESOTA' ~ 'MN',
                                                          state_prov == 'WISCONSIN' ~ 'WI',
                                                          state_prov == 'ILLINOIS' ~ 'IL',
                                                          state_prov == 'SASKATCHEWAN' ~ 'SK',
                                                          state_prov == 'UTAH' ~ 'UT',
                                                          state_prov == 'ALBERTA' ~ 'AB',
                                                          state_prov == 'ONTARIO' ~ 'ON',
                                                          state_prov == 'NEW HAMPSHIRE' ~ 'NH',
                                                          state_prov == 'KENTUCKY' ~ 'KY',
                                                          state_prov == 'VIRGINIA' ~ 'VA',
                                                          state_prov == 'ALASKA' ~ 'AK',
                                                          state_prov == 'ARKANSAS' ~ 'AR',
                                                          state_prov == 'NEW JERSEY' ~ 'NJ',
                                                          state_prov == 'MONTANA' ~ 'MT',
                                                          state_prov == 'MAINE' ~ 'ME',
                                                          state_prov == 'NEW MEXICO' ~ 'NM',
                                                          state_prov == 'COLORADO' ~ 'CO',
                                                          state_prov == 'LOUISIANA' ~ 'LA',
                                                          state_prov == 'MISSOURI' ~ 'MO',
                                                          state_prov == 'NOVA SCOTIA' ~ 'NS',
                                                          state_prov == 'MICHIGAN' ~ 'MI',
                                                          state_prov == 'SOUTH CAROLINA' ~ 'SC',
                                                          state_prov == 'VERMONT' ~ 'VT',
                                                          state_prov == 'WASHINGTON' ~ 'WA',
                                                          state_prov == 'CONNECTICUT' ~ 'CT',
                                                          state_prov == 'KANSAS' ~ 'KS',
                                                          state_prov == 'NEWFOUNDLAND' ~ 'NL',
                                                          state_prov == 'PRINCE EDWARD ISLAND' ~ 'PE',
                                                          state_prov == 'IDAHO' ~ 'ID',
                                                          state_prov == 'NORTH DAKOTA' ~ 'ND',
                                                          state_prov == 'MARYLAND' ~ 'MD',
                                                          state_prov == 'DELAWARE' ~ 'DE',
                                                          state_prov == 'NEW BRUNSWICK' ~ 'NB',
                                                          state_prov == 'PENNSYLVANIA' ~ 'PA',
                                                          state_prov == 'NUNAVUT' ~ 'NU',
                                                          state_prov == 'WEST VIRGINIA' ~ 'WV',
                                                          state_prov == 'FLORIDA' ~ 'FL'),
                                  Strata = paste0(country2, '-',state_prov2, '-',bcr))

# Restrict Strata by lat/long
# Reproject to get lat/longs
centroids2 <- st_transform(centroids, crs='WGS84')
coords <- data.frame(st_coordinates(centroids2$geometry))
centroids$Long <- coords$X
centroids$Lat <- coords$Y

centroids2 <- centroids %>% filter(bcr %in% c(5, 9:39), Lat < 53 & Long > -125.1 & Long < -67)

polys <- st_read("~/Desktop/MIG STRAT MS sub/Data/analytical_strata[78]/bcr_by_state_clean_dissolved_centroids_polys_102010.shp", quiet=T)

polys2 <- polys %>% filter(stratum %in% centroids2$stratum)

plot(polys2[,1], col=sf.colors(12, categorical = T), main='Strata')

# Stratum-based annual abundance metrics derived from hierarchical Bayesian models were accessed from pre-existing sources. Annual centre of abundance metrics were calculated for as many of the 409 species analysed by Beauchamp(2011) as possible. Annual COA metrics were calculated for 338 species for the BBS, 230 species for CBC, with 174 species present in both.
# We are provided with annual median abundance estimates and 95th percentiles derived from posteriors for a species within each stratum. Standard deviation was derived for these metrics by: ```SD <- (Upper_Percentile_95-Lower_Percentile_95)/3.92```
# We used classical error propagation to incorporate standard error of annual indices into COA metrics and subsequent trend analyses of COAs. Error is propagated using the first-order Taylor series method from package ```methods```.

## 1.2 BBS and CBC data handling ----

# Read in raw data files
CBC <- read_csv('Data/cbc_species_by_stratum_by_year_abundance_indices.csv') # This is CBC stratum level population estimates by species and year
BBS <- read_csv('Data/all 2019 BBS indices.csv') # This is BBS stratum level population estimates by species and year
beauchamp <- read_csv('Data/MS_data.csv') # This is data from Beauchamp (2011) combined with data taken from Birds of the World Online and oldbird.org

# BBS

# Only take strata level BBS data
BBS2 <- BBS %>% filter(Region_type == 'stratum')
names(BBS2)[c(2,23)] <- c('Strata','AOU')
BBS2 <- BBS2 %>% filter(AOU %in% unique(beauchamp$AOU))
BBS2 <- plyr::join(BBS2, beauchamp[c('AOU', 'Genus','Species')])
length(unique(BBS2$AOU)) # 338 sp.
BBS_final <- BBS2[c("Genus","Species","AOU","Strata","Year","Index","Index_q_0.025","Index_q_0.975" )]
names(BBS_final)[5:8] <- c('Year','BBS_median', 'BBS_LP95', 'BBS_UP95')

# CBC
# Create Strata names in same format as BBS
CBC <- CBC %>% mutate(country = case_when(country == 'CAN' ~ 'CA', country == 'USA' ~'US'), Strata = paste0(country, '-',state,'-',bcr))
# Split species column into parts
CBC2 <- CBC %>% separate(species, c('Genus', 'Species','Species2','Species3'), ' ', remove=F)
# Remove species where multiple species are grouped
CBC3 <- CBC2 %>% filter(is.na(Species2))
# Add AOU number for CBC species
CBC3 <- plyr::join(CBC3, beauchamp[c('Genus','Species','AOU')])
# Only retain Beauchamp species
CBC3 <- CBC3 %>% filter(!is.na(AOU))
unique(CBC3$AOU) # 284 sp.
CBC_final <- CBC3[c('Genus','Species','AOU','Strata','count_year','estimate_median',"estimate_lcl95","estimate_ucl95")]
names(CBC_final)[5:8] <- c('Year','CBC_median', 'CBC_LP95', 'CBC_UP95')
CBC_final <- CBC_final %>% filter(Year < 2020 & Year>1969)

# CBC & BBS
Both_final <- full_join(BBS_final, CBC_final)
centroid_DF <- as.data.frame(centroids2[c('Strata', 'xcoord', 'ycoord', 'Lat', 'Long')])
Both_final2 <- plyr::join(Both_final, as.data.frame(centroid_DF)[1:5])
Both_final3 <- Both_final2 %>% filter(!is.na(xcoord))
unique(Both_final3$Strata)

# Estimate SDs
Both_final3$BBS_SD <- (Both_final3$BBS_UP95-Both_final3$BBS_LP95)/3.92
Both_final3$CBC_SD <- (Both_final3$CBC_UP95-Both_final3$CBC_LP95)/3.92

#write.csv(Both_final3, 'Annual_Strata_metrics_all.csv', row.names=F)

# 1.3 Propagate error and calculate COA and slope coefficients - errors method ----
# Please note: this section does not need to be run. A large number of files will be produced but a combined file is made available and can be read in the following section.

options(errors.notation = "plus-minus", digits = 15)

#Both_final3 <- read_csv('Data/Strata_metrics_all.csv')
AOUS <- read.csv('Data/AOU_list.csv') # This is just AOU species codes we are interested in. This list was created to save computational time, and is based upon known missing traits which means the excluded species cannot be used in later analysis.

#Both_final3 <- Both_final3 %>% filter(AOU %in% AOUS$AOUS)

CBC <- Both_final3 %>% filter(!(is.na(CBC_SD))) %>% dplyr::select(1:5, 9:15, 17)
#write.csv(CBC, 'Data/Strata_metrics_CBC.csv', row.names=F)
CBC<- read_csv('Data/Strata_metrics_CBC.csv')

# Assign errors to numeric values
errors(CBC$CBC_median) <- CBC$CBC_SD
errors(CBC$xcoord) <- 100
errors(CBC$ycoord) <- 100

BBS <- Both_final3 %>% filter(!(is.na(BBS_SD))) %>% dplyr::select(1:8, 12:16)
#write.csv(BBS, 'Data/Strata_metrics_BBS.csv', row.names=F)
BBS<- read_csv('Data/Strata_metrics_BBS.csv')

errors(BBS$BBS_median) <- BBS$BBS_SD
errors(BBS$xcoord) <- 100
errors(BBS$ycoord) <- 100

length(unique(BBS$AOU))
length(unique(CBC$AOU))
BBS <- BBS %>% filter(AOU %in% unique(CBC$AOU))
CBC <- CBC %>% filter(AOU %in% unique(BBS$AOU))
length(unique(BBS$AOU))
length(unique(CBC$AOU))

options(warn = -1)
counter <- 0
templist <- list()
for(i in unique(BBS$AOU)){
  counter <- counter+1
    temp <- CBC %>% filter(AOU == i)
    temp2 <- temp %>% group_by(Year) %>% dplyr::summarise(coa_x = sum(CBC_median*xcoord)/sum(CBC_median), coa_y = sum(CBC_median*ycoord)/sum(CBC_median)) %>% dplyr::summarise(CBC_Lat_slope = sum((Year - mean(Year))*(coa_y - mean(coa_y)))/sum((Year - mean(Year))^2), CBC_Lon_slope = sum((Year - mean(Year))*(coa_x - mean(coa_x)))/sum((Year - mean(Year))^2))
    
    temp2$CBC_Lon_slope_Error <- errors(temp2$CBC_Lon_slope)
    temp2$CBC_Lat_slope_Error <- errors(temp2$CBC_Lat_slope)
    
    temp2$AOU <- i
    temp2 <- temp2[c(5, 1:4)]
    temp2 <- drop_errors(temp2)
    templist[[counter]] <- temp2

    print(paste(round(counter/length(unique(BBS$AOU)),2)))
}

CBC_slope_metrics <- rbindlist(templist)
write.csv(CBC_slope_metrics, 'Data/CBC_slope_metrics_217.csv', row.names=F)

# Create centres of abundance and linear slope coeffcients
CBC_slope_metrics <- CBC %>% group_by(AOU, Year) %>% dplyr::summarise(coa_x = sum(CBC_median*xcoord)/sum(CBC_median), coa_y = sum(CBC_median*ycoord)/sum(CBC_median)) %>% dplyr::summarise(CBC_Lat_slope = sum((Year - mean(Year))*(coa_y - mean(coa_y)))/sum((Year - mean(Year))^2), CBC_Lon_slope = sum((Year - mean(Year))*(coa_x - mean(coa_x)))/sum((Year - mean(Year))^2))

# Make the associated error its own column
CBC_slope_metrics$CBC_Lon_slope_Error <- errors(CBC_slope_metrics$CBC_Lon_slope)
CBC_slope_metrics$CBC_Lat_slope_Error <- errors(CBC_slope_metrics$CBC_Lat_slope)

#write.csv(CBC_slope_metrics, 'CBC_slope_metrics.csv', row.names = F)
CBC_slope_metrics <- read_csv('Data/CBC_slope_metrics.csv')


counter <- 0
templist <- list()
for(i in unique(BBS$AOU)){
  counter <- counter+1
  temp <- BBS %>% filter(AOU == i)
  temp2 <- temp %>% group_by(Year) %>% dplyr::summarise(coa_x = sum(BBS_median*xcoord)/sum(BBS_median), coa_y = sum(BBS_median*ycoord)/sum(BBS_median)) %>% dplyr::summarise(BBS_Lat_slope = sum((Year - mean(Year))*(coa_y - mean(coa_y)))/sum((Year - mean(Year))^2), BBS_Lon_slope = sum((Year - mean(Year))*(coa_x - mean(coa_x)))/sum((Year - mean(Year))^2))
  
  temp2$BBS_Lon_slope_Error <- errors(temp2$BBS_Lon_slope)
  temp2$BBS_Lat_slope_Error <- errors(temp2$BBS_Lat_slope)
  
  temp2$AOU <- i
  temp2 <- temp2[c(5, 1:4)]
  temp2 <- drop_errors(temp2)
  templist[[counter]] <- temp2
  
  print(paste(round(counter/length(unique(BBS$AOU)),2)))
}

BBS_slope_metrics <- rbindlist(templist)
write.csv(BBS_slope_metrics, 'Data/BBS_slope_metrics_217.csv', row.names=F)


BBS_slope_metrics <- BBS %>% group_by(AOU, Year) %>% dplyr::summarise(coa_x = sum(BBS_median*xcoord)/sum(BBS_median), coa_y = sum(BBS_median*ycoord)/sum(BBS_median)) %>% dplyr::summarise(BBS_Lat_slope = sum((Year - mean(Year))*(coa_y - mean(coa_y)))/sum((Year - mean(Year))^2), BBS_Lon_slope = sum((Year - mean(Year))*(coa_x - mean(coa_x)))/sum((Year - mean(Year))^2))

BBS_slope_metrics$BBS_Lon_slope_Error <- errors(BBS_slope_metrics$BBS_Lon_slope)
BBS_slope_metrics$BBS_Lat_slope_Error <- errors(BBS_slope_metrics$BBS_Lat_slope)

# write.csv(BBS_slope_metrics, 'BBS_slope_metrics.csv', row.names = F)
BBS_slope_metrics <- read_csv('Data/BBS_slope_metrics.csv')

# 2. Production of trait database ----

# Beauchamp (2011) dataset is available: https://royalsocietypublishing.org/doi/suppl/10.1098/rsbl.2011.0243

# Migration (long/short), Mig_distance (short/long), and Male_mass_g were taken from Beauchamp (2011). FlocksWithAdults is a three level trait (Solo, Flocks without adults, and Flocks with adults) was derived from Bird of the World Online literature search. Habitat_specialism_score and	Diet_specialism_score were derived from the Elton Traits database. GenFLS (modeled generation length) was taken from  Bird et al. 2020. Trend (loglinear BBS population trend for core range) was taken from USGS BBS data portal.

# Dataset 'beauchamp' created earlier was combined with the slope metrics datasets to produce:
trait_slopes <- read_csv('Data/trait_slopes.csv')
CBC_slope_metrics <- read_csv('Data/CBC_slope_metrics_217.csv')
BBS_slope_metrics <- read_csv('Data/BBS_slope_metrics_217.csv')

trait_slopes <- trait_slopes %>% filter(!Migration %in% c('Mostly Resident','Dispersive'))
trait_slopes <- trait_slopes %>% filter(!Combined_Migration_Timing %in% c('Both'))
Slope_metrics <- join(BBS_slope_metrics, CBC_slope_metrics)

# Tidy the data:

trait_slopes$Combined_Cohort_Timing <- as.factor(trait_slopes$Combined_Cohort_Timing)
trait_slopes$Combined_Cohort_Timing[grepl("^\\s*$", trait_slopes$Combined_Cohort_Timing)] <- NA
trait_slopes$Combined_Cohort_Timing <- droplevels(trait_slopes$Combined_Cohort_Timing)

trait_slopes$Migration <- as.factor(trait_slopes$Migration)
trait_slopes$Migration[grepl("^\\s*$", trait_slopes$Migration)] <- NA
trait_slopes$Migration <- droplevels(trait_slopes$Migration)

trait_slopes$Combined_Flock_Size <- as.factor(trait_slopes$Combined_Flock_Size)
trait_slopes$Combined_Flock_Size[grepl("^\\s*$", trait_slopes$Combined_Flock_Size)] <- NA
trait_slopes$Combined_Flock_Size <- droplevels(trait_slopes$Combined_Flock_Size)

trait_slopes$Flock_structure <- as.factor(trait_slopes$Flock_structure)
trait_slopes$Flock_structure[grepl("^\\s*$", trait_slopes$Flock_structure)] <- NA
trait_slopes$Flock_structure <- droplevels(trait_slopes$Flock_structure)

trait_slopes$Family_grouping <- as.factor(trait_slopes$Family_grouping)
trait_slopes$Family_grouping[grepl("^\\s*$", trait_slopes$Family_grouping)] <- NA
trait_slopes$Family_grouping <- droplevels(trait_slopes$Family_grouping)

trait_slopes$Altitudinal <- as.factor(trait_slopes$Altitudinal)
trait_slopes$Altitudinal[grepl("^\\s*$", trait_slopes$Altitudinal)] <- NA
trait_slopes$Altitudinal <- droplevels(trait_slopes$Altitudinal)

trait_slopes$Combined_Migration_Timing <- as.factor(trait_slopes$Combined_Migration_Timing)
trait_slopes$Combined_Migration_Timing[grepl("^\\s*$", trait_slopes$Combined_Migration_Timing)] <- NA
trait_slopes$Combined_Migration_Timing <- droplevels(trait_slopes$Combined_Migration_Timing)

trait_slopes$Combined_MFC <- as.factor(trait_slopes$Combined_MFC)
trait_slopes$Combined_MFC[grepl("^\\s*$", trait_slopes$Combined_MFC)] <- NA
trait_slopes$Combined_MFC <- droplevels(trait_slopes$Combined_MFC)

trait_slopes$Mig_distance <- as.factor(trait_slopes$Mig_distance)
trait_slopes$Mig_distance[grepl("^\\s*$", trait_slopes$Mig_distance)] <- NA
trait_slopes$Mig_distance <- droplevels(trait_slopes$Mig_distance)

# Assign errors to slope coefficients

errors(Slope_metrics$BBS_Lat_slope) <- Slope_metrics$BBS_Lat_slope_Error
errors(Slope_metrics$BBS_Lon_slope) <- Slope_metrics$BBS_Lon_slope_Error
errors(Slope_metrics$CBC_Lat_slope) <- Slope_metrics$CBC_Lat_slope_Error
errors(Slope_metrics$CBC_Lon_slope) <- Slope_metrics$CBC_Lon_slope_Error

# Use trig to calculate singular vector slope coefficient 

Slope_metrics$BBS_dist <- sqrt((Slope_metrics$BBS_Lon_slope)^2 + (Slope_metrics$BBS_Lat_slope)^2)
Slope_metrics$CBC_dist <- sqrt((Slope_metrics$CBC_Lon_slope)^2 + (Slope_metrics$CBC_Lat_slope)^2)

Slope_metrics$BBS_dist_error <- errors(Slope_metrics$BBS_dist)
Slope_metrics$CBC_dist_error <- errors(Slope_metrics$CBC_dist)
Slope_metrics$BBS_dist <- drop_errors(Slope_metrics$BBS_dist)
Slope_metrics$CBC_dist <- drop_errors(Slope_metrics$CBC_dist)
Slope_metrics <- tibble(Slope_metrics)
trait_slopes <- join(trait_slopes, Slope_metrics[c('AOU',"BBS_dist","CBC_dist","BBS_dist_error","CBC_dist_error")])

trait_slopes <- trait_slopes %>% mutate(WithAdults = case_when(Combined_Flock_Size == 'Solo'|Combined_Cohort_Timing %in% c('Adults first','Juveniles first') ~ 'No',
                                                               Combined_Flock_Size %in% c('Small','Medium','Large') & Combined_Cohort_Timing %in% c('Concurrent') ~ 'Yes'),
                                        Flocks = case_when(Combined_Flock_Size == 'Solo' ~ 'No', 
                                                           Combined_Flock_Size %in% c('Small','Medium','Large') ~ 'Yes'),
                                        FlocksWithAdults = case_when(Combined_Flock_Size == 'Solo' ~ 'Solo',
                                                                     Flocks=='Yes' & WithAdults=='Yes'~'FlocksWithAdults',
                                                                     Flocks=='Yes' & WithAdults=='No'~'FlocksNoAdults'))


# Population trends

trends <- read.csv('Data/BBS_trends.csv') 
trends <- trends[c('AOU','Region.Name','Trend','Relative.Abundance')]

trends <- trends %>% filter(Region.Name == 'Survey-wide                             ')
trends <- trends[c('AOU','Trend','Relative.Abundance')]

trait_slopes <- join(trait_slopes, trends)



BBS_trait_slopes <- trait_slopes[c('AOU','Genus','Species','BBS_dist','BBS_dist_error', 'Migration', 'Mig_distance','WithAdults', 'Flocks', 'FlocksWithAdults', 'Male_mass_g', 'Trend','Relative.Abundance', 'Combined_Flock_Size', 'Combined_Migration_Timing')]
CBC_trait_slopes <- trait_slopes[c('AOU','Genus','Species','CBC_dist','CBC_dist_error', 'Migration', 'Mig_distance','WithAdults', 'Flocks', 'FlocksWithAdults', 'Male_mass_g', 'Trend','Relative.Abundance','Combined_Flock_Size', 'Combined_Migration_Timing')]

#BBS_trait_slopes <- na.omit(BBS_trait_slopes) # 234 or 215
#CBC_trait_slopes <- na.omit(CBC_trait_slopes) # 166 or 122

BBS_trait_slopes$BBS_species <- paste(BBS_trait_slopes$Genus,BBS_trait_slopes$Species)
BBS_trait_slopes$BBS_Tree_species <- BBS_trait_slopes$BBS_species
BBS_trait_slopes$BBS_Tree_species <- sub(' ', '_',BBS_trait_slopes$BBS_Tree_species)

# 2.1 Species name fixes ----

# Due to nomenclature differences between the Hackett tree used for PGLS and that used currently by the AOU, some names need to be altered:

BBS_trait_slopes$BBS_Tree_species <- gsub('Mareca_americana','Anas_americana',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Centronyx_bairdii','Ammodramus_bairdii',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Artemisiospiza_belli','Amphispiza_belli',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Spatula_discors','Anas_discors',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Cardellina_pusilla','Wilsonia_pusilla',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Thalasseus_maximus','Sterna_maxima',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Tringa_semipalmata','Catoptrophorus_semipalmatus',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Chroicocephalus_philadelphia','Larus_philadelphia',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Leiothlypis_virginiae','Vermivora_virginiae',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Leiothlypis_peregrina','Vermivora_peregrina',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Charadrius_nivosus','Charadrius_alexandrinus',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Rhynchophanes_mccownii','Calcarius_mccownii',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Hydroprogne_caspia','Sterna_caspia',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Urile_pelagicus','Phalacrocorax_pelagicus',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Antigone_canadensis','Grus_canadensis',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Melanitta_deglandi','Melanitta_fusca',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Corthylio_calendula','Regulus_calendula',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Spinus_pinus','Carduelis_pinus',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Leiothlypis_celata','Vermivora_celata',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Leiothlypis_ruficapilla','Vermivora_ruficapilla',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Leiothlypis_luciae','Vermivora_luciae',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Spatula_cyanoptera','Anas_cyanoptera',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Mareca_strepera','Anas_strepera',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Ammospiza_leconteii','Ammodramus_leconteii',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Circus_hudsonius','Circus_cyaneus',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Sternula_antillarum','Sterna_antillarum',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Parkesia_noveboracensis','Seiurus_noveboracensis',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Parkesia_motacilla','Seiurus_motacilla',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Leucophaeus_atricilla','Larus_atricilla',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Ardea_alba','Casmerodius_albus',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Antrostomus_carolinensis','Caprimulgus_carolinensis',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Cardellina_canadensis','Wilsonia_canadensis',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Nannopterum_auritus','Phalacrocorax_auritus',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Leucophaeus_pipixcan','Larus_pipixcan',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Gallinula_galeata','Gallinula_chloropus',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Antrostomus_vociferus','Caprimulgus_vociferus',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Geothlypis_philadelphia','Oporornis_philadelphia',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Geothlypis_tolmiei','Oporornis_tolmiei',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Geothlypis_formosa','Oporornis_formosus',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Setophaga_palmarum','Dendroica_palmarum',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Setophaga_pinus','Dendroica_pinus',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Setophaga_discolor','Dendroica_discolor',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Setophaga_townsendi','Dendroica_townsendi',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Setophaga_dominica','Dendroica_dominica',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Setophaga_magnolia','Dendroica_magnolia',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Setophaga_caerulescens','Dendroica_caerulescens',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Setophaga_striata','Dendroica_striata',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Setophaga_tigrina','Dendroica_tigrina',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Setophaga_pensylvanica','Dendroica_pensylvanica',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Setophaga_occidentalis','Dendroica_occidentalis',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Setophaga_citrina','Wilsonia_citrina',BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$BBS_Tree_species <- gsub('Setophaga_americana','Parula_americana',BBS_trait_slopes$BBS_Tree_species)


CBC_trait_slopes$CBC_species <- paste(CBC_trait_slopes$Genus,CBC_trait_slopes$Species)
CBC_trait_slopes$CBC_Tree_species <- CBC_trait_slopes$CBC_species
CBC_trait_slopes$CBC_Tree_species <- sub(' ', '_',CBC_trait_slopes$CBC_Tree_species)

# Species name fixes ----
CBC_trait_slopes$CBC_Tree_species <- gsub('Mareca_americana','Anas_americana',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Centronyx_bairdii','Ammodramus_bairdii',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Spatula_discors','Anas_discors',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Cardellina_pusilla','Wilsonia_pusilla',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Thalasseus_maximus','Sterna_maxima',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Tringa_semipalmata','Catoptrophorus_semipalmatus',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Leiothlypis_peregrina','Vermivora_peregrina',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Charadrius_nivosus','Charadrius_alexandrinus',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Rhynchophanes_mccownii','Calcarius_mccownii',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Hydroprogne_caspia','Sterna_caspia',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Urile_pelagicus','Phalacrocorax_pelagicus',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Antigone_canadensis','Grus_canadensis',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Melanitta_deglandi','Melanitta_fusca',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Corthylio_calendula','Regulus_calendula',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Leiothlypis_celata','Vermivora_celata',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Leiothlypis_ruficapilla','Vermivora_ruficapilla',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Spatula_cyanoptera','Anas_cyanoptera',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Mareca_strepera','Anas_strepera',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Ammospiza_leconteii','Ammodramus_leconteii',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Parkesia_noveboracensis','Seiurus_noveboracensis',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Leucophaeus_atricilla','Larus_atricilla',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Ardea_alba','Casmerodius_albus',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Antrostomus_carolinensis','Caprimulgus_carolinensis',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Leucophaeus_pipixcan','Larus_pipixcan',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Antrostomus_vociferus','Caprimulgus_vociferus',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Setophaga_palmarum','Dendroica_palmarum',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Setophaga_pinus','Dendroica_pinus',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Setophaga_discolor','Dendroica_discolor',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Setophaga_townsendi','Dendroica_townsendi',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Setophaga_occidentalis','Dendroica_occidentalis',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Setophaga_americana','Parula_americana',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Acanthis_flammea','Carduelis_flammea',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Acanthis_hornemanni','Carduelis_hornemanni',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Hydrocoloeus_minutus','Larus_minutus',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Anser_rossii','Chen_rossii',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Bubo_scandiacus','Bubo_scandiaca',CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$CBC_Tree_species <- gsub('Calidris_virgata','Aphriza_virgata',CBC_trait_slopes$CBC_Tree_species)

# Take data from elton traits database
data(elton_birds)
elton_birds <- elton_birds
elton_birds <- elton_birds %>% rowwise() %>% mutate(Habitat_specialism_score = max(across(.cols = 25:31), na.rm=T), Diet_specialism_score = max(across(.cols = 11:20), na.rm=T))
elton_birds <- elton_birds %>% filter(English != 'Blue-winged Warbler')

elton_birds2 <- elton_birds[c('scientificNameStd','Habitat_specialism_score','Diet_specialism_score','Diet.5Cat')]
elton_birds2$scientificNameStd <- sub(' ', '_',elton_birds2$scientificNameStd)

elton_specs <- as.character(elton_birds2$scientificNameStd)
elton_specs <- sub(' ', '_',elton_specs)
elton_specs <- data.frame(elton_specs)

# Elton name changes
elton_birds2$scientificNameStd <- gsub('Hydroprogne_caspia','Sterna_caspia',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Acanthis_flammea','Carduelis_flammea',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Leucophaeus_pipixcan','Larus_pipixcan',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Ardea_alba','Casmerodius_albus',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Leiothlypis_ruficapilla','Vermivora_ruficapilla',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Pyrola_americana','Parula_americana',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Parkesia_noveboracensis','Seiurus_noveboracensis',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Leiothlypis_celata','Vermivora_celata',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Setophaga_palmarum','Dendroica_palmarum',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Setophaga_pinus','Dendroica_pinus',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Setophaga_discolor','Dendroica_discolor',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Bubo_scandiacus','Bubo_scandiaca',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Leiothlypis_peregrina','Vermivora_peregrina',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Rhynchophanes_mccownii','Calcarius_mccownii',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Setophaga_townsendi','Dendroica_townsendi',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Cardellina_pusilla','Wilsonia_pusilla',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Setophaga_occidentalis','Dendroica_occidentalis',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Acanthis_hornemanni','Carduelis_hornemanni',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Artemisiospiza_belli','Amphispiza_belli',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Setophaga_caerulescens','Dendroica_caerulescens',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Setophaga_striata','Dendroica_striata',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Chroicocephalus_philadelphia','Larus_philadelphia',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Cardellina_canadensis','Wilsonia_canadensis',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Setophaga_tigrina','Dendroica_tigrina',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Setophaga_pensylvanica','Dendroica_pensylvanica',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Setophaga_citrina','Wilsonia_citrina',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Geothlypis_formosa','Oporornis_formosus',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Sternula_antillarum','Sterna_antillarum',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Parkesia_motacilla','Seiurus_motacilla',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Leiothlypis_luciae','Vermivora_luciae',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Geothlypis_tolmiei','Oporornis_tolmiei',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Setophaga_magnolia','Dendroica_magnolia',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Geothlypis_philadelphia','Oporornis_philadelphia',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Spinus_pinus','Carduelis_pinus',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Leiothlypis_virginiae','Vermivora_virginiae',elton_birds2$scientificNameStd)
elton_birds2$scientificNameStd <- gsub('Setophaga_dominica','Dendroica_dominica',elton_birds2$scientificNameStd)

# Generation length

gen_lengths <- read.csv('Data/bird_gen_lengths.csv') 

gen_specs <- as.character(gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- sub(' ', '_',gen_lengths$Scientific.name)
gen_specs <- sub(' ', '_',gen_specs)
gen_specs <- data.frame(gen_specs)

gen_lengths$Scientific.name <- gsub('Mareca_americana','Anas_americana', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Spatula_discors','Anas_discors', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Mareca_strepera','Anas_strepera', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Spatula_cyanoptera','Anas_cyanoptera', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Passerculus_bairdii','Ammodramus_bairdii', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Ammospiza_leconteii','Ammodramus_leconteii', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Haematopus_palliatus','Haematopus_bachmani', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Himantopus_himantopus','Himantopus_mexicanus', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Hydroprogne_caspia','Sterna_caspia', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Thalasseus_maximus','Sterna_maxima', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Antrostomus_carolinensis','Caprimulgus_carolinensis', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Antrostomus_vociferus','Caprimulgus_vociferus', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Acanthis_flammea','Carduelis_flammea', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Ardea_alba','Casmerodius_albus', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Butorides_striata','Butorides_virescens', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Setophaga_occidentalis','Dendroica_occidentalis', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Setophaga_palmarum','Dendroica_palmarum', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Setophaga_pinus','Dendroica_pinus', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Setophaga_discolor','Dendroica_discolor', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Setophaga_townsendi','Dendroica_townsendi', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Hydrocoloeus_minutus','Larus_minutus', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Leiothlypis_ruficapilla','Vermivora_ruficapilla', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Leiothlypis_celata','Vermivora_celata', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Leiothlypis_peregrina','Vermivora_peregrina', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Setophaga_americana','Parula_americana', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Parkesia_noveboracensis','Seiurus_noveboracensis', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Urile_pelagicus','Phalacrocorax_pelagicus', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Anser_rossii','Chen_rossii', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Antigone_canadensis','Grus_canadensis', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Bubo_scandiacus','Bubo_scandiaca', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Calidris_virgata','Aphriza_virgata', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Rhynchophanes_mccownii','Calcarius_mccownii', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Charadrius_semipalmatus','Catoptrophorus_semipalmatus', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Cardellina_pusilla','Wilsonia_pusilla', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Artemisiospiza_belli','Amphispiza_belli', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Setophaga_caerulescens','Dendroica_caerulescens', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Setophaga_striata','Dendroica_striata', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Setophaga_tigrina','Dendroica_tigrina', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Setophaga_pensylvanica','Dendroica_pensylvanica', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Setophaga_dominica','Dendroica_dominica', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Setophaga_magnolia','Dendroica_magnolia', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Icterus_bullockiorum','Icterus_bullockii', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Selasphorus_calliope','Stellula_calliope', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Cardellina_canadensis','Wilsonia_canadensis', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Nannopterum_auritus','Phalacrocorax_auritus', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Setophaga_citrina','Wilsonia_citrina', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Geothlypis_formosa','Oporornis_formosus', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Sternula_antillarum','Sterna_antillarum', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Parkesia_motacilla','Seiurus_motacilla', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Leiothlypis_luciae','Vermivora_luciae', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Geothlypis_tolmiei','Oporornis_tolmiei', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Geothlypis_philadelphia','Oporornis_philadelphia', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Spatula_clypeata','Anas_clypeata', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Spinus_pinus','Carduelis_pinus', gen_lengths$Scientific.name)
gen_lengths$Scientific.name <- gsub('Leiothlypis_virginiae','Vermivora_virginiae', gen_lengths$Scientific.name)

names(elton_birds2)[1] <- names(BBS_trait_slopes)[17]
names(gen_lengths)[7] <- names(BBS_trait_slopes)[17]
gen_lengths2 <- gen_lengths[c('BBS_Tree_species','GenFLS','Order')]

traits_birds <- join(elton_birds2, gen_lengths2)

BBS_trait_slopes <- join(BBS_trait_slopes, traits_birds)

names(traits_birds)[1] <- 'CBC_Tree_species'
CBC_trait_slopes <- join(CBC_trait_slopes, traits_birds)


BBS_trait_slopes_flocks <- BBS_trait_slopes[c(-8,-10)]
CBC_trait_slopes_flocks <- CBC_trait_slopes[c(-8,-10)]

#BBS_trait_slopes <- na.omit(BBS_trait_slopes)
#CBC_trait_slopes <- na.omit(CBC_trait_slopes)
#BBS_trait_slopes_flocks <- na.omit(BBS_trait_slopes_flocks)
#CBC_trait_slopes_flocks <- na.omit(CBC_trait_slopes_flocks)

BBS_trait_slopes <- BBS_trait_slopes[!is.na(BBS_trait_slopes$BBS_dist),]
BBS_trait_slopes_flocks <- BBS_trait_slopes_flocks[!is.na(BBS_trait_slopes_flocks$BBS_dist),]
CBC_trait_slopes <- CBC_trait_slopes[!is.na(CBC_trait_slopes$CBC_dist),]
CBC_trait_slopes_flocks <- CBC_trait_slopes_flocks[!is.na(CBC_trait_slopes_flocks$CBC_dist),]

# We only assess species available in both seasons:

BBS_trait_slopes_flocks <- BBS_trait_slopes_flocks %>% filter(AOU %in% CBC_trait_slopes_flocks$AOU)

## Banding for migratory cohort timing overlap 
# The results of this section are already within the 'Beauchamp' dataset read in during an earlier section. 

band_check <- trait_slopes[!is.na(trait_slopes$Banding_cohort_timing),]$Common_name
band_check2 <- gsub(' ','_',band_check)
trait_slopes$Overlap_timing <- NA

for(p in 1:length(band_check)){
  
  data <- read_csv(paste0("~/Desktop/MIG STRAT MS sub/Data/Banding_data/",band_check2[p], '_banding_data.csv')) # replace with correct working directory in your system
  
  hy <- data %>% filter(Age == 'Hatch-year', EVENT_MONTH > 5)
  ahy <- data %>% filter(Age == 'After hatch-year')
  data <- rbind(hy, ahy)
  
  start <- beauchamp %>% filter(AOU ==  trait_slopes[trait_slopes$Common_name == band_check[p],]$AOU) %>% dplyr::select(Start)
  end <- beauchamp %>% filter(AOU ==  trait_slopes[trait_slopes$Common_name == band_check[p],]$AOU) %>% dplyr::select(End)
  
  data <- data %>% filter(JulianDay >= start$Start, JulianDay <= end$End)
  
  mod1 <- gam(LAT_DD ~ s(JulianDay) + s(JulianDay, by=factor(Age)) + Age, data=data)
  
  newdata=expand.grid(JulianDay=seq(min(data$JulianDay, na.rm=T), max(data$JulianDay, na.rm=T),1), Age=unique(data$Age))
  newdata$prediction <- predict(mod1, newdata = newdata)
  
  # Calculate the amount of overlap between two distributions in R
  x <- seq(min(data$JulianDay, na.rm=T), max(data$JulianDay, na.rm=T),1)
  y1 <- (newdata[newdata$Age=='Hatch-year',]$prediction-min(newdata$prediction))/(max(newdata$prediction)-min(newdata$prediction))
  y2 <- (newdata[newdata$Age=='After hatch-year',]$prediction-min(newdata$prediction))/(max(newdata$prediction)-min(newdata$prediction))
  
  dat <- data.frame(JulianDay=x,Hatch_year=y1, After_hatch_year=y2)
  dat <- melt(dat, id='JulianDay')
  dat$variable <- relevel(dat$variable, ref='After_hatch_year')
  
  # vectors to store the area under the curves and the total
  # under = min of the two vlues
  # total = max of the two values
  auc <- rep(NA, length(x))
  total <- rep(NA, length(x))
  for(i in 1:length(auc)){
    auc[i] <- min(c(y1[i], y2[i]))
    total[i] <- max(c(y1[i], y2[i]))
  }
  # add zeroes to start and end of both vectors
  # so that the y values and x values meet at the
  # bottom
  auc <- c(0,auc,0)
  total <- c(0, total, 0)
  # repeat the min and max of x for this reason
  x_area <- c(min(x), x, max(x))
  
  # a function to calculate area from xy coords
  # via a contour integral.
  area<-function(X){
    X<-rbind(X,X[1,])
    x<-X[,1]
    y<-X[,2] 
    lx<-length(x)
    abs(sum((x[2:lx]-x[1:lx-1])*(y[2:lx]+y[1:lx-1]))/2)
  }
  
  # area under the two curves
  auc_area <- area(cbind(x_area, auc))
  # equals 0.5929751
  
  # total area
  total_area <- area(cbind(x_area, total))
  # equals 1.201519
  
  trait_slopes[trait_slopes$Common_name == band_check[p],]$Overlap_timing <- auc_area/total_area
  
  print(round(p/length(band_check),3))
}

trait_slopes <- trait_slopes %>% mutate(Banding_cohort_timing2 = case_when(Overlap_timing >= 0.85 ~ 'Concurrent',Overlap_timing < 0.85 ~ 'Non-concurrent'))

trait_slopes$Overlap_timing <- round(trait_slopes$Overlap_timing,2)

trait_slopes <- trait_slopes %>% mutate(Combined_Cohort_Timing2 = case_when(
  !is.na(BOTW_cohort_timing) ~ BOTW_cohort_timing,
  is.na(BOTW_cohort_timing) ~ Banding_cohort_timing2))

compare_banding <- trait_slopes[!is.na(trait_slopes$Banding_cohort_timing) & !is.na(trait_slopes$BOTW_cohort_timing), ]

trait_slopes <- trait_slopes %>% mutate(WithAdults = case_when(Combined_Flock_Size == 'Solo'|Combined_Cohort_Timing2 %in% c('Adults first','Juveniles first', 'Non-concurrent') ~ 'No',
                                                                         Combined_Flock_Size %in% c('Small','Medium','Large') & Combined_Cohort_Timing2 %in% c('Concurrent') ~ 'Yes'),
                                                  Flocks = case_when(Combined_Flock_Size == 'Solo' ~ 'No', 
                                                                     Combined_Flock_Size %in% c('Small','Medium','Large') ~ 'Yes'),
                                                  FlocksWithAdults = case_when(Combined_Flock_Size == 'Solo' ~ 'Solo',
                                                                               Flocks=='Yes' & WithAdults=='Yes'~'FlocksWithAdults',
                                                                               Flocks=='Yes' & WithAdults=='No'~'FlocksNoAdults'))

write.csv(trait_slopes, 'Data/trait_slopes_217.csv', row.names=F)

trait_slopes <- trait_slopes %>% filter(AOU %in% BBS_trait_slopes_flocks$AOU)

# 3. PGLS trait analysis ----

BBS_trait_slopes_flocks <- read.csv('Data/BBS_trait_slopes_flocks.csv')
BBS_trait_slopes <- read.csv('Data/BBS_trait_slopes.csv')
CBC_trait_slopes_flocks <- read.csv('Data/CBC_trait_slopes_flocks.csv')
CBC_trait_slopes <- read.csv('Data/CBC_trait_slopes.csv')

# Range size and mig distance ----

#set_ebirdst_access_key("mq9v44oricfm", overwrite = T)
ebirdst_runs <- ebirdst_runs

specs <- paste(BBS_trait_slopes_flocks$Genus, BBS_trait_slopes_flocks$Species)

ebirdst_runs2 <- ebirdst_runs %>% filter(scientific_name %in% specs)

for (i in 1:nrow(ebirdst_runs2)){
  path <- ebirdst_download(species = ebirdst_runs2$scientific_name[i], path = '~/Desktop/MIG STRAT MS sub/Data/eBird', pattern='range')}
# seasonal, low res, smoothed ranges


polys3 <- st_union(polys2)

ranges <- load_ranges(path = paste0('~/Desktop/MIG STRAT MS sub/Data/eBird/2021/',ebirdst_runs2$species_code[1]), resolution = "lr")
range_breeding <- filter(ranges, season == "breeding")
polys3 <- st_transform(polys3, crs=st_crs(range_breeding))
polys4 <- st_make_valid(polys3)

ebirdst_runs2$Breeding_rangesize <- NA
ebirdst_runs2$Nonbreeding_rangesize <- NA
ebirdst_runs2$Breeding_overlap_size <- NA
ebirdst_runs2$Nonbreeding_overlap_size <- NA
ebirdst_runs2$Mig_dist_numeric <- NA

Poly_Coord_df <- data.frame(X=c(-130, -6), Y=c(-90,90))

pol = st_polygon(
  list(
    cbind(
      Poly_Coord_df$X[c(1,2,2,1,1)], 
      Poly_Coord_df$Y[c(1,1,2,2,1)])
  )
)
polc = st_sfc(pol, crs=st_crs(polys4))

BL_ranges <- st_read('Data/BirdLife_ranges.shp')
BL_ranges2 <- st_read('Data/BirdLife_ranges_2.shp')
BL_ranges <- rbind(BL_ranges, BL_ranges2)
rm(BL_ranges2)
BL_ranges <- filter(BL_ranges, SCINAME %in% c(miss_specs))

for (i in 1:nrow(ebirdst_runs2)){
  ranges <- load_ranges(path = paste0('~/Desktop/MIG STRAT MS sub/Data/eBird/2021/',ebirdst_runs2$species_code[i]), resolution = "lr")
  
  if(nrow(ranges) ==1){
    range_breeding <- range_nonbreeding <- ranges
  }else{
    range_breeding <- filter(ranges, season == "breeding")
    range_nonbreeding <- filter(ranges, season == "nonbreeding")}
  
  range_breeding <- st_intersection(range_breeding, polc)
  range_nonbreeding <- st_intersection(range_nonbreeding, polc)
  
#  ggplot(polys3) +
#    geom_sf(alpha=0.25, fill='black') +
#    geom_sf(range_breeding,  mapping=aes(), alpha=0.5, fill='red') +
#    geom_sf(range_nonbreeding, mapping=aes(), alpha=0.25, fill='blue') +
#    geom_sf(breeding_overlap, mapping=aes(), alpha=0.25, fill='darkred') #+
#    geom_sf(nonbreeding_overlap, mapping=aes(), alpha=0.25, fill#='darkblue') +
#    geom_sf(polc, mapping=aes(), alpha=0.25, fill='darkblue') 
  
  if(ebirdst_runs2$breeding_range_modeled[i] == TRUE){
    breeding_overlap <- st_intersection(polys4, range_breeding)
    ebirdst_runs2$Breeding_rangesize[i] <- st_area(range_breeding)
    ebirdst_runs2$Breeding_overlap_size[i] <- st_area(breeding_overlap)}else{
      sp_ranges <- BL_ranges %>% filter(SCINAME == ebirdst_runs2$scientific_name[i])
      range_breeding <- filter(sp_ranges, SEASONA %in% c(1,2))
      # SEASONA 1 = native resident, 2 = native breeding, 3 = native non-breeding, 4 = passage
      range_breeding <- st_make_valid(range_breeding)
      range_breeding <- st_union(range_breeding)
      range_breeding <- st_make_valid(range_breeding)
      range_breeding <- st_intersection(range_breeding, polc)
      range_breeding <- st_make_valid(range_breeding)
      range_breeding <- range_breeding[!st_is_empty(range_breeding),drop=F]
      breeding_overlap <- st_intersection(polys4, range_breeding)
      breeding_overlap <- st_make_valid(breeding_overlap)
      breeding_overlap <- breeding_overlap[!st_is_empty(breeding_overlap),drop=F]
      
      ebirdst_runs2$Breeding_rangesize[i] <- st_area(range_breeding)
      ebirdst_runs2$Breeding_overlap_size[i] <- st_area(breeding_overlap)
    }
  
  if(ebirdst_runs2$nonbreeding_range_modeled[i] == TRUE){
    nonbreeding_overlap <- st_intersection(polys4, range_nonbreeding)
    ebirdst_runs2$Nonbreeding_rangesize[i] <- st_area(range_nonbreeding)
    ebirdst_runs2$Nonbreeding_overlap_size[i] <- st_area(nonbreeding_overlap)}else{
      sp_ranges <- BL_ranges %>% filter(SCINAME == ebirdst_runs2$scientific_name[i])
      range_nonbreeding <- filter(sp_ranges, SEASONA %in% c(1,3))
      range_nonbreeding <- st_make_valid(range_nonbreeding)
      range_nonbreeding <- st_union(range_nonbreeding)
      range_nonbreeding <- st_make_valid(range_nonbreeding)
      range_nonbreeding <- st_intersection(range_nonbreeding, polc)
      range_nonbreeding <- st_make_valid(range_nonbreeding)
      range_nonbreeding <- range_nonbreeding[!st_is_empty(range_nonbreeding),drop=F]
      nonbreeding_overlap <- st_intersection(polys4, range_nonbreeding)
      nonbreeding_overlap <- st_make_valid(nonbreeding_overlap)
      ebirdst_runs2$Nonbreeding_rangesize[i] <- st_area(range_nonbreeding)
      ebirdst_runs2$Nonbreeding_overlap_size[i] <- st_area(nonbreeding_overlap)
    }
  

    ebirdst_runs2$Mig_dist_numeric[i] <- st_distance(st_centroid(range_breeding), st_centroid(range_nonbreeding))/1000
  
  print(i/nrow(ebirdst_runs2))
  
}

ebirdst_runs2$Breeding_overlap_prop <- ebirdst_runs2$Breeding_overlap_size/ebirdst_runs2$Breeding_rangesize
ebirdst_runs2$Nonbreeding_overlap_prop <- ebirdst_runs2$Nonbreeding_overlap_size/ebirdst_runs2$Nonbreeding_rangesize

trait_slopes <- read_csv('Data/trait_slopes.csv')
trait_slopes$scientific_name <- paste(trait_slopes$Genus, trait_slopes$Species)
ebirdst_runs2 <- left_join(ebirdst_runs2, trait_slopes[c('scientific_name','Combined_Migration_Timing')], )

names(ebirdst_runs2)[2] <- names(BBS_trait_slopes_flocks)[12]
BBS_trait_slopes_flocks <- left_join(BBS_trait_slopes_flocks, ebirdst_runs2[c(2, 24:31)], by=names(ebirdst_runs2)[2])
BBS_trait_slopes <- left_join(BBS_trait_slopes, ebirdst_runs2[c(2, 24:31)], by=names(ebirdst_runs2)[2])

names(ebirdst_runs2)[2] <- names(CBC_trait_slopes_flocks)[12]
CBC_trait_slopes_flocks <- left_join(CBC_trait_slopes_flocks, ebirdst_runs2[c(2, 24:31)], by=names(ebirdst_runs2)[2])
CBC_trait_slopes <- left_join(CBC_trait_slopes, ebirdst_runs2[c(2, 24:31)], by=names(ebirdst_runs2)[2])


write.csv(BBS_trait_slopes_flocks,'Data/BBS_trait_slopes_flocks.csv', row.names = F)
write.csv(BBS_trait_slopes,'Data/BBS_trait_slopes.csv', row.names = F)
write.csv(CBC_trait_slopes_flocks,'Data/CBC_trait_slopes_flocks.csv', row.names = F)
write.csv(CBC_trait_slopes,'Data/CBC_trait_slopes.csv', row.names = F)

## 3.1 BBS (just flocks) ----


BBS_trait_slopes_flocks <- read.csv('Data/BBS_trait_slopes_flocks.csv')
BBS_trait_slopes <- read.csv('Data/BBS_trait_slopes.csv')
CBC_trait_slopes_flocks <- read.csv('Data/CBC_trait_slopes_flocks.csv')
CBC_trait_slopes <- read.csv('Data/CBC_trait_slopes.csv')

#### Missing migration timing
# Used full 401 sp. dataset to fill gaps

### Bewick's Wren Thryomanes bewickii 
## All wrens nocturnal 
# Nocturnal

### Long-billed Curlew Numenius americanus
## All Numenius nocturnal
# Nocturnal

### Marbled Godwit Limosa fedoa
## Limosa mostly nocturnal
# Nocturnal

### Townsend's Solitaire Myadestes townsendi
## No other Myadestes. All other Turdidae nocturnal
# Nocturnal

### Thick-billed Longspur Rhynchophanes mccownii
## Longspurs mostly nocturnal
# Nocturnal

### Short-eared Owl Asio flammeus
## Owls mostly Nocturnal
# Nocturnal

### Little Blue Heron Egretta caerulea
## Other Egretta nocturnal 
# Nocturnal

CBC_trait_slopes[is.na(CBC_trait_slopes$Combined_Migration_Timing),]$Combined_Migration_Timing <- 'night'
CBC_trait_slopes_flocks[is.na(CBC_trait_slopes_flocks$Combined_Migration_Timing),]$Combined_Migration_Timing <- 'night'
BBS_trait_slopes[is.na(BBS_trait_slopes$Combined_Migration_Timing),]$Combined_Migration_Timing <- 'night'
BBS_trait_slopes_flocks[is.na(BBS_trait_slopes_flocks$Combined_Migration_Timing),]$Combined_Migration_Timing <- 'night'

CBC_trait_slopes <- na.omit(CBC_trait_slopes)
CBC_trait_slopes_flocks <- na.omit(CBC_trait_slopes_flocks)
BBS_trait_slopes <- na.omit(BBS_trait_slopes)
BBS_trait_slopes_flocks <- na.omit(BBS_trait_slopes_flocks)

### 3.1.1 Full model ----

species_list <- unique(CBC_trait_slopes_flocks$CBC_Tree_species)
species_list <- sub(' ', '_',species_list)

BBS_trait_slopes_flocks <- BBS_trait_slopes_flocks %>% filter(AOU %in% CBC_trait_slopes_flocks$AOU)

# Load bird supertree based on Hackett's backbone 
bird_tree_hackett <- read.tree("Data/Hackett.tre")

bird_tree_hackett <- collapse.singles(bird_tree_hackett)

# Get only bird species from the supertree that are also included in collected data
# Extract list of all species from the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)
bird_tree_species <- data.frame(Name = bird_tree_species)

bird_tree_species2 <- bird_tree_species %>% tidyr::separate(Name, c('Genus', 'Species'))
bird_tree_species$Genus <- bird_tree_species2$Genus
bird_tree_species$Species <- bird_tree_species2$Species

# Check the overlap of species names between collected data file and the supertree
# All species should be present. If not, species names may not match with names in the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)

#intersect(bird_tree_species, species_list)
#setdiff(BBS_trait_slopes$BBS_Tree_species,bird_tree_species)

# Prune supertree to the list of taxa included in the data
pruned_birds_stree <- drop.tip(bird_tree_hackett, bird_tree_hackett$tip.label[-match(BBS_trait_slopes_flocks$BBS_Tree_species, bird_tree_hackett$tip.label)])

#is.binary(pruned_birds_stree) # TRUE
#is.ultrametric(pruned_birds_stree) # TRUE

treedat <- data.frame(BBS_Tree_species=pruned_birds_stree$tip.label)
treedat <- join(treedat, BBS_trait_slopes_flocks[c('BBS_Tree_species','Flocks')])
names(treedat)[1] <- 'label'
pruned_birds_stree2 <- full_join(pruned_birds_stree,treedat, by='label')

pruned_birds_stree4 <- pruned_birds_stree2
pruned_birds_stree4@phylo$tip.label <- gsub('_',' ',pruned_birds_stree4@phylo$tip.label)

SE <- BBS_trait_slopes_flocks$BBS_dist_error
SE<-setNames(SE,BBS_trait_slopes_flocks$BBS_Tree_species)
BBS_trait_slopes_flocks$SE <- SE

BBS_trait_slopes_flocks <- BBS_trait_slopes_flocks[match(pruned_birds_stree$tip.label, BBS_trait_slopes_flocks$BBS_Tree_species),]

rownames(BBS_trait_slopes_flocks) <- BBS_trait_slopes_flocks$BBS_Tree_species

BBS_trait_slopes_flocks2 <- BBS_trait_slopes_flocks %>% mutate(
  Habitat_specialism_score = scale(Habitat_specialism_score),
  Diet_specialism_score = scale(Diet_specialism_score),
  GenFLS = scale(GenFLS),
  Trend = scale(abs(Trend)),
  Relative.Abundance = scale(Relative.Abundance),
  Breeding_rangesize = scale(Breeding_rangesize),
  Breeding_overlap_size = scale(Breeding_overlap_size),
  Mig_dist_numeric = scale(Mig_dist_numeric)
)

BBSfit<-pgls.SEy((BBS_dist) ~ Migration
                 + Mig_dist_numeric
                 + Flocks
                 + Habitat_specialism_score
                 + Diet_specialism_score
                 + GenFLS
                 + Trend
                 + Breeding_rangesize
                 + Breeding_overlap_size
                 + Combined_Migration_Timing
                 ,data=BBS_trait_slopes_flocks2,se=(SE),tree=pruned_birds_stree,method="ML")

car::vif(BBSfit) # Check for vif > 5. None so here.

BBScoefs <- data.frame(confint(BBSfit))
BBScoefs$Coef <- coef(BBSfit)
BBScoefs$Variable <- rownames(BBScoefs)
BBScoefs <- BBScoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  Variable == '(Intercept)' ~ 'Intercept',
  Variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  Variable == 'Mig_dist_numeric' ~ 'Migratory distance',
  Variable == 'FlocksNo' ~ 'Solo',
  Variable == 'FlocksYes' ~ 'Flocking\n(ref: solo)',
  Variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  Variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  Variable == 'GenFLS' ~ 'Generation length',
  Variable == 'Trend' ~ 'Absolute population trend',
  Variable == 'Breeding_rangesize' ~ 'Range size km2',
  Variable == 'Combined_Migration_Timingnight' ~ 'Migration timing\n(ref: day)'
))

BBSresults1 <- ggplot(BBScoefs, aes(x=reorder(Variable,Coef), y=Coef, colour=Significance)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(alpha=0.75) +
  geom_point() +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0) +
  labs(y='Coefficient', x='', title='Breeding model A', subtitle='122 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="none")

BBSresults1
BBScoefs_flocks <- BBScoefs

### 3.1.2 Backwards stepwise deletion ----

BBSfitBSD<-pgls.SEy((BBS_dist) ~ Habitat_specialism_score
                 + Trend
                 + Breeding_rangesize
                 ,data=BBS_trait_slopes_flocks2,se=(SE),tree=pruned_birds_stree,method="ML")

BBScoefs <- data.frame(confint(BBSfitBSD))
BBScoefs$Coef <- coef(BBSfitBSD)
BBScoefs$Variable <- rownames(BBScoefs)
BBScoefs <- BBScoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  Variable == '(Intercept)' ~ 'Intercept',
  Variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  Variable == 'Mig_distanceshort' ~ 'Short-distance migrant\n(ref: long-distance)',
  Variable == 'FlocksNo' ~ 'Solo',
  Variable == 'FlocksYes' ~ 'Flocks\n(ref: solo)',
  Variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  Variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  Variable == 'GenFLS' ~ 'Generation length',
  Variable == 'Trend' ~ 'Absolute population trend'
))

BBSresults1 <- ggplot(BBScoefs, aes(x=reorder(Variable,Coef), y=Coef, colour=Significance)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(alpha=0.75) +
  geom_point() +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0) +
  labs(y='Coefficient', x='', title='Breeding model A', subtitle='122 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','Significant' = '#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="none")

BBSresults1

### 3.1.3 Model averaging (just flocks) ----

model.set <- MuMIn::dredge(BBSfit, m.lim=c(0,nrow(BBS_trait_slopes_flocks2)/10))

top.models <- MuMIn::get.models(model.set, subset = delta <2.0)

averaged.model <- MuMIn::model.avg(top.models) 

#summary(averaged.model)
BBSAVGcoefs <- data.frame(averaged.model$coefficients)
BBSAVGcoefs$Avg <- rownames(BBSAVGcoefs)
BBSAVGcoefs <- melt(BBSAVGcoefs)
levels(BBSAVGcoefs$variable)[1] <-  '(Intercept)'

Conf <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                           level = 0.95,    # lets you adjust the confidence interval
                           full = T))    
Conf$variable <- rownames(Conf)
Conf$Avg <- 'full' 

Conf2 <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                            level = 0.95,    # lets you adjust the confidence interval
                            full = F))    
Conf2$variable <- rownames(Conf2)
Conf2$Avg <- 'subset'

BBSAVGcoefs <- join(BBSAVGcoefs, Conf)

BBSAVGcoefs[BBSAVGcoefs$Avg == 'subset',]$X2.5.. <- Conf2$X2.5..
BBSAVGcoefs[BBSAVGcoefs$Avg == 'subset',]$X97.5.. <- Conf2$X97.5..

BBSAVGcoefs <- BBSAVGcoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  variable == '(Intercept)' ~ 'Intercept',
  variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  variable == 'Mig_dist_numeric' ~ 'Migratory distance',
  variable == 'FlocksNo' ~ 'Solo',
  variable == 'FlocksYes' ~ 'Flocking\n(ref: solo)',
  variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  variable == 'GenFLS' ~ 'Generation length',
  variable == 'Trend' ~ 'Absolute population trend',
  variable == 'Breeding_rangesize' ~ 'Total range size',
  variable == 'Breeding_overlap_size' ~ 'Overlap range size',
  variable == 'Combined_Migration_Timingnight' ~ 'Migration timing\n(ref: day)'
))

BBSAVGplot <- ggplot(BBSAVGcoefs, aes(x=reorder(Variable, value), y=value, colour=Significance, shape = Avg)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(posiiton=position_dodge(), position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.25)) +
  labs(y='Coefficient', x='', shape='Average', title='Breeding model A', subtitle='122 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="right")

BBSAVGplot

BBSAVGcoefs_flocks <- BBSAVGcoefs
BBSAVGcoefs_flocks$n_mod <- length(top.models)


### 3.2 BBS (FlockWithAdults) ----

### 3.2.1 Full model ----
species_list <- unique(CBC_trait_slopes$CBC_Tree_species)
species_list <- sub(' ', '_',species_list)

BBS_trait_slopes <- BBS_trait_slopes %>% filter(AOU %in% CBC_trait_slopes$AOU)

# Load bird supertree based on Hackett's backbone 
bird_tree_hackett <- read.tree("Data/Hackett.tre")
bird_tree_hackett <- collapse.singles(bird_tree_hackett)

# Get only bird species from the supertree that are also included in collected data
# Extract list of all species from the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)
bird_tree_species <- data.frame(Name = bird_tree_species)

bird_tree_species2 <- bird_tree_species %>% tidyr::separate(Name, c('Genus', 'Species'))
bird_tree_species$Genus <- bird_tree_species2$Genus
bird_tree_species$Species <- bird_tree_species2$Species

# Check the overlap of species names between collected data file and the supertree
# All species should be present. If not, species names may not match with names in the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)

#intersect(bird_tree_species, species_list)
#setdiff(BBS_trait_slopes$BBS_Tree_species,bird_tree_species)

# Prune supertree to the list of taxa included in the data
pruned_birds_stree <- drop.tip(bird_tree_hackett, bird_tree_hackett$tip.label[-match(BBS_trait_slopes$BBS_Tree_species, bird_tree_hackett$tip.label)])

#is.binary(pruned_birds_stree) # TRUE
#is.ultrametric(pruned_birds_stree) # TRUE

treedat <- data.frame(BBS_Tree_species=pruned_birds_stree$tip.label)
treedat <- join(treedat, BBS_trait_slopes[c('BBS_Tree_species','FlocksWithAdults')])
names(treedat)[1] <- 'label'
pruned_birds_stree2 <- full_join(pruned_birds_stree,treedat, by='label')

pruned_birds_stree3 <- pruned_birds_stree2
pruned_birds_stree3@phylo$tip.label <- gsub('_',' ',pruned_birds_stree3@phylo$tip.label)

SE <- BBS_trait_slopes$BBS_dist_error
SE<-setNames(SE,BBS_trait_slopes$BBS_Tree_species)
BBS_trait_slopes$SE <- SE

BBS_trait_slopes <- BBS_trait_slopes[match(pruned_birds_stree$tip.label, BBS_trait_slopes$BBS_Tree_species),]

rownames(BBS_trait_slopes) <- BBS_trait_slopes$BBS_Tree_species

BBS_trait_slopes2 <- BBS_trait_slopes %>% mutate(
  Habitat_specialism_score = scale(Habitat_specialism_score),
  Diet_specialism_score = scale(Diet_specialism_score),
  GenFLS = scale(GenFLS),
  Trend = scale(abs(Trend)),
  Relative.Abundance = scale(Relative.Abundance),
  Breeding_rangesize = scale(Breeding_rangesize),
  Breeding_overlap_size = scale(Breeding_overlap_size),
  Mig_dist_numeric = scale(Mig_dist_numeric)
)
BBS_trait_slopes2$FlocksWithAdults <- factor(BBS_trait_slopes2$FlocksWithAdults)
BBS_trait_slopes2$FlocksWithAdults <- relevel(BBS_trait_slopes2$FlocksWithAdults, ref='Solo')
BBSfit<-pgls.SEy((BBS_dist) ~ Migration
                 + Mig_dist_numeric
                 + FlocksWithAdults
                 + Habitat_specialism_score
                 + Diet_specialism_score
                 + GenFLS
                 + Trend
                 + Breeding_rangesize
                 + Breeding_overlap_size
                 + Combined_Migration_Timing
                 ,data=BBS_trait_slopes2,se=(SE),tree=pruned_birds_stree,method="ML")
#summary(BBSfit)
car::vif(BBSfit)
#anova(BBSfit)
#coef(BBSfit)
#MuMIn::dredge(BBSfit)

BBScoefs <- data.frame(confint(BBSfit))
BBScoefs$Coef <- coef(BBSfit)
BBScoefs$Variable <- rownames(BBScoefs)
BBScoefs <- BBScoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  Variable == '(Intercept)' ~ 'Intercept',
  Variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  Variable == 'Mig_dist_numeric' ~ 'Migratory distance',
  Variable == 'FlocksWithAdultsFlocksNoAdults' ~ 'Age-separated flocks\n(ref: solo)',
  Variable == 'FlocksWithAdultsFlocksWithAdults' ~ 'Mixed-age flocks\n(ref: solo)',
  Variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  Variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  Variable == 'GenFLS' ~ 'Generation length',
  Variable == 'Trend' ~ 'Absolute population trend',
  Variable == 'Breeding_rangesize' ~ 'Total range size',
  Variable == 'Breeding_overlap_size' ~ 'Overlap range size',
  Variable == 'Combined_Migration_Timingnight' ~ 'Migration timing\n(ref: day)'
))

BBSresults2 <- ggplot(BBScoefs, aes(x=reorder(Variable,Coef), y=Coef, colour=Significance)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(alpha=0.75) +
  geom_point() +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0) +
  labs(y='Coefficient', x='', title='Breeding model B', subtitle='81 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="none")

BBSresults2
BBScoefs_FWA <- BBScoefs

### 3.2.2 Backwards stepwise deletion ----

BBSfitBSD<-pgls.SEy((BBS_dist) ~ FlocksWithAdults
                                          + Habitat_specialism_score
                                          + Trend
                                          ,data=BBS_trait_slopes2,se=(SE),tree=pruned_birds_stree,method="ML")

#summary(BBSfitBSD)
#car::vif(BBSfitBSD)
anova(BBSfitBSD)

BBScoefs <- data.frame(confint(BBSfitBSD))
BBScoefs$Coef <- coef(BBSfitBSD)
BBScoefs$Variable <- rownames(BBScoefs)
BBScoefs <- BBScoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  Variable == '(Intercept)' ~ 'Intercept',
  Variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  Variable == 'Mig_distanceshort' ~ 'Short-distance migrant\n(ref: long-distance)',
  Variable == 'FlocksWithAdultsFlocksNoAdults' ~ 'Age-separated flocks\n(ref: solo)',
  Variable == 'FlocksWithAdultsFlocksWithAdults' ~ 'Mixed-age flocks\n(ref: solo)',
  Variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  Variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  Variable == 'GenFLS' ~ 'Generation length',
  Variable == 'Trend' ~ 'Absolute population trend'
))

BBSresults2 <- ggplot(BBScoefs, aes(x=reorder(Variable,Coef), y=Coef, colour=Significance)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(alpha=0.75) +
  geom_point() +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0) +
  labs(y='Coefficient', x='', title='Breeding model B', subtitle='81 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','Significant' = '#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="none")

BBSresults2

### 3.2.3 Model averaging ----

model.set <- MuMIn::dredge(BBSfit, m.lim=c(0,nrow(BBS_trait_slopes2)/10))

top.models <- MuMIn::get.models(model.set, subset = delta <2.0)

averaged.model <- MuMIn::model.avg(top.models) 

#summary(averaged.model)
BBSAVGcoefs <- data.frame(averaged.model$coefficients)
BBSAVGcoefs$Avg <- rownames(BBSAVGcoefs)
BBSAVGcoefs <- melt(BBSAVGcoefs)
levels(BBSAVGcoefs$variable)[1] <-  '(Intercept)'

Conf <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                           level = 0.95,    # lets you adjust the confidence interval
                           full = T))    
Conf$variable <- rownames(Conf)
Conf$Avg <- 'full' 

Conf2 <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                            level = 0.95,    # lets you adjust the confidence interval
                            full = F))    
Conf2$variable <- rownames(Conf2)
Conf2$Avg <- 'subset'

BBSAVGcoefs <- join(BBSAVGcoefs, Conf)

BBSAVGcoefs[BBSAVGcoefs$Avg == 'subset',]$X2.5.. <- Conf2$X2.5..
BBSAVGcoefs[BBSAVGcoefs$Avg == 'subset',]$X97.5.. <- Conf2$X97.5..

BBSAVGcoefs <- BBSAVGcoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  variable == '(Intercept)' ~ 'Intercept',
  variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  variable == 'Mig_dist_numeric' ~ 'Migratory distance',
  variable == 'FlocksWithAdultsFlocksNoAdults' ~ 'Age-separated flocks\n(ref: solo)',
  variable == 'FlocksWithAdultsFlocksWithAdults' ~ 'Mixed-age flocks\n(ref: solo)',
  variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  variable == 'GenFLS' ~ 'Generation length',
  variable == 'Trend' ~ 'Absolute population trend',
  variable == 'Breeding_rangesize' ~ 'Total range size',
  variable == 'Breeding_overlap_size' ~ 'Overlap range size',
  variable == 'Combined_Migration_Timingnight' ~ 'Migration timing\n(ref: day)'
))

BBSAVGplot <- ggplot(BBSAVGcoefs, aes(x=reorder(Variable, value), y=value, colour=Significance, shape = Avg)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(posiiton=position_dodge(), position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.25)) +
  labs(y='Coefficient', x='', shape='Average', title='Breeding model B', subtitle='81 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="right")

BBSAVGplot

BBSAVGcoefs_fwa <- BBSAVGcoefs
BBSAVGcoefs_fwa$n_mod <- length(top.models)

## 3.3 CBC (just flocks) ----

### 3.3.1 CBC full model ----
species_list <- unique(CBC_trait_slopes_flocks$CBC_Tree_species)
species_list <- sub(' ', '_',species_list)

CBC_trait_slopes_flocks <- CBC_trait_slopes_flocks %>% filter(AOU %in% CBC_trait_slopes_flocks$AOU)

# Load bird supertree based on Hackett's backbone 
bird_tree_hackett <- read.tree("Data/Hackett.tre")

bird_tree_hackett <- collapse.singles(bird_tree_hackett)

# Get only bird species from the supertree that are also included in collected data
# Extract list of all species from the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)
bird_tree_species <- data.frame(Name = bird_tree_species)

bird_tree_species2 <- bird_tree_species %>% tidyr::separate(Name, c('Genus', 'Species'))
bird_tree_species$Genus <- bird_tree_species2$Genus
bird_tree_species$Species <- bird_tree_species2$Species

# Check the overlap of species names between collected data file and the supertree
# All species should be present. If not, species names may not match with names in the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)

#intersect(bird_tree_species, species_list)
#setdiff(CBC_trait_slopes$CBC_Tree_species,bird_tree_species)

# Prune supertree to the list of taxa included in the data
pruned_birds_stree <- drop.tip(bird_tree_hackett, bird_tree_hackett$tip.label[-match(CBC_trait_slopes_flocks$CBC_Tree_species, bird_tree_hackett$tip.label)])

#is.binary(pruned_birds_stree) # TRUE
#is.ultrametric(pruned_birds_stree) # TRUE

treedat <- data.frame(CBC_Tree_species=pruned_birds_stree$tip.label)
treedat <- join(treedat, CBC_trait_slopes_flocks[c('CBC_Tree_species','Flocks')])
names(treedat)[1] <- 'label'
pruned_birds_stree2 <- full_join(pruned_birds_stree,treedat, by='label')

pruned_birds_stree3 <- pruned_birds_stree2
pruned_birds_stree3@phylo$tip.label <- gsub('_',' ',pruned_birds_stree3@phylo$tip.label)

SE <- CBC_trait_slopes_flocks$CBC_dist_error
SE<-setNames(SE,CBC_trait_slopes_flocks$CBC_Tree_species)
CBC_trait_slopes_flocks$SE <- SE

CBC_trait_slopes_flocks <- CBC_trait_slopes_flocks[match(pruned_birds_stree$tip.label, CBC_trait_slopes_flocks$CBC_Tree_species),]

rownames(CBC_trait_slopes_flocks) <- CBC_trait_slopes_flocks$CBC_Tree_species

CBC_trait_slopes_flocks2 <- CBC_trait_slopes_flocks %>% mutate(
  Habitat_specialism_score = scale(Habitat_specialism_score),
  Diet_specialism_score = scale(Diet_specialism_score),
  GenFLS = scale(GenFLS),
  Trend = scale(abs(Trend)),
  Relative.Abundance = scale(Relative.Abundance),
  Nonbreeding_overlap_size = scale(Nonbreeding_overlap_size),
  Nonbreeding_rangesize = scale(Nonbreeding_rangesize),
  Mig_dist_numeric = scale(Mig_dist_numeric)
)

CBCfit<-pgls.SEy((CBC_dist) ~ Migration
                 + Mig_dist_numeric
                 + Flocks
                 + Habitat_specialism_score
                 + Diet_specialism_score
                 + GenFLS
                 + Trend
                 + Nonbreeding_rangesize
                 + Nonbreeding_overlap_size
                 + Combined_Migration_Timing
                 ,data=CBC_trait_slopes_flocks2,se=(SE),tree=pruned_birds_stree,method="ML")
#summary(CBCfit)
car::vif(CBCfit)
#anova(CBCfit)
#coef(CBCfit)
#MuMIn::dredge(CBCfit)

CBCcoefs <- data.frame(confint(CBCfit))
CBCcoefs$Coef <- coef(CBCfit)
CBCcoefs$Variable <- rownames(CBCcoefs)
CBCcoefs <- CBCcoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  Variable == '(Intercept)' ~ 'Intercept',
  Variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  Variable == 'Mig_dist_numeric' ~ 'Migratory distance',
  Variable == 'FlocksNo' ~ 'Solo',
  Variable == 'FlocksYes' ~ 'Flocking\n(ref: solo)',
  Variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  Variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  Variable == 'GenFLS' ~ 'Generation length',
  Variable == 'Trend' ~ 'Absolute population trend',
  Variable == 'Nonbreeding_rangesize' ~ 'Total range size',
  Variable == 'Nonbreeding_overlap_size' ~ 'Overlap range size',
  Variable == 'Combined_Migration_Timingnight' ~ 'Migration timing\n(ref: day)'
))

NonBreedingresults1 <- ggplot(CBCcoefs, aes(x=reorder(Variable,Coef), y=Coef, colour=Significance)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(alpha=0.75) +
  geom_point() +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0) +
  labs(y='Coefficient', x='', title='Non-breeding model A', subtitle='122 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="none")

NonBreedingresults1

CBCcoefs_flocks <- CBCcoefs

### 3.3.2 Backwards stepwise deletion ----

CBCfitBSD<-pgls.SEy((CBC_dist) ~ Flocks
                    + Diet_specialism_score
                    ,data=CBC_trait_slopes_flocks2,se=(SE),tree=pruned_birds_stree,method="ML")
#summary(CBCfitBSD)
#car::vif(CBCfitBSD)
#anova(CBCfitBSD)

CBCcoefs <- data.frame(confint(CBCfitBSD))
CBCcoefs$Coef <- coef(CBCfitBSD)
CBCcoefs$Variable <- rownames(CBCcoefs)
CBCcoefs <- CBCcoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  Variable == '(Intercept)' ~ 'Intercept',
  Variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  Variable == 'Mig_distanceshort' ~ 'Short-distance migrant\n(ref: long-distance)',
  Variable == 'FlocksNo' ~ 'Solo',
  Variable == 'FlocksYes' ~ 'Flocking\n(ref: solo)',
  Variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  Variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  Variable == 'GenFLS' ~ 'Generation length',
  Variable == 'Trend' ~ 'Absolute population trend'
))

NonBreedingresults1 <- ggplot(CBCcoefs, aes(x=reorder(Variable,Coef), y=Coef, colour=Significance)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(alpha=0.75) +
  geom_point() +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0) +
  labs(y='Coefficient', x='', title='Non-breeding model A', subtitle='122 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','Significant' = '#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="none")

NonBreedingresults1

### 3.3.3 Model averaging (just flocks) ----

model.set <- MuMIn::dredge(CBCfit, m.lim=c(0,nrow(CBC_trait_slopes_flocks2)/10))

top.models <- MuMIn::get.models(model.set, subset = delta <2.0)

averaged.model <- MuMIn::model.avg(top.models) 

#summary(averaged.model)
CBCAVGcoefs <- data.frame(averaged.model$coefficients)
CBCAVGcoefs$Avg <- rownames(CBCAVGcoefs)
CBCAVGcoefs <- melt(CBCAVGcoefs)
levels(CBCAVGcoefs$variable)[1] <-  '(Intercept)'

Conf <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                           level = 0.95,    # lets you adjust the confidence interval
                           full = T))    
Conf$variable <- rownames(Conf)
Conf$Avg <- 'full' 

Conf2 <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                            level = 0.95,    # lets you adjust the confidence interval
                            full = F))    
Conf2$variable <- rownames(Conf2)
Conf2$Avg <- 'subset'

CBCAVGcoefs <- join(CBCAVGcoefs, Conf)

CBCAVGcoefs[CBCAVGcoefs$Avg == 'subset',]$X2.5.. <- Conf2$X2.5..
CBCAVGcoefs[CBCAVGcoefs$Avg == 'subset',]$X97.5.. <- Conf2$X97.5..

CBCAVGcoefs <- CBCAVGcoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  variable == '(Intercept)' ~ 'Intercept',
  variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  variable == 'Mig_dist_numeric' ~ 'Migratory distance',
  variable == 'FlocksNo' ~ 'Solo',
  variable == 'FlocksYes' ~ 'Flocking\n(ref: solo)',
  variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  variable == 'GenFLS' ~ 'Generation length',
  variable == 'Trend' ~ 'Absolute population trend',
  variable == 'Nonbreeding_rangesize' ~ 'Total range size',
  variable == 'Nonbreeding_overlap_size' ~ 'Overlap range size',
  variable == 'Combined_Migration_Timingnight' ~ 'Migration timing\n(ref: day)'
))

CBCAVGplot <- ggplot(CBCAVGcoefs, aes(x=reorder(Variable, value), y=value, colour=Significance, shape = Avg)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(posiiton=position_dodge(), position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.25)) +
  labs(y='Coefficient', x='', shape='Average', title='Non-breeding model A', subtitle='122 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="right")

CBCAVGplot

CBCAVGcoefs_flocks <- CBCAVGcoefs
CBCAVGcoefs_flocks$n_mod <- length(top.models)

### 3.4 CBC (FlockWithAdults) ----

# 3.4.1 Full model ----

species_list <- unique(CBC_trait_slopes$CBC_Tree_species)
species_list <- sub(' ', '_',species_list)

# Check the overlap of species names between collected data file and the supertree
# All species should be present. If not, species names may not match with names in the supertree

#intersect(bird_tree_species, species_list)

#setdiff(CBC_trait_slopes$CBC_Tree_species,bird_tree_species)

# Prune supertree to the list of taxa included in the data
CBC_pruned_birds_stree <- drop.tip(bird_tree_hackett, bird_tree_hackett$tip.label[-na.omit(match(CBC_trait_slopes$CBC_Tree_species, bird_tree_hackett$tip.label))])

#is.binary(pruned_birds_stree) # TRUE
#is.ultrametric(pruned_birds_stree) # TRUE

treedat <- data.frame(CBC_Tree_species=CBC_pruned_birds_stree$tip.label)
treedat <- join(treedat, CBC_trait_slopes[c('CBC_Tree_species','FlocksWithAdults','Order')])
names(treedat)[1] <- 'label'
pruned_birds_stree2 <- full_join(CBC_pruned_birds_stree,treedat, by='label')

pruned_birds_stree3 <- pruned_birds_stree2
pruned_birds_stree3@phylo$tip.label <- gsub('_',' ',pruned_birds_stree3@phylo$tip.label)

CBC_SE <- CBC_trait_slopes$CBC_dist_error
CBC_SE<-setNames(CBC_SE,CBC_trait_slopes$CBC_Tree_species)
CBC_trait_slopes$SE <- CBC_SE

CBC_trait_slopes <- CBC_trait_slopes[match(CBC_pruned_birds_stree$tip.label, CBC_trait_slopes$CBC_Tree_species),]

rownames(CBC_trait_slopes) <- CBC_trait_slopes$CBC_Tree_species

CBC_trait_slopes2 <- CBC_trait_slopes %>% mutate(
  Habitat_specialism_score = scale(Habitat_specialism_score),
  Diet_specialism_score = scale(Diet_specialism_score),
  GenFLS = scale(GenFLS),
  Trend = scale(abs(Trend)),
  Relative.Abundance = scale(Relative.Abundance),
  Mig_dist_numeric = scale(Mig_dist_numeric),
  Nonbreeding_overlap_size = scale(Nonbreeding_overlap_size),
  Nonbreeding_rangesize = scale(Nonbreeding_rangesize)
)

CBC_trait_slopes2$FlocksWithAdults <- factor(CBC_trait_slopes2$FlocksWithAdults)
CBC_trait_slopes2$FlocksWithAdults <- relevel(CBC_trait_slopes2$FlocksWithAdults, ref='Solo')

CBCfit<-pgls.SEy((CBC_dist) ~ Migration
                 + Mig_dist_numeric
                 + FlocksWithAdults
                 + Habitat_specialism_score
                 + Diet_specialism_score
                 + GenFLS
                 #+ Trend
                 #+ Nonbreeding_rangesize
                 + Nonbreeding_overlap_size
                 + Combined_Migration_Timing
                 ,data=CBC_trait_slopes2,se=CBC_SE,tree=CBC_pruned_birds_stree,method="ML")

#summary(CBCfit)

car::vif(CBCfit)
#anova(CBCfit)
#coef(CBCfit)
#confint(CBCfit)
#MuMIn::dredge(CBCfit)
#plot(CBCfit)

CBCcoefs <- data.frame(confint(CBCfit))
CBCcoefs$Coef <- coef(CBCfit)
CBCcoefs$Variable <- rownames(CBCcoefs)
CBCcoefs <- CBCcoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  Variable == '(Intercept)' ~ 'Intercept',
  Variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  Variable == 'Mig_dist_numeric' ~ 'Migratory distance',
  Variable == 'FlocksWithAdultsFlocksNoAdults' ~ 'Age-separated flocks\n(ref: solo)',
  Variable == 'FlocksWithAdultsFlocksWithAdults' ~ 'Mixed-age flocks\n(ref: solo)',
  Variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  Variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  Variable == 'GenFLS' ~ 'Generation length',
  Variable == 'Trend' ~ 'Absolute population trend',
  Variable == 'Nonbreeding_rangesize' ~ 'Total range size',
  Variable == 'Nonbreeding_overlap_size' ~ 'Overlap range size',
  Variable == 'Combined_Migration_Timingnight' ~ 'Migration timing\n(ref: day)'
))

NonBreedingresults2 <- ggplot(CBCcoefs, aes(x=reorder(Variable,Coef), y=Coef, colour=Significance)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point() +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0) +
  labs(y='Coefficient', x='', title='Non-breeding model B', subtitle='81 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="none")

NonBreedingresults2

test_results <- data.frame(emmeans::emmeans(CBCfit, specs = "FlocksWithAdults", plotIt=FALSE))
test_results <- test_results %>% mutate(FlocksWithAdults = case_when(FlocksWithAdults == 'FlocksWithAdults' ~ 'Mixed age flocks',
                                                                     FlocksWithAdults == 'Solo' ~ 'Solo',
                                                                     FlocksWithAdults == 'FlocksNoAdults' ~ 'Age separated flocks'))
test_results$FlocksWithAdults <- as.factor(test_results$FlocksWithAdults)
test_results$FlocksWithAdults <- relevel(test_results$FlocksWithAdults, ref='Solo')

CBCcoefs_FWA <- CBCcoefs

### 3.4.2 Backwards stepwise deletion ---- 

CBCfitBSD<-pgls.SEy((CBC_dist) ~ FlocksWithAdults
                    ,data=CBC_trait_slopes2,se=CBC_SE,tree=CBC_pruned_birds_stree,method="ML")

#summary(CBCfitBSD)
#car::vif(CBCfitBSD)
#anova(CBCfitBSD)

CBCcoefs <- data.frame(confint(CBCfitBSD))
CBCcoefs$Coef <- coef(CBCfitBSD)
CBCcoefs$Variable <- rownames(CBCcoefs)
CBCcoefs <- CBCcoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  Variable == '(Intercept)' ~ 'Intercept',
  Variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  Variable == 'Mig_distanceshort' ~ 'Short-distance migrant\n(ref: long-distance)',
  Variable == 'FlocksWithAdultsFlocksNoAdults' ~ 'Age-separated flocks\n(ref: solo)',
  Variable == 'FlocksWithAdultsFlocksWithAdults' ~ 'Mixed-age flocks\n(ref: solo)',
  Variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  Variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  Variable == 'GenFLS' ~ 'Generation length',
  Variable == 'Trend' ~ 'Absolute population trend'
))

NonBreedingresults2 <- ggplot(CBCcoefs, aes(x=reorder(Variable,Coef), y=Coef, colour=Significance)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point() +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0) +
  labs(y='Coefficient', x='', title='Non-breeding model B', subtitle='81 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','Significant' = '#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="none")

NonBreedingresults2

### 3.4.3 Model averaging ----

model.set <- MuMIn::dredge(CBCfit, m.lim=c(0,nrow(CBC_trait_slopes2)/10))

top.models <- MuMIn::get.models(model.set, subset = delta <2.0)

averaged.model <- MuMIn::model.avg(top.models) 

#summary(averaged.model)
CBCAVGcoefs <- data.frame(averaged.model$coefficients)
CBCAVGcoefs$Avg <- rownames(CBCAVGcoefs)
CBCAVGcoefs <- melt(CBCAVGcoefs)
levels(CBCAVGcoefs$variable)[1] <-  '(Intercept)'

Conf <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                           level = 0.95,    # lets you adjust the confidence interval
                           full = T))    
Conf$variable <- rownames(Conf)
Conf$Avg <- 'full' 

Conf2 <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                            level = 0.95,    # lets you adjust the confidence interval
                            full = F))    
Conf2$variable <- rownames(Conf2)
Conf2$Avg <- 'subset'

CBCAVGcoefs <- join(CBCAVGcoefs, Conf)

CBCAVGcoefs[CBCAVGcoefs$Avg == 'subset',]$X2.5.. <- Conf2$X2.5..
CBCAVGcoefs[CBCAVGcoefs$Avg == 'subset',]$X97.5.. <- Conf2$X97.5..

CBCAVGcoefs <- CBCAVGcoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  variable == '(Intercept)' ~ 'Intercept',
  variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  variable == 'Mig_dist_numeric' ~ 'Migratory distance',
  variable == 'FlocksWithAdultsFlocksNoAdults' ~ 'Age-separated flocks\n(ref: solo)',
  variable == 'FlocksWithAdultsFlocksWithAdults' ~ 'Mixed-age flocks\n(ref: solo)',
  variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  variable == 'GenFLS' ~ 'Generation length',
  variable == 'Trend' ~ 'Absolute population trend',
  variable == 'Nonbreeding_rangesize' ~ 'Total range size',
  variable == 'Nonbreeding_overlap_size' ~ 'Overlap range size',
  variable == 'Combined_Migration_Timingnight' ~ 'Migration timing\n(ref: day)'
))

CBCAVGplot <- ggplot(CBCAVGcoefs, aes(x=reorder(Variable, value), y=value, colour=Significance, shape = Avg)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(posiiton=position_dodge(), position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.25)) +
  labs(y='Coefficient', x='', shape='Average', title='Non-breeding model B', subtitle='81 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="right")

CBCAVGplot
CBCAVGcoefs_fwa <- CBCAVGcoefs
CBCAVGcoefs_fwa$n_mod <- length(top.models)

test_results <- data.frame(emmeans::emmeans(averaged.model, specs = "FlocksWithAdults", plotIt=FALSE, data=CBC_trait_slopes2))
test_results <- test_results %>% mutate(FlocksWithAdults = case_when(FlocksWithAdults == 'FlocksWithAdults' ~ 'Mixed age flocks',
                                                                     FlocksWithAdults == 'Solo' ~ 'Solo',
                                                                     FlocksWithAdults == 'FlocksNoAdults' ~ 'Age separated flocks'))
test_results$FlocksWithAdults <- as.factor(test_results$FlocksWithAdults)
test_results$FlocksWithAdults <- relevel(test_results$FlocksWithAdults, ref='Solo')


BBSAVGcoefs_flocks
BBSAVGcoefs_fwa
CBCAVGcoefs_flocks
CBCAVGcoefs_fwa

# 4. Manuscript figures ----

CBC_movements <- data.frame(AOU=rep(CBC_trait_slopes2$AOU,each=2), Lat=NA,Long=NA)
setwd("~/Desktop/MIG STRAT MS sub/Data/COAs")

for(i in unique(CBC_movements$AOU)){
  
  CBC_COA <- read.csv(paste0(i,'_CBC_COAs.csv'))
  
  ## we need to change from EPSG 102010 to WGS84 for easy plotting
  
  CBC_COA_1 <- st_as_sf(CBC_COA, coords = c('coa_x', 'coa_y'))
  
  st_crs(CBC_COA_1) <- st_crs(centroids)

  CBC_COA_2 <- st_transform(CBC_COA_1, crs='WGS84')
  
  coords <- data.frame(st_coordinates(CBC_COA_2$geometry))
  CBC_COA$coa_x <- coords$X
  CBC_COA$coa_y <- coords$Y

  CBCmod_lat <- lm(coa_y ~ year, data=CBC_COA)
  CBCmod_lon <- lm(coa_x ~ year, data=CBC_COA)
  
  CBC_pred_lat <- ggpredict(CBCmod_lat, terms = c('year[1970,2019]'))
  CBC_pred_lon <- ggpredict(CBCmod_lon, terms = c('year[1970,2019]'))
  
  CBC_movements[CBC_movements$AOU==i,]$Lat <- CBC_pred_lat$predicted
  CBC_movements[CBC_movements$AOU==i,]$Long <- CBC_pred_lon$predicted
  
}

CBC_movements2 <- join(CBC_movements, CBC_trait_slopes2[c('AOU','FlocksWithAdults','CBC_dist_error')])
CBC_movements2$CBC_dist_error <- 1/(CBC_movements2$CBC_dist_error)

CBC_trait_slopes2 <- CBC_trait_slopes2 %>% mutate(FlocksWithAdults = case_when(FlocksWithAdults == 'FlocksWithAdults' ~ 'Mixed age flocks',
                                                                               FlocksWithAdults == 'Solo' ~ 'Solo',
                                                                               FlocksWithAdults == 'FlocksNoAdults' ~ 'Age separated flocks'))
CBC_trait_slopes2$FlocksWithAdults <- as.factor(CBC_trait_slopes2$FlocksWithAdults)

## Figure 1 ----

centroids <- st_read("~/Desktop/MIG STRAT MS sub/Data/analytical_strata[78]/bcr_by_state_clean_dissolved_centroid_points_102010.shp", quiet=T)
world2 <- ne_countries(scale = "large", continent = 'north america',returnclass = "sf")
lakes <- st_read('~/Desktop/MIG STRAT MS sub/Data/ne_50m_lakes/ne_50m_lakes.shp')
st_crs(world2) <- st_crs(centroids)
st_crs(lakes) <- st_crs(centroids)

CBC_movements3 <- CBC_movements2 %>% mutate(FlocksWithAdults = case_when(FlocksWithAdults == 'FlocksWithAdults' ~ 'Mixed-age flocks',
                                                                         FlocksWithAdults == 'Solo' ~ 'Solo',
                                                                         FlocksWithAdults == 'FlocksNoAdults' ~ 'Age-separated flocks'))

cols <- c("Solo" = "#ED6B5A", "Age-separated flocks" = "#3083DC", "Mixed-age flocks" = "#1F008F")

crs_string <- "+proj=ortho +lat_0=20 +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs"
world3 <- world2 %>% st_transform(crs = crs_string)
world3 <- map_data('world')



polys5 <- st_transform(polys3, 4326)

polys5 <- tibble(lon = list(st_coordinates(polys5)[, 1]),
            lat = list(st_coordinates(polys5)[, 2])) %>% 
  unnest(c(lon, lat)) %>% 
  st_drop_geometry()
polys6 <- lidR::concaveman(as.matrix(polys5))
hull <- polys5 %>%
  slice(chull(lon,lat))

mapplot <- ggplot() +  
  geom_polygon(data=world3, aes(long, lat, group=group), fill='grey85', size=0.2, colour='grey65') +
  coord_map('ortho',  orientation = c(-0,-100,0),xlim = c(-118,-79), ylim = c(25,45)) +
  geom_path(CBC_movements3, mapping=aes(x=Long, y=Lat, group=AOU, colour=FlocksWithAdults,alpha=CBC_dist_error),linewidth=.6,arrow=arrow(length=unit(0.1,'cm'), type='closed')) +
  #geom_text(CBC_movements3[seq(1,nrow(CBC_movements3), 2),], mapping=aes(x=Long, y=Lat, label=AOU), size=2) +
  theme_void() +
  theme(panel.background = element_blank(),panel.border = element_blank(), legend.position = c(0.5,0.075), legend.text = element_text(size=10), legend.title = element_text(size=10), plot.background = element_blank()) +
  scale_color_manual(values=cols, guide = 'none') +
  scale_alpha_continuous(range = c(0.2, 1), n.breaks=4, breaks=seq(0.001, 0.004, 0.001), labels=c('<0.001', seq(0.002, 0.004, 0.001))) +
  labs(fill = "", colour = "", alpha='Reciprocal error:') +
  guides(alpha = guide_legend(nrow=1,byrow=TRUE))
mapplot
cols2 <- c("FlocksWithAdults" = "#1F008F","Solo" = "#ED6B5A", "FlocksNoAdults" = "#3083DC")



effectsplot <- ggplot(test_results, aes(x=FlocksWithAdults, y=emmean/1000, colour=FlocksWithAdults)) +
  geom_errorbar(aes(ymin=(emmean-SE)/1000, ymax=(emmean+SE)/1000), width=0.1, size=1, position=position_dodge(width=0.75)) +
  geom_point(size=3, position=position_dodge(width=0.75)) +
  scale_color_manual(values=c("#ED6B5A", "#3083DC", "#1F008F"), guide=F) +
  theme_half_open() +
  scale_x_discrete(labels=c('Solo','Age-separated\nflocks','Mixed-age\nflocks')) +
  labs(y='Estimated annual shift in non-breeding\ncentre of abundance (km)') +
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size=7),axis.text.x = element_text(size=10), axis.title.y = element_text(size=10), legend.position = 'bottom', legend.title = element_blank())


treeplot <- ggtree::ggtree(pruned_birds_stree3, layout='circular', aes(colour=FlocksWithAdults)) + 
  geom_tippoint(aes(colour=FlocksWithAdults), size=0.2) + 
  geom_strip(taxa1='Dendrocygna autumnalis',taxa2= 'Anas americana', barsize=0.3, color='black', label="Anseriformes", offset=4, offset.text = 5, fontsize=2) +
  geom_strip('Dendroica townsendi', 'Myiarchus cinerascens', barsize=0.3, color='black', label="Passeriformes", offset=4,offset.text = 7, fontsize=2) +
  geom_strip('Falco sparverius', 'Falco peregrinus', barsize=0.3, color='black', label="Falconiformes", offset=4, offset.text = 55, fontsize=2) +
  geom_strip('Larus marinus', 'Tringa flavipes', barsize=0.3, color='black', label="Charadriiformes", offset=4,offset.text = 10, hjust=1, fontsize=2) +
  geom_strip(72, 72, barsize=0.3, color='black', label="Coraciiformes", offset=4,offset.text = 50, fontsize=2, hjuust=-0.1) +
  geom_strip(71, 71, barsize=0.3, color='black', label="Strigiformes", offset=4,offset.text = 45, fontsize=2) +
  geom_strip('Pandion haliaetus', 'Buteo regalis', barsize=0.3, color='black', label="Accipitriformes", offset=4,offset.text = 55, hjust=0.05, fontsize=2) +
  geom_strip('Zenaida asiatica', 'Zenaida asiatica', barsize=0.3, color='black', label="Columbiformes", offset=4,offset.text = 6, fontsize=2) +
  geom_strip('Plegadis falcinellus', 'Plegadis falcinellus', barsize=0.3, color='black', label="Pelecaniformes", offset=4,offset.text = 11,hjust=0.1, fontsize=2) +
  geom_strip('Rallus limicola', 'Grus canadensis', barsize=0.3, color='black', label="Gruiformes", offset=4,offset.text = 16, fontsize=2, hjust=0.2) +
  geom_strip('Archilochus alexandri', 'Selasphorus rufus', barsize=0.3, color='black', label="Apodiformes", offset=4,offset.text = 7, fontsize=2) +
  labs(color='') +
  scale_color_manual(values=cols2,labels=c('Age-separated\nflocks','Mixed-age\nflocks','Solo')) +
  theme(legend.text = element_text(size=10),legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10), legend.background = element_blank(), legend.position = c(0.5,0.19), legend.box = "horizontal") +
  guides(colour=guide_legend(nrow=1,byrow=TRUE))
c(1,0.5)

(main_fig <- ggdraw() +
    theme(plot.background = element_rect(fill="white", color = NA)) +
    #draw_plot(treeplot, -0.08,-0.02, 0.6,0.6) +
    draw_plot(treeplot, -0.2,-0.15, 0.78,0.85) +
    draw_plot(mapplot, 0, 0.55, 0.5,0.45) +
    draw_plot(effectsplot,0.5,0,0.5,.95) +
    draw_plot_label(
      c("A", "B", "C"),
      c(0.01, 0.01, 0.51),
      c(1, 0.55, 1),
      size = 10
    ))

ggsave(main_fig, file='Figure_1_ortho.pdf', width=180, height=100, units='mm', dpi=600)

(main_fig <- ggdraw() +
    theme(plot.background = element_rect(fill="white", color = NA)) +
    #draw_plot(treeplot, -0.08,-0.02, 0.6,0.6) +
    draw_plot(treeplot, -0.06,-0.1, 0.65,0.9) +
    draw_plot(mapplot,0, 0.57, 1,0.41) +
    draw_plot(effectsplot,0.52,0.1,0.48,0.45) +
    draw_plot_label(
      c("A", "B", "C"),
      c(0.02, 0.02, 0.51),
      c(1, 0.55, 0.55),
      size = 10
    ))

ggsave(main_fig, file='Figure_1_alt_ortho.pdf', width=180, height=180, units='mm', dpi=600)

# Figure 2 ----

BBSAVGcoefs_flocks$Model <- 'Model A'
BBSAVGcoefs_fwa$Model <- 'Model B'
CBCAVGcoefs_flocks$Model <- 'Model A'
CBCAVGcoefs_fwa$Model <- 'Model B'

BBSAVGcoefs_flocks$Season <- 'Breeding'
BBSAVGcoefs_fwa$Season <- 'Breeding'
CBCAVGcoefs_flocks$Season <- 'Non-breeding'
CBCAVGcoefs_fwa$Season <- 'Non-breeding'

AVGcoefs <- rbind(BBSAVGcoefs_flocks,
      BBSAVGcoefs_fwa,
      CBCAVGcoefs_flocks,
      CBCAVGcoefs_fwa)
AVGcoefs <- AVGcoefs %>% filter(Avg == 'full')

 reorder(AVGcoefs$Variable,AVGcoefs$value)
unique(AVGcoefs$Variable)

str(AVGcoefs)
AVGcoefs <- AVGcoefs %>% mutate(Variable = forcats::fct_relevel(Variable, 
                          'Intercept','Flocking\n(ref: solo)','Mixed-age flocks\n(ref: solo)','Age-separated flocks\n(ref: solo)','Absolute population trend','Migratory distance','Habitat specialism score','Total range size','Diet specialism score','Overlap range size','Generation length','Partial migrant\n(ref: complete)','Migration timing\n(ref: day)')) 

levels(AVGcoefs$Variable)[levels(AVGcoefs$Variable) == 'Overlap range size'] <- 'Sampled range size'
AVGcoefs[is.na(AVGcoefs$Variable),]$Variable <- 'Sampled range size'

results_plot <- ggplot(AVGcoefs, aes(x=Variable, y=value, colour=Significance, shape=Model)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.4) +
  geom_point(size=2,position=position_dodge(width=0.3)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.3)) +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(AVGcoefs$Variable))) +
  facet_grid(~ Season) +
  labs(y='Coefficient') +
  theme(axis.title.x = element_text(size=13), axis.text.y = element_text(size=10), legend.position="bottom", axis.title.y = element_blank(), strip.background = element_blank(), legend.title = element_blank(), legend.text = element_text(size=10), plot.background = element_rect(fill='white'))

results_plot

ggsave(results_plot, file='Figure_2.pdf', width=180, height=150, units='mm', dpi=600)


#### Sensitivity analysis - flocking metric ----

BBS_trait_slopes_flocks <- read.csv('Data/BBS_trait_slopes_flocks.csv')
CBC_trait_slopes_flocks <- read.csv('Data/CBC_trait_slopes_flocks.csv')



CBC_trait_slopes_flocks[is.na(CBC_trait_slopes_flocks$Combined_Migration_Timing),]$Combined_Migration_Timing <- 'night'
BBS_trait_slopes_flocks[is.na(BBS_trait_slopes_flocks$Combined_Migration_Timing),]$Combined_Migration_Timing <- 'night'

data <- read.csv('Data/MS_data.csv')
data <- data %>% dplyr::select(Genus, Species, Flock.size, Maximum.flock.size.while.travelling)

BBS_trait_slopes_flocks <- left_join(BBS_trait_slopes_flocks, data)


BBS_trait_slopes_flocks[BBS_trait_slopes_flocks$Genus == 'Recurvirostra' & BBS_trait_slopes_flocks$Species == 'americana',]$Flock.size <- 'Medium'
BBS_trait_slopes_flocks[BBS_trait_slopes_flocks$Genus == 'Hydroprogne' & BBS_trait_slopes_flocks$Species == 'caspia',]$Flock.size <- 'Medium'
BBS_trait_slopes_flocks[BBS_trait_slopes_flocks$Genus == 'Tringa' & BBS_trait_slopes_flocks$Species == 'melanoleuca',]$Flock.size <- 'Medium'
BBS_trait_slopes_flocks[BBS_trait_slopes_flocks$Genus == 'Porzana' & BBS_trait_slopes_flocks$Species == 'carolina',]$Flock.size <- 'Medium'
BBS_trait_slopes_flocks[BBS_trait_slopes_flocks$Genus == 'Tyrannus' & BBS_trait_slopes_flocks$Species == 'verticalis',]$Flock.size <- 'Medium'
BBS_trait_slopes_flocks[BBS_trait_slopes_flocks$Genus == 'Bucephala' & BBS_trait_slopes_flocks$Species == 'albeola',]$Flock.size <- 'Small'




#BBS_trait_slopes_flocks <- BBS_trait_slopes_flocks %>% mutate(Flock_size = case_when(
#  Flocks == 'No' ~ 'Solo',
#  Flock.size == 'Small' ~ 'Small',
#  Flock.size %in% c('Medium', 'Large') ~ 'Large',
#  Flock.size == '' & Maximum.flock.size.while.travelling < 11 ~ 'Small',
#  Flock.size == '' & Maximum.flock.size.while.travelling > 10 ~ 'Large'
#))

BBS_trait_slopes_flocks <- BBS_trait_slopes_flocks %>% mutate(Flock_size = case_when(
  Flocks == 'No' ~ 'Solo',
  Flock.size == 'Small' ~ 'Small flocks',
  Flock.size %in% c('Medium', 'Large') ~ 'Large flocks',
  Flock.size == '' ~ NA,
  Flock.size == '' ~ NA
))

CBC_trait_slopes_flocks <- left_join(CBC_trait_slopes_flocks, data)
#CBC_trait_slopes_flocks <- CBC_trait_slopes_flocks %>% mutate(Flock_size = case_when(
#  Flocks == 'No' ~ 'Solo',
#  Flock.size == 'Small' ~ 'Small',
#  Flock.size %in% c('Medium', 'Large') ~ 'Large',
#  Flock.size == '' & Maximum.flock.size.while.travelling < 11 ~ 'Small',
#  Flock.size == '' & Maximum.flock.size.while.travelling > 10 ~ 'Large'
#))

CBC_trait_slopes_flocks[CBC_trait_slopes_flocks$Genus == 'Recurvirostra' & CBC_trait_slopes_flocks$Species == 'americana',]$Flock.size <- 'Medium'
CBC_trait_slopes_flocks[CBC_trait_slopes_flocks$Genus == 'Hydroprogne' & CBC_trait_slopes_flocks$Species == 'caspia',]$Flock.size <- 'Medium'
CBC_trait_slopes_flocks[CBC_trait_slopes_flocks$Genus == 'Tringa' & CBC_trait_slopes_flocks$Species == 'melanoleuca',]$Flock.size <- 'Medium'
CBC_trait_slopes_flocks[CBC_trait_slopes_flocks$Genus == 'Porzana' & CBC_trait_slopes_flocks$Species == 'carolina',]$Flock.size <- 'Medium'
CBC_trait_slopes_flocks[CBC_trait_slopes_flocks$Genus == 'Tyrannus' & CBC_trait_slopes_flocks$Species == 'verticalis',]$Flock.size <- 'Medium'
CBC_trait_slopes_flocks[CBC_trait_slopes_flocks$Genus == 'Bucephala' & CBC_trait_slopes_flocks$Species == 'albeola',]$Flock.size <- 'Small'


CBC_trait_slopes_flocks <- CBC_trait_slopes_flocks %>% mutate(Flock_size = case_when(
  Flocks == 'No' ~ 'Solo',
  Flock.size == 'Small' ~ 'Small flocks',
  Flock.size %in% c('Medium', 'Large') ~ 'Large flocks',
  Flock.size == '' ~ NA,
  Flock.size == '' ~ NA
))

BBS_trait_slopes_flocks$Flock_size <- as.factor(BBS_trait_slopes_flocks$Flock_size)
BBS_trait_slopes_flocks$Flock_size <- relevel(BBS_trait_slopes_flocks$Flock_size, ref='Small flocks')
BBS_trait_slopes_flocks$Flock_size <- relevel(BBS_trait_slopes_flocks$Flock_size, ref='Solo')

CBC_trait_slopes_flocks$Flock_size <- as.factor(CBC_trait_slopes_flocks$Flock_size)
CBC_trait_slopes_flocks$Flock_size <- relevel(CBC_trait_slopes_flocks$Flock_size, ref='Small flocks')
CBC_trait_slopes_flocks$Flock_size <- relevel(CBC_trait_slopes_flocks$Flock_size, ref='Solo')

BBS_BEAU_TSF <- BBS_trait_slopes_flocks
CBC_BEAU_TSF <- CBC_trait_slopes_flocks

BBS_BEAU_TSF_10 <- BBS_BEAU_TSF %>% filter(!is.na(Maximum.flock.size.while.travelling)) %>% mutate(Flock_size = case_when(
  Maximum.flock.size.while.travelling < 11 ~ 'Solo/Small',
  Maximum.flock.size.while.travelling > 10 ~ 'Flocks'
)) 
BBS_BEAU_TSF_5 <- BBS_BEAU_TSF %>% filter(!is.na(Maximum.flock.size.while.travelling)) %>% mutate(Flock_size = case_when(
  Maximum.flock.size.while.travelling < 6 ~ 'Solo/Small',
  Maximum.flock.size.while.travelling > 5 ~ 'Flocks'
)) 

CBC_BEAU_TSF_10 <- CBC_BEAU_TSF %>% filter(!is.na(Maximum.flock.size.while.travelling)) %>% mutate(Flock_size = case_when(
  Maximum.flock.size.while.travelling < 11 ~ 'Solo/Small',
  Maximum.flock.size.while.travelling > 10 ~ 'Flocks'
)) 
CBC_BEAU_TSF_5 <- CBC_BEAU_TSF %>% filter(!is.na(Maximum.flock.size.while.travelling)) %>% mutate(Flock_size = case_when(
  Maximum.flock.size.while.travelling < 6 ~ 'Solo/Small',
  Maximum.flock.size.while.travelling > 5 ~ 'Flocks'
)) 

BBS_BEAU_TSF_5$Flock_size <- as.factor(BBS_BEAU_TSF_5$Flock_size)
BBS_BEAU_TSF_5$Flock_size <- relevel(BBS_BEAU_TSF_5$Flock_size, ref='Solo/Small')

BBS_BEAU_TSF_10$Flock_size <- as.factor(BBS_BEAU_TSF_10$Flock_size)
BBS_BEAU_TSF_10$Flock_size <- relevel(BBS_BEAU_TSF_10$Flock_size, ref='Solo/Small')

CBC_BEAU_TSF_5$Flock_size <- as.factor(CBC_BEAU_TSF_5$Flock_size)
CBC_BEAU_TSF_5$Flock_size <- relevel(CBC_BEAU_TSF_5$Flock_size, ref='Solo/Small')

CBC_BEAU_TSF_10$Flock_size <- as.factor(CBC_BEAU_TSF_10$Flock_size)
CBC_BEAU_TSF_10$Flock_size <- relevel(CBC_BEAU_TSF_10$Flock_size, ref='Solo/Small')

BBS_trait_slopes_flocks <- BBS_trait_slopes_flocks[!is.na(BBS_trait_slopes_flocks$Flock_size),]
CBC_trait_slopes_flocks <- CBC_trait_slopes_flocks[!is.na(CBC_trait_slopes_flocks$Flock_size),]

#### 4.0.0 BBS ----
### 4.1.1 Full model (just flocks) ----

species_list <- unique(CBC_trait_slopes_flocks$CBC_Tree_species)
species_list <- sub(' ', '_',species_list)

BBS_trait_slopes_flocks <- BBS_trait_slopes_flocks %>% filter(AOU %in% CBC_trait_slopes_flocks$AOU)

# Load bird supertree based on Hackett's backbone 
bird_tree_hackett <- read.tree("Data/Hackett.tre")

bird_tree_hackett <- collapse.singles(bird_tree_hackett)

# Get only bird species from the supertree that are also included in collected data
# Extract list of all species from the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)
bird_tree_species <- data.frame(Name = bird_tree_species)

bird_tree_species2 <- bird_tree_species %>% tidyr::separate(Name, c('Genus', 'Species'))
bird_tree_species$Genus <- bird_tree_species2$Genus
bird_tree_species$Species <- bird_tree_species2$Species

# Check the overlap of species names between collected data file and the supertree
# All species should be present. If not, species names may not match with names in the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)

#intersect(bird_tree_species, species_list)
#setdiff(BBS_trait_slopes$BBS_Tree_species,bird_tree_species)

# Prune supertree to the list of taxa included in the data
pruned_birds_stree <- drop.tip(bird_tree_hackett, bird_tree_hackett$tip.label[-match(BBS_trait_slopes_flocks$BBS_Tree_species, bird_tree_hackett$tip.label)])

#is.binary(pruned_birds_stree) # TRUE
#is.ultrametric(pruned_birds_stree) # TRUE

treedat <- data.frame(BBS_Tree_species=pruned_birds_stree$tip.label)
treedat <- join(treedat, BBS_trait_slopes_flocks[c('BBS_Tree_species','Flocks')])
names(treedat)[1] <- 'label'
pruned_birds_stree2 <- full_join(pruned_birds_stree,treedat, by='label')

pruned_birds_stree4 <- pruned_birds_stree2
pruned_birds_stree4@phylo$tip.label <- gsub('_',' ',pruned_birds_stree4@phylo$tip.label)

SE <- BBS_trait_slopes_flocks$BBS_dist_error
SE<-setNames(SE,BBS_trait_slopes_flocks$BBS_Tree_species)
BBS_trait_slopes_flocks$SE <- SE

BBS_trait_slopes_flocks <- BBS_trait_slopes_flocks[match(pruned_birds_stree$tip.label, BBS_trait_slopes_flocks$BBS_Tree_species),]

rownames(BBS_trait_slopes_flocks) <- BBS_trait_slopes_flocks$BBS_Tree_species

BBS_trait_slopes_flocks2 <- BBS_trait_slopes_flocks %>% mutate(
  Habitat_specialism_score = scale(Habitat_specialism_score),
  Diet_specialism_score = scale(Diet_specialism_score),
  GenFLS = scale(GenFLS),
  Trend = scale(abs(Trend)),
  Relative.Abundance = scale(Relative.Abundance),
  Breeding_rangesize = scale(Breeding_rangesize),
  Breeding_overlap_size = scale(Breeding_overlap_size),
  Mig_dist_numeric = scale(Mig_dist_numeric)
)

BBSfit<-pgls.SEy((BBS_dist) ~ Migration
                 + Mig_dist_numeric
                 + Flock_size
                 + Habitat_specialism_score
                 + Diet_specialism_score
                 + GenFLS
                 + Trend
                 + Breeding_rangesize
                 + Breeding_overlap_size
                 + Combined_Migration_Timing
                 ,data=BBS_trait_slopes_flocks2,se=(SE),tree=pruned_birds_stree,method="ML")

car::vif(BBSfit) # Check for vif > 5. None so here.

### 4.1.2 Model averaging (just flocks) ----

model.set <- MuMIn::dredge(BBSfit, m.lim=c(0,nrow(BBS_trait_slopes_flocks2)/10))

top.models <- MuMIn::get.models(model.set, subset = delta <2.0)

averaged.model <- MuMIn::model.avg(top.models) 

#summary(averaged.model)
BBSAVGcoefs <- data.frame(averaged.model$coefficients)
BBSAVGcoefs$Avg <- rownames(BBSAVGcoefs)
BBSAVGcoefs <- melt(BBSAVGcoefs)
levels(BBSAVGcoefs$variable)[1] <-  '(Intercept)'

Conf <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                           level = 0.95,    # lets you adjust the confidence interval
                           full = T))    
Conf$variable <- rownames(Conf)
Conf$Avg <- 'full' 

Conf2 <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                            level = 0.95,    # lets you adjust the confidence interval
                            full = F))    
Conf2$variable <- rownames(Conf2)
Conf2$Avg <- 'subset'

BBSAVGcoefs <- join(BBSAVGcoefs, Conf)

BBSAVGcoefs[BBSAVGcoefs$Avg == 'subset',]$X2.5.. <- Conf2$X2.5..
BBSAVGcoefs[BBSAVGcoefs$Avg == 'subset',]$X97.5.. <- Conf2$X97.5..

BBSAVGcoefs <- BBSAVGcoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  variable == '(Intercept)' ~ 'Intercept',
  variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  variable == 'Mig_dist_numeric' ~ 'Migratory distance',
  variable == 'Flock_sizeSmall.flocks' ~ 'Small flocks\n(ref: Solo)',
  variable == 'Flock_sizeLarge.flocks' ~ 'Large flocks\n(ref: Solo)',
  variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  variable == 'GenFLS' ~ 'Generation length',
  variable == 'Trend' ~ 'Absolute population trend',
  variable == 'Breeding_rangesize' ~ 'Total range size',
  variable == 'Breeding_overlap_size' ~ 'Overlap range size',
  variable == 'Combined_Migration_Timingnight' ~ 'Migration timing\n(ref: day)'
))

BBSAVGplot <- ggplot(BBSAVGcoefs, aes(x=reorder(variable, value), y=value, colour=Significance, shape = Avg)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.25)) +
  labs(y='Coefficient', x='', shape='Average', title='Breeding model A', subtitle='122 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="right")

BBSAVGplot

BBSAVGcoefs_A <- BBSAVGcoefs
BBSAVGcoefs_A$n_mod <- length(top.models)

### 4.1.3 Full model (BEAU_TSF_5) ----

species_list <- unique(CBC_BEAU_TSF_5$CBC_Tree_species)
species_list <- sub(' ', '_',species_list)

BBS_BEAU_TSF_5 <- BBS_BEAU_TSF_5 %>% filter(AOU %in% CBC_BEAU_TSF_5$AOU)

# Load bird supertree based on Hackett's backbone 
bird_tree_hackett <- read.tree("Data/Hackett.tre")

bird_tree_hackett <- collapse.singles(bird_tree_hackett)

# Get only bird species from the supertree that are also included in collected data
# Extract list of all species from the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)
bird_tree_species <- data.frame(Name = bird_tree_species)

bird_tree_species2 <- bird_tree_species %>% tidyr::separate(Name, c('Genus', 'Species'))
bird_tree_species$Genus <- bird_tree_species2$Genus
bird_tree_species$Species <- bird_tree_species2$Species

# Check the overlap of species names between collected data file and the supertree
# All species should be present. If not, species names may not match with names in the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)

#intersect(bird_tree_species, species_list)
#setdiff(BBS_trait_slopes$BBS_Tree_species,bird_tree_species)

# Prune supertree to the list of taxa included in the data
pruned_birds_stree <- drop.tip(bird_tree_hackett, bird_tree_hackett$tip.label[-match(BBS_BEAU_TSF_5$BBS_Tree_species, bird_tree_hackett$tip.label)])

#is.binary(pruned_birds_stree) # TRUE
#is.ultrametric(pruned_birds_stree) # TRUE

treedat <- data.frame(BBS_Tree_species=pruned_birds_stree$tip.label)
treedat <- join(treedat,BBS_BEAU_TSF_5[c('BBS_Tree_species','Flocks')])
names(treedat)[1] <- 'label'
pruned_birds_stree2 <- full_join(pruned_birds_stree,treedat, by='label')

pruned_birds_stree4 <- pruned_birds_stree2
pruned_birds_stree4@phylo$tip.label <- gsub('_',' ',pruned_birds_stree4@phylo$tip.label)

SE <- BBS_BEAU_TSF_5$BBS_dist_error
SE<-setNames(SE,BBS_BEAU_TSF_5$BBS_Tree_species)
BBS_BEAU_TSF_5$SE <- SE

BBS_BEAU_TSF_5 <- BBS_BEAU_TSF_5[match(pruned_birds_stree$tip.label, BBS_BEAU_TSF_5$BBS_Tree_species),]

rownames(BBS_BEAU_TSF_5) <- BBS_BEAU_TSF_5$BBS_Tree_species

BBS_BEAU_TSF_5_2 <- BBS_BEAU_TSF_5 %>% mutate(
  Habitat_specialism_score = scale(Habitat_specialism_score),
  Diet_specialism_score = scale(Diet_specialism_score),
  GenFLS = scale(GenFLS),
  Trend = scale(abs(Trend)),
  Relative.Abundance = scale(Relative.Abundance),
  Breeding_rangesize = scale(Breeding_rangesize),
  Breeding_overlap_size = scale(Breeding_overlap_size),
  Mig_dist_numeric = scale(Mig_dist_numeric)
)

BBSfit<-pgls.SEy((BBS_dist) ~ Migration
                 + Mig_dist_numeric
                 + Flock_size
                 + Habitat_specialism_score
                 + Diet_specialism_score
                 + GenFLS
                 + Trend
                 #+ Breeding_rangesize
                 + Breeding_overlap_size
                 + Combined_Migration_Timing
                 ,data=BBS_BEAU_TSF_5_2,se=(SE),tree=pruned_birds_stree,method="ML")

car::vif(BBSfit) # Check for vif > 5. None so here.

### 4.1.4 Model averaging (BEAU_TSF_5) ----

model.set <- MuMIn::dredge(BBSfit, m.lim=c(0,nrow(BBS_BEAU_TSF_5_2)/10))

top.models <- MuMIn::get.models(model.set, subset = delta <2.0)

averaged.model <- MuMIn::model.avg(top.models) 

#summary(averaged.model)
BBSAVGcoefs <- data.frame(averaged.model$coefficients)
BBSAVGcoefs$Avg <- rownames(BBSAVGcoefs)
BBSAVGcoefs <- melt(BBSAVGcoefs)
levels(BBSAVGcoefs$variable)[1] <-  '(Intercept)'

Conf <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                           level = 0.95,    # lets you adjust the confidence interval
                           full = T))    
Conf$variable <- rownames(Conf)
Conf$Avg <- 'full' 

Conf2 <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                            level = 0.95,    # lets you adjust the confidence interval
                            full = F))    
Conf2$variable <- rownames(Conf2)
Conf2$Avg <- 'subset'

BBSAVGcoefs <- join(BBSAVGcoefs, Conf)

BBSAVGcoefs[BBSAVGcoefs$Avg == 'subset',]$X2.5.. <- Conf2$X2.5..
BBSAVGcoefs[BBSAVGcoefs$Avg == 'subset',]$X97.5.. <- Conf2$X97.5..

BBSAVGcoefs <- BBSAVGcoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  variable == '(Intercept)' ~ 'Intercept',
  variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  variable == 'Mig_dist_numeric' ~ 'Migratory distance',
  variable == 'Flock_sizeLarge' ~ 'Large flocks\n(ref: Solo)',
  variable == 'Flock_sizeSmall' ~ 'Small flocks\n(ref: Solo)',
  variable == 'Flock_sizeFlocks' ~ 'Flocks\n(ref: Solo/Small)',
  variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  variable == 'GenFLS' ~ 'Generation length',
  variable == 'Trend' ~ 'Absolute population trend',
  variable == 'Breeding_rangesize' ~ 'Total range size',
  variable == 'Breeding_overlap_size' ~ 'Overlap range size',
  variable == 'Combined_Migration_Timingnight' ~ 'Migration timing\n(ref: day)'
))

BBSAVGplot <- ggplot(BBSAVGcoefs, aes(x=reorder(Variable, value), y=value, colour=Significance, shape = Avg)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(posiiton=position_dodge(), position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.25)) +
  labs(y='Coefficient', x='', shape='Average', title='Breeding model A', subtitle='122 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="right")

BBSAVGplot

BBSAVGcoefs_B <- BBSAVGcoefs
BBSAVGcoefs_B$n_mod <- length(top.models)



### 4.1.3 Full model (BEAU_TSF_10) ----

species_list <- unique(CBC_BEAU_TSF_10$CBC_Tree_species)
species_list <- sub(' ', '_',species_list)

BBS_BEAU_TSF_10 <- BBS_BEAU_TSF_10 %>% filter(AOU %in% CBC_BEAU_TSF_10$AOU)

# Load bird supertree based on Hackett's backbone 
bird_tree_hackett <- read.tree("Data/Hackett.tre")

bird_tree_hackett <- collapse.singles(bird_tree_hackett)

# Get only bird species from the supertree that are also included in collected data
# Extract list of all species from the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)
bird_tree_species <- data.frame(Name = bird_tree_species)

bird_tree_species2 <- bird_tree_species %>% tidyr::separate(Name, c('Genus', 'Species'))
bird_tree_species$Genus <- bird_tree_species2$Genus
bird_tree_species$Species <- bird_tree_species2$Species

# Check the overlap of species names between collected data file and the supertree
# All species should be present. If not, species names may not match with names in the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)

#intersect(bird_tree_species, species_list)
#setdiff(BBS_trait_slopes$BBS_Tree_species,bird_tree_species)

# Prune supertree to the list of taxa included in the data
pruned_birds_stree <- drop.tip(bird_tree_hackett, bird_tree_hackett$tip.label[-match(BBS_BEAU_TSF_10$BBS_Tree_species, bird_tree_hackett$tip.label)])

#is.binary(pruned_birds_stree) # TRUE
#is.ultrametric(pruned_birds_stree) # TRUE

treedat <- data.frame(BBS_Tree_species=pruned_birds_stree$tip.label)
treedat <- join(treedat, BBS_BEAU_TSF_10[c('BBS_Tree_species','Flocks')])
names(treedat)[1] <- 'label'
pruned_birds_stree2 <- full_join(pruned_birds_stree,treedat, by='label')

pruned_birds_stree4 <- pruned_birds_stree2
pruned_birds_stree4@phylo$tip.label <- gsub('_',' ',pruned_birds_stree4@phylo$tip.label)

SE <- BBS_BEAU_TSF_10$BBS_dist_error
SE<-setNames(SE,BBS_BEAU_TSF_10$BBS_Tree_species)
BBS_BEAU_TSF_10$SE <- SE

BBS_BEAU_TSF_10 <- BBS_BEAU_TSF_10[match(pruned_birds_stree$tip.label, BBS_BEAU_TSF_10$BBS_Tree_species),]

rownames(BBS_BEAU_TSF_10) <- BBS_BEAU_TSF_10$BBS_Tree_species

BBS_BEAU_TSF_10_2 <- BBS_BEAU_TSF_10 %>% mutate(
  Habitat_specialism_score = scale(Habitat_specialism_score),
  Diet_specialism_score = scale(Diet_specialism_score),
  GenFLS = scale(GenFLS),
  Trend = scale(abs(Trend)),
  Relative.Abundance = scale(Relative.Abundance),
  Breeding_rangesize = scale(Breeding_rangesize),
  Breeding_overlap_size = scale(Breeding_overlap_size),
  Mig_dist_numeric = scale(Mig_dist_numeric)
)

BBSfit<-pgls.SEy((BBS_dist) ~ Migration
                 + Mig_dist_numeric
                 + Flock_size
                 + Habitat_specialism_score
                 + Diet_specialism_score
                 + GenFLS
                 + Trend
                 #+ Breeding_rangesize
                 + Breeding_overlap_size
                 + Combined_Migration_Timing
                 ,data=BBS_BEAU_TSF_10_2,se=(SE),tree=pruned_birds_stree,method="ML")

car::vif(BBSfit) # Check for vif > 5. None so here.

### 4.1.4 Model averaging (BEAU_TSF_10) ----

model.set <- MuMIn::dredge(BBSfit, m.lim=c(0,nrow(BBS_BEAU_TSF_10_2)/10))

top.models <- MuMIn::get.models(model.set, subset = delta <2.0)

averaged.model <- MuMIn::model.avg(top.models) 

#summary(averaged.model)
BBSAVGcoefs <- data.frame(averaged.model$coefficients)
BBSAVGcoefs$Avg <- rownames(BBSAVGcoefs)
BBSAVGcoefs <- melt(BBSAVGcoefs)
levels(BBSAVGcoefs$variable)[1] <-  '(Intercept)'

Conf <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                           level = 0.95,    # lets you adjust the confidence interval
                           full = T))    
Conf$variable <- rownames(Conf)
Conf$Avg <- 'full' 

Conf2 <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                            level = 0.95,    # lets you adjust the confidence interval
                            full = F))    
Conf2$variable <- rownames(Conf2)
Conf2$Avg <- 'subset'

BBSAVGcoefs <- join(BBSAVGcoefs, Conf)

BBSAVGcoefs[BBSAVGcoefs$Avg == 'subset',]$X2.5.. <- Conf2$X2.5..
BBSAVGcoefs[BBSAVGcoefs$Avg == 'subset',]$X97.5.. <- Conf2$X97.5..

BBSAVGcoefs <- BBSAVGcoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  variable == '(Intercept)' ~ 'Intercept',
  variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  variable == 'Mig_dist_numeric' ~ 'Migratory distance',
  variable == 'Flock_sizeLarge' ~ 'Large flocks\n(ref: Solo)',
  variable == 'Flock_sizeSmall' ~ 'Small flocks\n(ref: Solo)',
  variable == 'Flock_sizeFlocks' ~ 'Flocks\n(ref: Solo/Small)',
  variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  variable == 'GenFLS' ~ 'Generation length',
  variable == 'Trend' ~ 'Absolute population trend',
  variable == 'Breeding_rangesize' ~ 'Total range size',
  variable == 'Breeding_overlap_size' ~ 'Overlap range size',
  variable == 'Combined_Migration_Timingnight' ~ 'Migration timing\n(ref: day)'
))

BBSAVGplot <- ggplot(BBSAVGcoefs, aes(x=reorder(Variable, value), y=value, colour=Significance, shape = Avg)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(posiiton=position_dodge(), position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.25)) +
  labs(y='Coefficient', x='', shape='Average', title='Breeding model A', subtitle='122 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="right")

BBSAVGplot

BBSAVGcoefs_C <- BBSAVGcoefs
BBSAVGcoefs_C$n_mod <- length(top.models)


#### 4.0.0 CBC ----
### 4.1.1 Full model (just flocks) ----

species_list <- unique(CBC_trait_slopes_flocks$CBC_Tree_species)
species_list <- sub(' ', '_',species_list)

CBC_trait_slopes_flocks <- CBC_trait_slopes_flocks %>% filter(AOU %in% CBC_trait_slopes_flocks$AOU)

# Load bird supertree based on Hackett's backbone 
bird_tree_hackett <- read.tree("Data/Hackett.tre")

bird_tree_hackett <- collapse.singles(bird_tree_hackett)

# Get only bird species from the supertree that are also included in collected data
# Extract list of all species from the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)
bird_tree_species <- data.frame(Name = bird_tree_species)

bird_tree_species2 <- bird_tree_species %>% tidyr::separate(Name, c('Genus', 'Species'))
bird_tree_species$Genus <- bird_tree_species2$Genus
bird_tree_species$Species <- bird_tree_species2$Species

# Check the overlap of species names between collected data file and the supertree
# All species should be present. If not, species names may not match with names in the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)

#intersect(bird_tree_species, species_list)
#setdiff(CBC_trait_slopes$CBC_Tree_species,bird_tree_species)

# Prune supertree to the list of taxa included in the data
pruned_birds_stree <- drop.tip(bird_tree_hackett, bird_tree_hackett$tip.label[-match(CBC_trait_slopes_flocks$CBC_Tree_species, bird_tree_hackett$tip.label)])

#is.binary(pruned_birds_stree) # TRUE
#is.ultrametric(pruned_birds_stree) # TRUE

treedat <- data.frame(CBC_Tree_species=pruned_birds_stree$tip.label)
treedat <- join(treedat, CBC_trait_slopes_flocks[c('CBC_Tree_species','Flocks')])
names(treedat)[1] <- 'label'
pruned_birds_stree2 <- full_join(pruned_birds_stree,treedat, by='label')

pruned_birds_stree4 <- pruned_birds_stree2
pruned_birds_stree4@phylo$tip.label <- gsub('_',' ',pruned_birds_stree4@phylo$tip.label)

SE <- CBC_trait_slopes_flocks$CBC_dist_error
SE<-setNames(SE,CBC_trait_slopes_flocks$CBC_Tree_species)
CBC_trait_slopes_flocks$SE <- SE

CBC_trait_slopes_flocks <- CBC_trait_slopes_flocks[match(pruned_birds_stree$tip.label, CBC_trait_slopes_flocks$CBC_Tree_species),]

rownames(CBC_trait_slopes_flocks) <- CBC_trait_slopes_flocks$CBC_Tree_species

CBC_trait_slopes_flocks2 <- CBC_trait_slopes_flocks %>% mutate(
  Habitat_specialism_score = scale(Habitat_specialism_score),
  Diet_specialism_score = scale(Diet_specialism_score),
  GenFLS = scale(GenFLS),
  Trend = scale(abs(Trend)),
  Relative.Abundance = scale(Relative.Abundance),
  Breeding_rangesize = scale(Breeding_rangesize),
  Breeding_overlap_size = scale(Breeding_overlap_size),
  Mig_dist_numeric = scale(Mig_dist_numeric)
)

CBCfit<-pgls.SEy((CBC_dist) ~ Migration
                 + Mig_dist_numeric
                 + Flock_size
                 + Habitat_specialism_score
                 + Diet_specialism_score
                 + GenFLS
                 + Trend
                 + Breeding_rangesize
                 + Breeding_overlap_size
                 + Combined_Migration_Timing
                 ,data=CBC_trait_slopes_flocks2,se=(SE),tree=pruned_birds_stree,method="ML")

car::vif(CBCfit) # Check for vif > 5. None so here.

### 4.1.2 Model averaging (just flocks) ----

model.set <- MuMIn::dredge(CBCfit, m.lim=c(0,nrow(CBC_trait_slopes_flocks2)/10))

top.models <- MuMIn::get.models(model.set, subset = delta <2.0)

averaged.model <- MuMIn::model.avg(top.models) 

#summary(averaged.model)
CBCAVGcoefs <- data.frame(averaged.model$coefficients)
CBCAVGcoefs$Avg <- rownames(CBCAVGcoefs)
CBCAVGcoefs <- melt(CBCAVGcoefs)
levels(CBCAVGcoefs$variable)[1] <-  '(Intercept)'

Conf <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                           level = 0.95,    # lets you adjust the confidence interval
                           full = T))    
Conf$variable <- rownames(Conf)
Conf$Avg <- 'full' 

Conf2 <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                            level = 0.95,    # lets you adjust the confidence interval
                            full = F))    
Conf2$variable <- rownames(Conf2)
Conf2$Avg <- 'subset'

Conf$variable[Conf$variable == "Flock_sizeSmall flocks"] <- 'Flock_sizeSmall.flocks'
Conf$variable[Conf$variable == "Flock_sizeLarge flocks"] <- 'Flock_sizeLarge.flocks'

CBCAVGcoefs <- join(CBCAVGcoefs, Conf)

CBCAVGcoefs[CBCAVGcoefs$Avg == 'subset',]$X2.5.. <- Conf2$X2.5..
CBCAVGcoefs[CBCAVGcoefs$Avg == 'subset',]$X97.5.. <- Conf2$X97.5..
CBCAVGcoefs <- CBCAVGcoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  variable == '(Intercept)' ~ 'Intercept',
  variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  variable == 'Mig_dist_numeric' ~ 'Migratory distance',
  variable == 'Flock_sizeLarge.flocks' ~ 'Large flocks\n(ref: Solo)',
  variable == 'Flock_sizeSmall.flocks' ~ 'Small flocks\n(ref: Solo)',
  variable == 'Flock_sizeFlocks' ~ 'Flocks\n(ref: Solo/Small)',
  variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  variable == 'GenFLS' ~ 'Generation length',
  variable == 'Trend' ~ 'Absolute population trend',
  variable == 'Breeding_rangesize' ~ 'Total range size',
  variable == 'Breeding_overlap_size' ~ 'Overlap range size',
  variable == 'Combined_Migration_Timingnight' ~ 'Migration timing\n(ref: day)'
))

CBCAVGplot <- ggplot(CBCAVGcoefs, aes(x=reorder(variable, value), y=value, colour=Significance, shape = Avg)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.25)) +
  labs(y='Coefficient', x='', shape='Average', title='Breeding model A', subtitle='122 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="right")

CBCAVGplot

CBCAVGcoefs_A <- CBCAVGcoefs
CBCAVGcoefs_A$n_mod <- length(top.models)




### 4.1.3 Full model (BEAU_TSF_5) ----

species_list <- unique(CBC_BEAU_TSF_5$CBC_Tree_species)
species_list <- sub(' ', '_',species_list)

CBC_BEAU_TSF_5 <- CBC_BEAU_TSF_5 %>% filter(AOU %in% CBC_BEAU_TSF_5$AOU)

# Load bird supertree based on Hackett's backbone 
bird_tree_hackett <- read.tree("Data/Hackett.tre")

bird_tree_hackett <- collapse.singles(bird_tree_hackett)

# Get only bird species from the supertree that are also included in collected data
# Extract list of all species from the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)
bird_tree_species <- data.frame(Name = bird_tree_species)

bird_tree_species2 <- bird_tree_species %>% tidyr::separate(Name, c('Genus', 'Species'))
bird_tree_species$Genus <- bird_tree_species2$Genus
bird_tree_species$Species <- bird_tree_species2$Species

# Check the overlap of species names between collected data file and the supertree
# All species should be present. If not, species names may not match with names in the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)

#intersect(bird_tree_species, species_list)
#setdiff(CBC_trait_slopes$CBC_Tree_species,bird_tree_species)

# Prune supertree to the list of taxa included in the data
pruned_birds_stree <- drop.tip(bird_tree_hackett, bird_tree_hackett$tip.label[-match(CBC_BEAU_TSF_5$CBC_Tree_species, bird_tree_hackett$tip.label)])

#is.binary(pruned_birds_stree) # TRUE
#is.ultrametric(pruned_birds_stree) # TRUE

treedat <- data.frame(CBC_Tree_species=pruned_birds_stree$tip.label)
treedat <- join(treedat,CBC_BEAU_TSF_5[c('CBC_Tree_species','Flocks')])
names(treedat)[1] <- 'label'
pruned_birds_stree2 <- full_join(pruned_birds_stree,treedat, by='label')

pruned_birds_stree4 <- pruned_birds_stree2
pruned_birds_stree4@phylo$tip.label <- gsub('_',' ',pruned_birds_stree4@phylo$tip.label)

SE <- CBC_BEAU_TSF_5$CBC_dist_error
SE<-setNames(SE,CBC_BEAU_TSF_5$CBC_Tree_species)
CBC_BEAU_TSF_5$SE <- SE

CBC_BEAU_TSF_5 <- CBC_BEAU_TSF_5[match(pruned_birds_stree$tip.label, CBC_BEAU_TSF_5$CBC_Tree_species),]

rownames(CBC_BEAU_TSF_5) <- CBC_BEAU_TSF_5$CBC_Tree_species

CBC_BEAU_TSF_5_2 <- CBC_BEAU_TSF_5 %>% mutate(
  Habitat_specialism_score = scale(Habitat_specialism_score),
  Diet_specialism_score = scale(Diet_specialism_score),
  GenFLS = scale(GenFLS),
  Trend = scale(abs(Trend)),
  Relative.Abundance = scale(Relative.Abundance),
  Breeding_rangesize = scale(Breeding_rangesize),
  Breeding_overlap_size = scale(Breeding_overlap_size),
  Mig_dist_numeric = scale(Mig_dist_numeric)
)

CBCfit<-pgls.SEy((CBC_dist) ~ Migration
                 + Mig_dist_numeric
                 + Flock_size
                 + Habitat_specialism_score
                 + Diet_specialism_score
                 + GenFLS
                 + Trend
                 + Breeding_rangesize
                 + Breeding_overlap_size
                 + Combined_Migration_Timing
                 ,data=CBC_BEAU_TSF_5_2,se=(SE),tree=pruned_birds_stree,method="ML")

car::vif(CBCfit) # Check for vif > 5. None so here.

### 4.1.4 Model averaging (BEAU_TSF_5) ----

model.set <- MuMIn::dredge(CBCfit, m.lim=c(0,nrow(CBC_BEAU_TSF_5_2)/10))

top.models <- MuMIn::get.models(model.set, subset = delta <2.0)

averaged.model <- MuMIn::model.avg(top.models) 

#summary(averaged.model)
CBCAVGcoefs <- data.frame(averaged.model$coefficients)
CBCAVGcoefs$Avg <- rownames(CBCAVGcoefs)
CBCAVGcoefs <- melt(CBCAVGcoefs)
levels(CBCAVGcoefs$variable)[1] <-  '(Intercept)'

Conf <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                           level = 0.95,    # lets you adjust the confidence interval
                           full = T))    
Conf$variable <- rownames(Conf)
Conf$Avg <- 'full' 

Conf2 <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                            level = 0.95,    # lets you adjust the confidence interval
                            full = F))    
Conf2$variable <- rownames(Conf2)
Conf2$Avg <- 'subset'

CBCAVGcoefs <- join(CBCAVGcoefs, Conf)

CBCAVGcoefs[CBCAVGcoefs$Avg == 'subset',]$X2.5.. <- Conf2$X2.5..
CBCAVGcoefs[CBCAVGcoefs$Avg == 'subset',]$X97.5.. <- Conf2$X97.5..

CBCAVGcoefs <- CBCAVGcoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  variable == '(Intercept)' ~ 'Intercept',
  variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  variable == 'Mig_dist_numeric' ~ 'Migratory distance',
  variable == 'Flock_sizeFlocks' ~ 'Flocking\n(ref: Solo/Small)',
  variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  variable == 'GenFLS' ~ 'Generation length',
  variable == 'Trend' ~ 'Absolute population trend',
  variable == 'Breeding_rangesize' ~ 'Total range size',
  variable == 'Breeding_overlap_size' ~ 'Overlap range size',
  variable == 'Combined_Migration_Timingnight' ~ 'Migration timing\n(ref: day)'
))

CBCAVGplot <- ggplot(CBCAVGcoefs, aes(x=reorder(Variable, value), y=value, colour=Significance, shape = Avg)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(posiiton=position_dodge(), position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.25)) +
  labs(y='Coefficient', x='', shape='Average', title='Breeding model A', subtitle='122 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="right")

CBCAVGplot

CBCAVGcoefs_B <- CBCAVGcoefs
CBCAVGcoefs_B$n_mod <- length(top.models)

### 4.1.3 Full model (BEAU_TSF_10) ----

species_list <- unique(CBC_BEAU_TSF_10$CBC_Tree_species)
species_list <- sub(' ', '_',species_list)

CBC_BEAU_TSF_10 <- CBC_BEAU_TSF_10 %>% filter(AOU %in% CBC_BEAU_TSF_10$AOU)

# Load bird supertree based on Hackett's backbone 
bird_tree_hackett <- read.tree("Data/Hackett.tre")

bird_tree_hackett <- collapse.singles(bird_tree_hackett)

# Get only bird species from the supertree that are also included in collected data
# Extract list of all species from the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)
bird_tree_species <- data.frame(Name = bird_tree_species)

bird_tree_species2 <- bird_tree_species %>% tidyr::separate(Name, c('Genus', 'Species'))
bird_tree_species$Genus <- bird_tree_species2$Genus
bird_tree_species$Species <- bird_tree_species2$Species

# Check the overlap of species names between collected data file and the supertree
# All species should be present. If not, species names may not match with names in the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)

#intersect(bird_tree_species, species_list)
#setdiff(CBC_trait_slopes$CBC_Tree_species,bird_tree_species)

# Prune supertree to the list of taxa included in the data
pruned_birds_stree <- drop.tip(bird_tree_hackett, bird_tree_hackett$tip.label[-match(CBC_BEAU_TSF_10$CBC_Tree_species, bird_tree_hackett$tip.label)])

#is.binary(pruned_birds_stree) # TRUE
#is.ultrametric(pruned_birds_stree) # TRUE

treedat <- data.frame(CBC_Tree_species=pruned_birds_stree$tip.label)
treedat <- join(treedat, CBC_BEAU_TSF_10[c('CBC_Tree_species','Flocks')])
names(treedat)[1] <- 'label'
pruned_birds_stree2 <- full_join(pruned_birds_stree,treedat, by='label')

pruned_birds_stree4 <- pruned_birds_stree2
pruned_birds_stree4@phylo$tip.label <- gsub('_',' ',pruned_birds_stree4@phylo$tip.label)

SE <- CBC_BEAU_TSF_10$CBC_dist_error
SE<-setNames(SE,CBC_BEAU_TSF_10$CBC_Tree_species)
CBC_BEAU_TSF_10$SE <- SE

CBC_BEAU_TSF_10 <- CBC_BEAU_TSF_10[match(pruned_birds_stree$tip.label, CBC_BEAU_TSF_10$CBC_Tree_species),]

rownames(CBC_BEAU_TSF_10) <- CBC_BEAU_TSF_10$CBC_Tree_species

CBC_BEAU_TSF_10_2 <- CBC_BEAU_TSF_10 %>% mutate(
  Habitat_specialism_score = scale(Habitat_specialism_score),
  Diet_specialism_score = scale(Diet_specialism_score),
  GenFLS = scale(GenFLS),
  Trend = scale(abs(Trend)),
  Relative.Abundance = scale(Relative.Abundance),
  Breeding_rangesize = scale(Breeding_rangesize),
  Breeding_overlap_size = scale(Breeding_overlap_size),
  Mig_dist_numeric = scale(Mig_dist_numeric)
)

CBCfit<-pgls.SEy((CBC_dist) ~ Migration
                 + Mig_dist_numeric
                 + Flock_size
                 + Habitat_specialism_score
                 + Diet_specialism_score
                 + GenFLS
                 + Trend
                 + Breeding_rangesize
                 + Breeding_overlap_size
                 + Combined_Migration_Timing
                 ,data=CBC_BEAU_TSF_10_2,se=(SE),tree=pruned_birds_stree,method="ML")

car::vif(CBCfit) # Check for vif > 5. None so here.

### 4.1.4 Model averaging (BEAU_TSF_10) ----

model.set <- MuMIn::dredge(CBCfit, m.lim=c(0,nrow(CBC_BEAU_TSF_10_2)/10))

top.models <- MuMIn::get.models(model.set, subset = delta <2.0)

averaged.model <- MuMIn::model.avg(top.models) 

#summary(averaged.model)
CBCAVGcoefs <- data.frame(averaged.model$coefficients)
CBCAVGcoefs$Avg <- rownames(CBCAVGcoefs)
CBCAVGcoefs <- melt(CBCAVGcoefs)
levels(CBCAVGcoefs$variable)[1] <-  '(Intercept)'

Conf <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                           level = 0.95,    # lets you adjust the confidence interval
                           full = T))    
Conf$variable <- rownames(Conf)
Conf$Avg <- 'full' 

Conf2 <- data.frame(confint(averaged.model,  # Confidence intervals for averaged parameter estimates 
                            level = 0.95,    # lets you adjust the confidence interval
                            full = F))    
Conf2$variable <- rownames(Conf2)
Conf2$Avg <- 'subset'

CBCAVGcoefs <- join(CBCAVGcoefs, Conf)

CBCAVGcoefs[CBCAVGcoefs$Avg == 'subset',]$X2.5.. <- Conf2$X2.5..
CBCAVGcoefs[CBCAVGcoefs$Avg == 'subset',]$X97.5.. <- Conf2$X97.5..

CBCAVGcoefs <- CBCAVGcoefs %>% mutate(Significance = case_when(
  X2.5.. > 0 & X97.5.. > 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. < 0 ~ 'Significant',
  X2.5.. < 0 & X97.5.. > 0 ~ 'Non-significant',
  X2.5.. > 0 & X97.5.. < 0 ~ 'Non-significant'
),
Variable = case_when(
  variable == '(Intercept)' ~ 'Intercept',
  variable == 'MigrationPartial' ~ 'Partial migrant\n(ref: complete)',
  variable == 'Mig_dist_numeric' ~ 'Migratory distance',
  variable == 'Flock_sizeFlocks' ~ 'Flocking\n(ref: Solo/Small)',
  variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  variable == 'GenFLS' ~ 'Generation length',
  variable == 'Trend' ~ 'Absolute population trend',
  variable == 'Breeding_rangesize' ~ 'Total range size',
  variable == 'Breeding_overlap_size' ~ 'Overlap range size',
  variable == 'Combined_Migration_Timingnight' ~ 'Migration timing\n(ref: day)'
))

CBCAVGplot <- ggplot(CBCAVGcoefs, aes(x=reorder(Variable, value), y=value, colour=Significance, shape = Avg)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(posiiton=position_dodge(), position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.25)) +
  labs(y='Coefficient', x='', shape='Average', title='Breeding model A', subtitle='122 species asssessed') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="right")

CBCAVGplot

CBCAVGcoefs_C <- CBCAVGcoefs
CBCAVGcoefs_C$n_mod <- length(top.models)

#model.set2 <- data.frame(model.set)
#model.set2 <- model.set2 %>% 
#  mutate_if(is.numeric, round, digits=2)
#model.set2 <- model.set2 %>% filter(delta < 2)
#write.csv(model.set2, 'table19.csv', row.names=F)
#
#tabset <- CBCAVGcoefs_C%>% 
#  mutate_if(is.numeric, round, digits=2)
#tabset <- tabset %>% filter(Avg =='full')
#names(tabset) <- c("Avg","null","Estimate","2.5%","97.5%","Significance","Variable","n_mod")
#tabset <- tabset %>% dplyr::select(Variable,Estimate,`2.5%`,`97.5%`,Significance,n_mod)
#write.csv(tabset, 'table20.csv', row.names=F)

#### plots ----

BBSAVGcoefs_A$Model <- 'Model A'
BBSAVGcoefs_B$Model <- 'Model A\n(Flocks > 5 inds.)'
BBSAVGcoefs_C$Model <- 'Model B\n(Flocks > 10 inds.)'
CBCAVGcoefs_A$Model <- 'Model A'
CBCAVGcoefs_B$Model <- 'Model A\n(Flocks > 5 inds.)'
CBCAVGcoefs_C$Model <- 'Model B\n(Flocks > 10 inds.)'

BBSAVGcoefs_A$Season <- 'Breeding'
BBSAVGcoefs_B$Season <- 'Breeding'
BBSAVGcoefs_C$Season <- 'Breeding'
CBCAVGcoefs_A$Season <- 'Non-breeding'
CBCAVGcoefs_B$Season <- 'Non-breeding'
CBCAVGcoefs_C$Season <- 'Non-breeding'


Sens_AVGcoefs <- rbind(BBSAVGcoefs_A,
                  BBSAVGcoefs_B,
                  BBSAVGcoefs_C,
                  CBCAVGcoefs_A,
                  CBCAVGcoefs_B,
                  CBCAVGcoefs_C)
Sens_AVGcoefs <- Sens_AVGcoefs %>% filter(Avg == 'full')
#Sens_AVGcoefs[is.na(Sens_AVGcoefs$Variable),]$Variable <- 'MigrationPartial'

Sens_AVGcoefsA <- Sens_AVGcoefs %>% filter(Model %in% c('Model A'))
Sens_AVGcoefsB <- Sens_AVGcoefs %>% filter(Model %in% c('Model A\n(Flocks > 5 inds.)', 'Model B\n(Flocks > 10 inds.)'))


results_plotA <- ggplot(Sens_AVGcoefsA, aes(x=reorder(Variable,value), y=value, colour=Significance)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.4) +
  geom_point(size=1.5,position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.5)) +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  facet_grid(~ Season) +
  labs(y='Coefficient') +
  theme(axis.title.x = element_text(size=13), axis.text.y = element_text(size=10), legend.position="bottom", axis.title.y = element_blank(), strip.background = element_blank(), legend.title = element_blank(), legend.text = element_text(size=10), plot.background = element_rect(fill='white'))

results_plotA

results_plotB <- ggplot(Sens_AVGcoefsB, aes(x=reorder(Variable,value), y=value, colour=Significance, shape=Model)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.4) +
  geom_point(size=1.5,position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.5)) +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  facet_grid(~ Season) +
  labs(y='Coefficient') +
  theme(axis.title.x = element_text(size=13), axis.text.y = element_text(size=10), legend.position="bottom", axis.title.y = element_blank(), strip.background = element_blank(), legend.title = element_blank(), legend.text = element_text(size=10), plot.background = element_rect(fill='white'))

results_plotB

ggsave(results_plotA, file='Sens_tetsing_figure_1.pdf', width=190, height=180, units='mm', dpi=600)
ggsave(results_plotB, file='Sens_tetsing_figure_2.pdf', width=190, height=180, units='mm', dpi=600)


CBC_trait_slopes_flocks2 %>% filter(Flock_size == 'Large flocks') %>% dplyr::select('Genus', 'Species', 'CBC_dist')
