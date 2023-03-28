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
               readr)
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
CBC <- read_csv('cbc_species_by_stratum_by_year_abundance_indices.csv') # This is CBC stratum level population estimates by species and year
BBS <- read_csv('all 2019 BBS indices.csv') # This is BBS stratum level population estimates by species and year
beauchamp <- read_csv('MS_data.csv') # This is data from Beauchamp (2011) combined with data taken from Birds of the World Online and oldbird.org

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

Both_final3 <- read_csv('Strata_metrics_all.csv')
AOUS <- read.csv('AOU_list.csv') # This is just AOU species codes we are interested in. This list was created to save computational time, and is based upon known missing traits which means the excluuded species cannot be used in later analysis.

Both_final3 <- Both_final3 %>% filter(AOU %in% AOUS$AOUS)

CBC <- Both_final3 %>% filter(!(is.na(CBC_SD))) %>% dplyr::select(1:5, 9:15, 17)

# Assign errors to numeric values
errors(CBC$CBC_median) <- CBC$CBC_SD
errors(CBC$xcoord) <- 100
errors(CBC$ycoord) <- 100

# Create centres of abundance and linear slope coeffcients
CBC_slope_metrics <- CBC %>% group_by(AOU, Year) %>% dplyr::summarise(coa_x = sum(CBC_median*xcoord)/sum(CBC_median), coa_y = sum(CBC_median*ycoord)/sum(CBC_median)) %>% dplyr::summarise(CBC_Lat_slope = sum((Year - mean(Year))*(coa_y - mean(coa_y)))/sum((Year - mean(Year))^2), CBC_Lon_slope = sum((Year - mean(Year))*(coa_x - mean(coa_x)))/sum((Year - mean(Year))^2))

# Make the associated error its own column
CBC_slope_metrics$CBC_Lon_slope_Error <- errors(CBC_slope_metrics$CBC_Lon_slope)
CBC_slope_metrics$CBC_Lat_slope_Error <- errors(CBC_slope_metrics$CBC_Lat_slope)

#write.csv(CBC_slope_metrics, 'CBC_slope_metrics.csv', row.names = F)
CBC_slope_metrics <- read_csv('CBC_slope_metrics.csv')

BBS <- Both_final3 %>% filter(!(is.na(BBS_SD))) %>% dplyr::select(1:8, 12:16)

errors(BBS$BBS_median) <- BBS$BBS_SD
errors(BBS$xcoord) <- 100
errors(BBS$ycoord) <- 100

BBS_slope_metrics <- BBS %>% group_by(AOU, Year) %>% dplyr::summarise(coa_x = sum(BBS_median*xcoord)/sum(BBS_median), coa_y = sum(BBS_median*ycoord)/sum(BBS_median)) %>% dplyr::summarise(BBS_Lat_slope = sum((Year - mean(Year))*(coa_y - mean(coa_y)))/sum((Year - mean(Year))^2), BBS_Lon_slope = sum((Year - mean(Year))*(coa_x - mean(coa_x)))/sum((Year - mean(Year))^2))

BBS_slope_metrics$BBS_Lon_slope_Error <- errors(BBS_slope_metrics$BBS_Lon_slope)
BBS_slope_metrics$BBS_Lat_slope_Error <- errors(BBS_slope_metrics$BBS_Lat_slope)


# write.csv(BBS_slope_metrics, 'BBS_slope_metrics.csv', row.names = F)
BBS_slope_metrics <- read_csv('BBS_slope_metrics.csv')

# 2. Production of trait database ----

# Beauchamp (2011) dataset is available: https://royalsocietypublishing.org/doi/suppl/10.1098/rsbl.2011.0243

# Migration (long/short), Mig_distance (short/long), and Male_mass_g were taken from Beauchamp (2011). FlocksWithAdults is a three level trait (Solo, Flocks without adults, and Flocks with adults) was derived from Bird of the World Online literature search. Habitat_specialism_score and	Diet_specialism_score were derived from the Elton Traits database. GenFLS (modeled generation length) was taken from  Bird et al. 2020. Trend (loglinear BBS population trend for core range) was taken from USGS BBS data portal.

# Dataset 'beauchamp' created earlier was combined with the slope metrics datasets to produce:
trait_slopes <- read_csv('trait_slopes.csv')
CBC_slope_metrics <- read_csv('CBC_slope_metrics.csv')
BBS_slope_metrics <- read_csv('BBS_slope_metrics.csv')

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

trait_slopes <- join(trait_slopes, Slope_metrics[c('AOU',"BBS_dist","CBC_dist","BBS_dist_error","CBC_dist_error")])

trait_slopes <- trait_slopes %>% mutate(WithAdults = case_when(Combined_Flock_Size == 'Solo'|Combined_Cohort_Timing %in% c('Adults first','Juveniles first') ~ 'No',
                                                               Combined_Flock_Size %in% c('Small','Medium','Large') & Combined_Cohort_Timing %in% c('Concurrent') ~ 'Yes'),
                                        Flocks = case_when(Combined_Flock_Size == 'Solo' ~ 'No', 
                                                           Combined_Flock_Size %in% c('Small','Medium','Large') ~ 'Yes'),
                                        FlocksWithAdults = case_when(Combined_Flock_Size == 'Solo' ~ 'Solo',
                                                                     Flocks=='Yes' & WithAdults=='Yes'~'FlocksWithAdults',
                                                                     Flocks=='Yes' & WithAdults=='No'~'FlocksNoAdults'))


# Population trends

trends <- read.csv('BBS_trends.csv') 
trends <- trends[c('AOU','Region.Name','Trend','Relative.Abundance')]

trends <- trends %>% filter(Region.Name == 'Survey-wide                             ')
trends <- trends[c('AOU','Trend','Relative.Abundance')]

trait_slopes <- join(trait_slopes, trends)



BBS_trait_slopes <- trait_slopes[c('AOU','Genus','Species','BBS_dist','BBS_dist_error', 'Migration', 'Mig_distance','WithAdults', 'Flocks', 'FlocksWithAdults', 'Male_mass_g', 'Trend','Relative.Abundance')]
CBC_trait_slopes <- trait_slopes[c('AOU','Genus','Species','CBC_dist','CBC_dist_error', 'Migration', 'Mig_distance','WithAdults', 'Flocks', 'FlocksWithAdults', 'Male_mass_g', 'Trend','Relative.Abundance')]

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

gen_lengths <- read.csv('bird_gen_lengths.csv') 

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

names(elton_birds2)[1] <- names(BBS_trait_slopes)[15]
names(gen_lengths)[7] <- names(BBS_trait_slopes)[15]
gen_lengths2 <- gen_lengths[c('BBS_Tree_species','GenFLS','Order')]

traits_birds <- join(elton_birds2, gen_lengths2)

BBS_trait_slopes <- join(BBS_trait_slopes, traits_birds)

names(traits_birds)[1] <- 'CBC_Tree_species'
CBC_trait_slopes <- join(CBC_trait_slopes, traits_birds)


BBS_trait_slopes_flocks <- BBS_trait_slopes[c(-8,-10)]
CBC_trait_slopes_flocks <- CBC_trait_slopes[c(-8,-10)]

BBS_trait_slopes <- na.omit(BBS_trait_slopes)
CBC_trait_slopes <- na.omit(CBC_trait_slopes)

BBS_trait_slopes_flocks <- na.omit(BBS_trait_slopes_flocks)
CBC_trait_slopes_flocks <- na.omit(CBC_trait_slopes_flocks)

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
  
}

trait_slopes <- trait_slopes %>% mutate(Banding_cohort_timing2 = case_when(Overlap_timing >= 0.85 ~ 'Concurrent',Overlap_timing < 0.85 ~ 'Non-concurrent'))

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


# 3. PGLS trait analysis ----

rm(list = ls())
BBS_trait_slopes_flocks <- read.csv('BBS_trait_slopes_flocks.csv')
BBS_trait_slopes <- read.csv('BBS_trait_slopes.csv')
CBC_trait_slopes_flocks <- read.csv('CBC_trait_slopes_flocks.csv')
CBC_trait_slopes <- read.csv('CBC_trait_slopes.csv')

## 3.1 BBS (just flocks) ----

### 3.1.1 Full model ----
species_list <- unique(CBC_trait_slopes_flocks$CBC_Tree_species)
species_list <- sub(' ', '_',species_list)

BBS_trait_slopes_flocks <- BBS_trait_slopes_flocks %>% filter(AOU %in% CBC_trait_slopes_flocks$AOU)

# Load bird supertree based on Hackett's backbone 
bird_tree_hackett <- read.tree("Hackett.tre")

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
  Relative.Abundance = scale(Relative.Abundance)
)

BBSfit<-pgls.SEy((BBS_dist) ~ Migration
                 + Mig_distance
                 + Flocks
                 + Habitat_specialism_score
                 + Diet_specialism_score
                 + GenFLS
                 + Trend
                 ,data=BBS_trait_slopes_flocks2,se=(SE),tree=pruned_birds_stree,method="ML")

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
  Variable == 'Mig_distanceshort' ~ 'Short-distance migrant\n(ref: long-distance)',
  Variable == 'FlocksNo' ~ 'Solo',
  Variable == 'FlocksYes' ~ 'Flocking\n(ref: solo)',
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
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="none")

BBSresults1
BBScoefs_flocks <- BBScoefs

### 3.1.2 Backwards stepwise deletion ----

BBSfitBSD<-pgls.SEy((BBS_dist) ~ Mig_distance
                    + Habitat_specialism_score
                    + Trend
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

model.set <- MuMIn::dredge(BBSfit)

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
  variable == 'MigrationPartial' ~ 'Partial migrant',
  variable == 'Mig_distanceshort' ~ 'Short-distance migrant',
  variable == 'FlocksNo' ~ 'Solo',
  variable == 'FlocksYes' ~ 'Flocking',
  variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  variable == 'GenFLS' ~ 'Generation length',
  variable == 'Trend' ~ 'Absolute population trend'
))

### 3.2 BBS (FlockWithAdults) ----

### 3.2.1 Full model ----
species_list <- unique(CBC_trait_slopes$CBC_Tree_species)
species_list <- sub(' ', '_',species_list)

BBS_trait_slopes <- BBS_trait_slopes %>% filter(AOU %in% CBC_trait_slopes$AOU)

# Load bird supertree based on Hackett's backbone 
bird_tree_hackett <- read.tree("Hackett.tre")
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
  Relative.Abundance = scale(Relative.Abundance)
)
BBS_trait_slopes2$FlocksWithAdults <- factor(BBS_trait_slopes2$FlocksWithAdults)
BBS_trait_slopes2$FlocksWithAdults <- relevel(BBS_trait_slopes2$FlocksWithAdults, ref='Solo')
BBSfit<-pgls.SEy((BBS_dist) ~ Migration
                 + Mig_distance
                 + FlocksWithAdults
                 + Habitat_specialism_score
                 + Diet_specialism_score
                 + GenFLS
                 + Trend
                 ,data=BBS_trait_slopes2,se=(SE),tree=pruned_birds_stree,method="ML")
#summary(BBSfit)
#car::vif(BBSfit)
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
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="none")

BBSresults2
BBScoefs_FWA <- BBScoefs

### 3.2.2 Backwards stepwise deletion ----

BBSfitBSD<-pgls.SEy((BBS_dist) ~ Mig_distance
                    + FlocksWithAdults
                    + Trend
                    ,data=BBS_trait_slopes2,se=(SE),tree=pruned_birds_stree,method="ML")
#summary(BBSfitBSD)
#car::vif(BBSfitBSD)
#anova(BBSfitBSD)

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

model.set <- MuMIn::dredge(BBSfit)

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
  variable == 'MigrationPartial' ~ 'Partial migrant',
  variable == 'Mig_distanceshort' ~ 'Short-distance migrant',
  variable == 'FlocksWithAdultsFlocksNoAdults' ~ 'Age-separated flocks',
  variable == 'FlocksWithAdultsFlocksWithAdults' ~ 'Mixed-age flocks',
  variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  variable == 'GenFLS' ~ 'Generation length',
  variable == 'Trend' ~ 'Absolute population trend'
))

## 3.3 CBC (just flocks) ----

### 3.3.1 CBC full model ----
species_list <- unique(CBC_trait_slopes_flocks$CBC_Tree_species)
species_list <- sub(' ', '_',species_list)

CBC_trait_slopes_flocks <- CBC_trait_slopes_flocks %>% filter(AOU %in% CBC_trait_slopes_flocks$AOU)

# Load bird supertree based on Hackett's backbone 
bird_tree_hackett <- read.tree("Hackett.tre")

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
  Relative.Abundance = scale(Relative.Abundance)
)

CBCfit<-pgls.SEy((CBC_dist) ~ Migration
                 + Mig_distance
                 + Flocks
                 + Habitat_specialism_score
                 + Diet_specialism_score
                 + GenFLS
                 + Trend
                 ,data=CBC_trait_slopes_flocks2,se=(SE),tree=pruned_birds_stree,method="ML")
#summary(CBCfit)
#car::vif(CBCfit)
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

model.set <- MuMIn::dredge(CBCfit)

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
  variable == 'Mig_distanceshort' ~ 'Short-distance migrant\n(ref: long-distance)',
  variable == 'FlocksNo' ~ 'Solo',
  variable == 'FlocksYes' ~ 'Flocking\n(ref: solo)',
  variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  variable == 'GenFLS' ~ 'Generation length',
  variable == 'Trend' ~ 'Absolute population trend'
))

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
  Relative.Abundance = scale(Relative.Abundance)
)

CBC_trait_slopes2$FlocksWithAdults <- factor(CBC_trait_slopes2$FlocksWithAdults)
CBC_trait_slopes2$FlocksWithAdults <- relevel(CBC_trait_slopes2$FlocksWithAdults, ref='Solo')

CBCfit<-pgls.SEy((CBC_dist) ~ Migration
                 + Mig_distance
                 + FlocksWithAdults
                 + Habitat_specialism_score
                 + Diet_specialism_score
                 + GenFLS
                 + Trend
                 ,data=CBC_trait_slopes2,se=CBC_SE,tree=CBC_pruned_birds_stree,method="ML")

#summary(CBCfit)
#car::vif(CBCfit)
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

model.set <- MuMIn::dredge(CBCfit)

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
  variable == 'Mig_distanceshort' ~ 'Short-distance migrant\n(ref: long-distance)',
  variable == 'FlocksWithAdultsFlocksNoAdults' ~ 'Age-separated flocks\n(ref: solo)',
  variable == 'FlocksWithAdultsFlocksWithAdults' ~ 'Mixed age flocks\n(ref: solo)',
  variable == 'Habitat_specialism_score' ~ 'Habitat specialism score',
  variable == 'Diet_specialism_score' ~ 'Diet specialism score',
  variable == 'GenFLS' ~ 'Generation length',
  variable == 'Trend' ~ 'Absolute population trend'
))

CBCAVGplot <- ggplot(CBCAVGcoefs, aes(x=reorder(Variable, value), y=value, colour=Significance, shape = Avg)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6) +
  geom_point(posiiton=position_dodge(), position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.25)) +
  labs(y='Coefficient', x='', shape='Average') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="right")

CBCAVGplot


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

mapplot <- ggplot() +  
  geom_sf(data=world2, size=0.2) +
  geom_sf(data=lakes, fill="white", size=0.1) +
  coord_sf(xlim = c(-122.5,-73.5), ylim = c(25,45), expand = TRUE) +
  geom_path(CBC_movements3, mapping=aes(x=Long, y=Lat, group=AOU, colour=FlocksWithAdults,alpha=CBC_dist_error),linewidth=.5,arrow=arrow(length=unit(0.06,'cm'), type='closed')) +
  theme_void() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.1), legend.position = "bottom", legend.text = element_text(size=7), legend.title = element_text(size=7), plot.background = element_blank(), panel.background = element_rect(fill='white')) +
  scale_color_manual(values=cols, guide = 'none') +
  scale_alpha_continuous(range = c(0.2, 1), n.breaks=4, breaks=seq(0.001, 0.004, 0.001), labels=c('<0.001', seq(0.002, 0.004, 0.001))) +
  labs(fill = "", colour = "", alpha='Reciprocal error:') +
  guides(alpha = guide_legend(nrow=1,byrow=TRUE))

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
  geom_strip(taxa1='Dendrocygna autumnalis',taxa2= 'Anas americana', barsize=0.3, color='black', label="Anseriformes", offset=4, offset.text = 5, fontsize=1.75) +
  geom_strip('Dendroica townsendi', 'Myiarchus cinerascens', barsize=0.3, color='black', label="Passeriformes", offset=4,offset.text = 7, fontsize=1.75) +
  geom_strip('Falco sparverius', 'Falco peregrinus', barsize=0.3, color='black', label="Falconiformes", offset=4, offset.text = 55, fontsize=1.75) +
  geom_strip('Larus marinus', 'Tringa flavipes', barsize=0.3, color='black', label="Charadriiformes", offset=4,offset.text = 10, hjust=1, fontsize=1.75) +
  geom_strip(72, 72, barsize=0.3, color='black', label="Coraciiformes", offset=4,offset.text = 50, fontsize=1.75, hjuust=-0.1) +
  geom_strip(71, 71, barsize=0.3, color='black', label="Strigiformes", offset=4,offset.text = 45, fontsize=1.75) +
  geom_strip('Pandion haliaetus', 'Buteo regalis', barsize=0.3, color='black', label="Accipitriformes", offset=4,offset.text = 55, hjust=0.05, fontsize=1.75) +
  geom_strip('Zenaida asiatica', 'Zenaida asiatica', barsize=0.3, color='black', label="Columbiformes", offset=4,offset.text = 6, fontsize=1.75) +
  geom_strip('Plegadis falcinellus', 'Plegadis falcinellus', barsize=0.3, color='black', label="Pelecaniformes", offset=4,offset.text = 11,hjust=0.1, fontsize=1.75) +
  geom_strip('Rallus limicola', 'Grus canadensis', barsize=0.3, color='black', label="Gruiformes", offset=4,offset.text = 16, fontsize=1.75, hjust=0.2) +
  geom_strip('Archilochus alexandri', 'Selasphorus rufus', barsize=0.3, color='black', label="Apodiformes", offset=4,offset.text = 7, fontsize=1.75) +
  labs(color='') +
  scale_color_manual(values=cols2,labels=c('Age-separated\nflocks','Mixed-age\nflocks','Solo')) +
  theme(legend.text = element_text(size=7),legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10), legend.background = element_blank(), legend.position = c(1,0.5))


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


# Figure 2 ----

BBSresults1 <- BBSresults1 + theme(axis.text = element_text(size=7), axis.title = element_text(size=9), plot.title = element_text(size=11), plot.subtitle = element_text(size=10))
BBSresults2 <- BBSresults2 + theme(axis.text = element_text(size=7), axis.title = element_text(size=9), plot.title = element_text(size=11), plot.subtitle = element_text(size=10))
NonBreedingresults1 <- NonBreedingresults1 + theme(axis.text = element_text(size=7), axis.title = element_text(size=9), plot.title = element_text(size=11), plot.subtitle = element_text(size=10))
NonBreedingresults2 <- NonBreedingresults2 + theme(axis.text = element_text(size=7), axis.title = element_text(size=9), plot.title = element_text(size=11), plot.subtitle = element_text(size=10))

(main_fig <- ggdraw() +
    draw_plot(BBSresults1, x=0,y=0, width=0.5,height=1) +
    draw_plot(BBSresults2, x=0.5,y=0, width=0.5,height=1))


BBScoefs_flocks$Model <- 'Model A'
BBScoefs_FWA$Model <- 'Model B'

BBScoefs_all <- rbind(BBScoefs_flocks,BBScoefs_FWA)


ggplot(BBScoefs_all, aes(x=reorder(Variable,Coef), y=Coef, colour=Significance, shape=Model)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6, position=position_dodge(width=0.75)) +
  geom_point(alpha=0.75, position=position_dodge(width=0.75)) +
  geom_point(position=position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.75)) +
  labs(y='Coefficient', x='', title='Breeding') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="none")


CBCcoefs_flocks$Model <- 'Model A'
CBCcoefs_FWA$Model <- 'Model B'

CBCcoefs_all <- rbind(CBCcoefs_flocks,CBCcoefs_FWA)


ggplot(CBCcoefs_all, aes(x=reorder(Variable,Coef), y=Coef, colour=Significance, shape=Model)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.6, position=position_dodge(width=0.75)) +
  geom_point(alpha=0.75, position=position_dodge(width=0.75)) +
  geom_point(position=position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.75)) +
  labs(y='Coefficient', x='', title='Non-breeding') +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=10), legend.position="none")

BBScoefs_all2 <- BBScoefs_all
CBCcoefs_all2 <- CBCcoefs_all
BBScoefs_all2$Season <- 'Breeding'
CBCcoefs_all2$Season <- 'Non-breeding'

coefs_all <- rbind(BBScoefs_all2, CBCcoefs_all2)

results_plot <- ggplot(coefs_all, aes(x=reorder(Variable,Coef), y=Coef, colour=Significance, shape=Model)) +
  geom_hline(aes(yintercept=0), linetype='dashed', alpha=.4) +
  geom_point(size=2,position=position_dodge(width=0.3)) +
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..), width=0, position=position_dodge(width=0.3)) +
  theme_half_open() +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  facet_grid(~ Season) +
  labs(y='Coefficient') +
  theme(axis.title.x = element_text(size=13), axis.text.y = element_text(size=10), legend.position="bottom", axis.title.y = element_blank(), strip.background = element_blank(), legend.title = element_blank(), legend.text = element_text(size=10), plot.background = element_rect(fill='white'))

results_plot
