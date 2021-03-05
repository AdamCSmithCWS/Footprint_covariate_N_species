####### Data load and preparation

library(patchwork)

library(tidyverse)
library(bbsBayes)
library(sf)

all = read.csv("data/combined_2004_2019_Canadian_trends_and_intercepts.csv")
#weird header missmatch fix
names(all) <- c(names(all)[-1],"geometry2")

sp_consider = unique(all$species)
rm("all")



# spatial data load -------------------------------------------------------

laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system

locat = system.file("maps",
                    package = "bbsBayes")
map.file = "BBS_USGS_strata" # the most standard BBS stratification, used here just to provide some spatial limits

# reading in the strata spatial file (GIS shapefile)
strata_map = read_sf(dsn = locat,
                     layer = map.file)
strata_map = st_transform(strata_map,crs = laea) #reprojecting the geographic coordinate file to an equal area projection



# bbsBayes data setup ----------------------------------------------------------


stratification = "bbs_usgs"

model = "slope"

bbs_strat = stratify(by = stratification)

# extracting the sampling event stuff --------------------------

bbs_samp_event <- bbs_strat$route_strat
head(bbs_samp_event,2) #shows the top 2 rows of the data frame
# countrynum statenum Route RouteName Active Latitude Longitude BCR RouteTypeID RouteTypeDetailID
# 7         124       11     1 CHEMAINUS      0 48.78529 -123.5991   5           1                 1
# 16        124       11     1 CHEMAINUS      0 48.78529 -123.5991   5           1                 1
# RouteDataID RPID Year Month Day    ObsN TotalSpp StartTemp EndTemp TempScale StartWind EndWind
# 7      6356273  101 2014     6  15 1160088       48      10      14           C         0       2
# 16     6352859  101 2013     7   7 1160088       43      14      16           C         0       0
# StartSky EndSky StartTime EndTime Assistant QualityCurrentID RunType            State St_Abrev
# 7         0      1       441     846      1                   1       1 British Columbia       BC
# 16        0      0       441     900      1                   1       1 British Columbia       BC
# Country strat_name rt.uni  rt.uni.y
# 7       CA    CA-BC-5   11-1 11-1-2014
# 16      CA    CA-BC-5   11-1 11-1-2013

# removing the non-Canadian data ------------------------------------------

names_strata <- get_composite_regions(strata_type = stratification) ## bbsBayes function that gets a list of the strata names and their composite regions (provinces, BCRs, etc.)

CAN_strata_keep <- names_strata[which(names_strata$national == "CA"),"region"] # character vector of the strata names to remove from data




# drop all other sampling event data for which we do not have the observer's age
bbs_samp_event <- bbs_samp_event[which(bbs_samp_event$strat_name %in% CAN_strata_keep),]
nrow(bbs_samp_event)
#[1] 17671 - there are 17671 BBS route's run over the 53 year history of the BBS, in Canada




# bird observations -------------------------------------------------------


bbs_birds <- bbs_strat$bird_strat
head(bbs_birds) # shows the top 6 rows of the data frame
# statenum Route countrynum RouteDataID RPID Year  AOU Count10 Count20 Count30 Count40 Count50
# 1       11     1        124     6227579  101 1997 5880      12       8       3       3       2
# 2       11     1        124     6227579  101 1997 5671       6       0       0       0       0
# 3       11     1        124     6227579  101 1997 6850       0       2       0       1       0
# 4       11     1        124     6214132  101 1991 6882       0      10      12       0       0
# 5       11     1        124     6344624  101 2012 3160       0       0       1       0       0
# 6       11     1        124     6214132  101 1991 4860       2       4       3       0       0
# StopTotal SpeciesTotal BCR rt.uni  rt.uni.y
# 1        16           28   5   11-1 11-1-1997
# 2         2            6   5   11-1 11-1-1997
# 3         3            3   5   11-1 11-1-1997
# 4         2           22   5   11-1 11-1-1991
# 5         1            1   5   11-1 11-1-2012
# 6         6            9   5   11-1 11-1-1991
# AOU = unique numerical species identifier, but not very user-friendly
# rt.uni.y = unique identifier of the route and year combination - matches with same column in bbs_samp_event
# SpeciesTotal = total number of individuals observed during the survey


# merge in the english species names --------------------------------------
bbs_species <- bbs_strat$species_strat
head(bbs_species)
# seq   aou                      english                             french
# 1   6 01770 Black-bellied Whistling-Duck          Dendrocygne à ventre noir
# 2   7 01780       Fulvous Whistling-Duck                  Dendrocygne fauve
# 3   8 01760                Emperor Goose                       Oie empereur
# 4   9 01690       Snow Goose (all forms) Oie des neiges (toutes les formes)
# 5  10 01691      (Blue Goose) Snow Goose       Oie des neiges (forme bleue)
# 6  11 01700                 Ross's Goose                        Oie de Ross
#                          spanish        order   family       genus                  species sp.bbs
# 1         Dendrocygna autumnalis Anseriformes Anatidae Dendrocygna               autumnalis   1770
# 2            Dendrocygna bicolor Anseriformes Anatidae Dendrocygna                  bicolor   1780
# 3                Anser canagicus Anseriformes Anatidae       Anser                canagicus   1760
# 4             Anser caerulescens Anseriformes Anatidae       Anser             caerulescens   1690
# 5 Anser caerulescens (blue form) Anseriformes Anatidae       Anser caerulescens (blue form)   1691
# 6                   Anser rossii Anseriformes Anatidae       Anser                   rossii   1700
bbs_species <- bbs_species[,c("sp.bbs","english")]
bbs_species <- rename(bbs_species,AOU = sp.bbs) #chaning the name of the aou column to match the column name in bbs_birds



# mering the species names with the observation data ----------------------
bbs_birds <- dplyr::left_join(bbs_birds,bbs_species,by = "AOU")


#### remainder of columns in bbs_birds are not necessary and so can be dropped for now
bbs_birds <- select(bbs_birds,
                    rt.uni.y,SpeciesTotal,english)

#unique list of routes by year sampling events
rts_inc <- unique(bbs_samp_event$rt.uni.y)

#dropping all bird observations from routes where we don't have observer age info
bbs_birds <- filter(bbs_birds,rt.uni.y %in% rts_inc,
                    english %in% sp_consider)

# list of the species observed on the routes where we have observer age
species_inc <- unique(bbs_birds$english)

if(length(species_inc) != length(sp_consider)){stop("Species Missing!!!")}

# Adding in the zeros -----------------------------------------------------

#make a complete set of sampling events by species
full_event_species <- expand.grid(rt.uni.y = rts_inc,
                                  english = sp_consider)

#join the complete list with the actual list of sampling events (just a trick to generate all the extra columns from the sampling event)
full_bbs <- left_join(full_event_species,bbs_samp_event,
                      by = "rt.uni.y")

#join the full list of sampling events by species with the actual observations
full_bbs <- left_join(full_bbs,bbs_birds,
                      by = c("rt.uni.y","english"))

length(which(is.na(full_bbs$SpeciesTotal)))
#[1] 2130999 there are 2.13 Million combinations of surveys by species that should be zeros (i.e., where the species was not observed)

#replace these NA values with zeros
full_bbs[which(is.na(full_bbs$SpeciesTotal)),"SpeciesTotal"] <- 0 
length(which(is.na(full_bbs$SpeciesTotal)))
#[1] 0  

length(which(full_bbs$SpeciesTotal == 0))
#[1] 2130999
#there are 2130999 Million zeros

#renaming the SpeciesTotal column to something simpler
full_bbs <- rename(full_bbs,count = SpeciesTotal)



save(list = c("full_bbs"),file = "data/full_bbs.RData")


