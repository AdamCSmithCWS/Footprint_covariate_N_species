### building a spatial route-level trend model with footprint covariate on abundance and slope


library(bbsBayes)
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE,javascript = FALSE)
library(shinystan)
library(sf)
library(spdep)
library(patchwork)
 library(ggforce)
 library(tidybayes)
source("functions/mungeCARdata4stan.R") ## function to modify the BUGS formatted spatial neighbourhood data to the required format for the Stan iCAR model
#source("functions/prepare-jags-data-alt.R") ## small alteration of the bbsBayes function
#source("functions/GAM_basis_function.R") ## simple function to produce the GAM basis function that enters the Stan model as data


firstYear = 2004 # first year to consider
# means 15 year trends


# Species selection -------------------------------------------------------
# species for which spatial 15 trend model has been run
all = read.csv(paste0("data/combined_",firstYear,"_2019_Canadian_trends_and_intercepts.csv"))
#weird header missmatch fix
names(all) <- c(names(all)[-1],"geometry2")

sp_consider = unique(all$species)
rm("all")


# load species groupings from State of Canada's Birds
socb = read.csv("data/SOCB data supplement.csv")

groups = c("Grassland.birds",
           "Forest.birds",
           "Other.birds",
           "Aerial.insectivores",
           "suburban",
           "other.wetland.birds")



# group loop if desired--------------------------------------------------------------

group_name = groups[4]


# species selection of trend data -----------------------------------------


socb_c = which(grepl(names(socb),pattern = group_name))

sp_sel1 = socb[which(socb[,socb_c] == "Included in group"),"species"]

species_list = sp_sel1[which(sp_sel1 %in% sp_consider)]

nspecies = 7

species_list = species_list[c(1:3,5,7,8,9)]

# nspecies = length(species_list)



# bird data load ----------------------------------------------------------

#reformat full_bbs into species count by observation event matrix ---------------------------------------------------------------

load("data/full_bbs.RData") ## this is the full set of BBS observations for observers with age information


species_df = data.frame(english = species_list,
                        sp = 1:nspecies)



sel_bbs = full_bbs %>% filter(Year >= firstYear,
                              english %in% species_list,
                              Country == "CA") %>% # selecting the years and species
  arrange(english,statenum,strat_name,Route,Year) #sorting the dataframe




# Canadian current footprint data at route-level --------------------------


load("data/compiled_footprint_data.RData")


buf_sel = buffer_sizes[3] #selecting the 4 km buffer

# fp_components
# [1] "cumulative"                   "built"                       
# [3] "crop"                         "dam_and_associated_reservoir"
# [5] "forestry_harvest"             "mines"                       
# [7] "nav_water"                    "night_lights"                
# [9] "oil_gas"                      "pasture"                     
# [11] "population_density"           "rail"                        
# [13] "roads" 

preds <- fp_components[c(1)]

cls_sel <- paste(preds,buf_sel,"mean",sep = "_")
cls_sel_i <- c("rt.uni",cls_sel)



fp_can_sel <- fp_can_by_route %>% select(all_of(cls_sel_i)) %>% na.exclude() %>% 
  distinct() %>% 
  mutate(pred = across(all_of(cls_sel),scale))

nm1 = (2+length(cls_sel))
nm2 = (nm1+(length(cls_sel)-1))
w_preds <- nm1:nm2
#names(fp_can_sel)[c(nm1:nm2)] <- paste0("pred_",1:length(cls_sel))

ncovs = length(w_preds)

rts_w_cov <- unique(fp_can_sel$rt.uni)

## weird clunky way of removing the formatting from the fp_can_sel tibble
fp_data <- data.frame(rt.uni = unlist(fp_can_sel$rt.uni),
                      pred = unlist(fp_can_sel[,w_preds]),
                      row.names = NULL)

#identify routes where at least one of the species has been observed at least once
## and where we have covariate data
routes_incl <- sel_bbs %>% 
  mutate(gt_z = ifelse(count > 0,1,0)) %>% 
  group_by(rt.uni) %>% 
  summarise(sum_gt_z = sum(gt_z),
            non_zero = ifelse(sum_gt_z > 0,1,0)) %>% 
  filter(non_zero == 1)
#Drop routes where none of the species have ever appeared
sel_bbs <- sel_bbs %>% 
  filter(rt.uni %in% unlist(routes_incl$rt.uni)) %>% 
  filter(rt.uni %in% rts_w_cov) %>% 
  mutate(observer = as.integer(factor(ObsN)),
         year = Year-(firstYear-1),
         rts_cov = as.integer(factor(rt.uni))) %>% 
  left_join(.,species_df) %>% 
  left_join(.,fp_data) %>% 
  arrange(sp)

nobservers = max(sel_bbs$observer)


rt_cov_mat = sel_bbs %>% 
  select(rt.uni,rts_cov,pred) %>% 
  distinct() %>% 
  arrange(rts_cov)

nrts_cov = max(rt_cov_mat$rts_cov)
### this covariate matrix includes the covariate values for every route included and a unique
### route indicator common across all species

nyears = max(sel_bbs$year)
fixedyear = floor(nyears/2) #middle year, used to center the years values so that the 
# alphas (intercepts) represent the mean abundance in the middle of the time-series






# spatial neighbourhood define --------------------------------------------
laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system

locat = system.file("maps",
                    package = "bbsBayes")
map.file = "BBS_USGS_strata"

strata_map = read_sf(dsn = locat,
                     layer = map.file)
strata_map = st_transform(strata_map,crs = laea)



rt_cov_dat <- as.numeric(rt_cov_mat[,w_preds])

stan_data = list(nspecies = nspecies,
                fixedyear = fixedyear,
                nyears = nyears,
                nobservers = nobservers,
                nrts_cov = nrts_cov,
                #ncovs = ncovs,
                rt_cov_dat = rt_cov_dat)




# Species data set-up -----------------------------------------------------

species_dfsplit <- group_split(sel_bbs, sp, .keep = TRUE)

#build the correct species-specific dataframes
species_dataframes <- vector("list",length = nspecies)
# nrts_cov_sp <- 0


#### build a species-specific covariate matrix that matches the route indicators
#### for each species. 
#### necessary because the covariate effects have to be combined with the intercepts
#### and slope parameters 

ncounts <- 0

for(i in 1:nspecies){
  spec = species_df[which(species_df$sp == i),"english"]
  tmp <- species_dfsplit[[i]]
  nobstmp = tmp %>% group_by(rt.uni) %>% 
    summarise(nobs = ifelse(sum(count)>0,1,0)) %>% 
    filter(nobs > 0)
  tmp <- filter(tmp,rt.uni %in% unlist(nobstmp$rt.uni))
  tmp$routeF <- as.integer(factor(tmp$rt.uni)) 

  #   tmp1 = tmp$rts_cov
  #   nrts_cov_sp[i] = length(tmp1)
  # rts_cov_sp[i,1:nrts_cov_sp[i]] <- tmp1
  # 
  

  species_dataframes[[i]] <- tmp
  
  write.csv(sel_bbs,file = paste0("Data/",group_name,"_",spec,"_input_data_frame.csv")) # exporting the sel_bbs dataframe for future-me
  






# Spatial boundaries set up --------------------

# the iCAR (intrinsic Conditional AutoRegressive) spatial model uses neighbourhood structure
# to share information on abundance and trend (intercept and slope) among BBS routes
# 


real_strata_map = filter(strata_map,ST_12 %in% unique(tmp$strat_name)) #dropping strata with no data for this study

strata_bounds <- st_union(real_strata_map) #union to provide a simple border of the realised strata
strata_bounds_buf = st_buffer(strata_bounds,dist = 300000) #buffering the realised strata by 300km




# spatial neighbourhood define --------------------------------------------

route_map = tmp %>% 
  select(rt.uni,
         routeF,
         strat_name,
         St_Abrev,
         Latitude,
         Longitude) #selecting the data on the routes and their start coordinates

route_map <- unique(route_map) #data frame with just one row for each route
  
# reconcile duplicate spatial locations for routes -----------------------------------
# adhoc way of separating different routes with the same starting coordinates
# this shifts the starting coordinates of teh duplicates by ~1.5km to the North East 
# ensures that the duplicates have a unique spatial location, but remain very close to
# their original location and retain the correct neighbourhood relationships
# these duplicates happen when a "new" route is established because some large proportion
# of the end of a route is changed, but the start-point remains the same

dups = which(duplicated(route_map[,c("Latitude","Longitude")]))
while(length(dups) > 0){
  route_map[dups,"Latitude"] <- route_map[dups,"Latitude"]+0.01 #=0.01 decimal degrees ~ 1km
  route_map[dups,"Longitude"] <- route_map[dups,"Longitude"]+0.01 #=0.01 decimal degrees ~ 1km
  dups = which(duplicated(route_map[,c("Latitude","Longitude")]))
  
}
dups = which(duplicated(route_map[,c("Latitude","Longitude")])) 
if(length(dups) > 0){stop(paste(spec,"ERROR - At least one duplicate route remains"))}

route_map = st_as_sf(route_map,coords = c("Longitude","Latitude"))
st_crs(route_map) <- 4269 #NAD83 coordinate reference system (crs) commonly used by US federal agencies to store geographic coordinate information

route_map = st_transform(route_map,crs = laea) #reproject to the same equal area projection

# Voronoi polygons from route locations -----------------------------------
# this creates one polygon for each route location, where the polygons define the
# space closest to each route and therefore the neighbouring polygons define the
# neighbourhood relationships for each route
box <- st_as_sfc(st_bbox(route_map)) # a bounding box to limit the extent of the neighbourhood polygons

v <- st_cast(st_voronoi(st_union(route_map), envelope = box))

vint = st_sf(st_cast(st_intersection(v,strata_bounds_buf),"POLYGON")) #making the voronoi polygons
vintj = st_join(vint,route_map,join = st_contains) #joining the polygons to teh routes
vintj = arrange(vintj,routeF) #sorting teh polygons by the numerical route indicator
# ensures that route 1 is route 1 in both the spatial and non-spatial data


nb_db = poly2nb(vintj,row.names = vintj$route,queen = FALSE) # creates the lists of neighbours for each route


#

# plotting the neighbourhoods to check ------------------------------------

cc = suppressWarnings(st_coordinates(st_centroid(vintj)))

ggp = ggplot(data = route_map)+
  geom_sf(data = vintj,alpha = 0.3)+ 
  geom_sf(aes(col = St_Abrev))+
 geom_sf_text(aes(label = routeF),size = 3,alpha = 0.3)+
  theme(legend.position = "none")
pdf(file = paste0("output/",group_name,"_",spec," route maps.pdf"),
    width = 11,
    height = 8.5)
plot(nb_db,cc,col = "pink")
text(labels = rownames(cc),cc ,pos = 2)
print(ggp)
dev.off()
save(list = c("route_map",
              "vintj",
              "nb_db",
              "cc"),
     file = paste0("data/speciesRouteData/",group_name,"_",spec,"_",nspecies,"_route_data.RData"))

## stop here and look at the maps (2 pages)
## in the first page each route location is plotted as a point and all neighbours are linked by red lines 
## in the second page all of the voronoi polygons with their route numbers are plotted
### assuming the above maps look reasonable 

nb_info = spdep::nb2WB(nb_db) #function int eh spdep package designed to produce neighbourhood data for GeoBUGS (and older Bayesian analysis program)


### re-arrange GEOBUGS formated nb_info into appropriate format for Stan model
car_stan_dat <- mungeCARdata4stan(adjBUGS = nb_info$adj,
                                  numBUGS = nb_info$num)

stan_data[[paste0("N_edges",i)]] <- car_stan_dat$N_edges
stan_data[[paste0("node1",i)]] <- car_stan_dat$node1
stan_data[[paste0("node2",i)]] <- car_stan_dat$node2

stan_data[[paste0("nroutes",i)]] <- max(tmp$routeF)
stan_data[[paste0("ncounts",i)]] <- nrow(tmp)
stan_data[[paste0("count",i)]] <- tmp$count
stan_data[[paste0("year",i)]] <- tmp$year
stan_data[[paste0("route",i)]] <- tmp$routeF
stan_data[[paste0("route_cov",i)]] <- tmp$rts_cov
stan_data[[paste0("observer",i)]] <- tmp$observer

ncounts = ncounts+nrow(tmp)

}#end species loop

stan_data[["ncounts"]] <- ncounts



if(length(stan_data) != 7+nspecies*10){stop("Something missing from stan_data object")}
ch_ns <- paste0("ncounts",1:nspecies)
ncs1 = 0
for(j in ch_ns){
  ncs1 = ncs1+stan_data[[j]]
}
if(ncs1 != ncounts){stop(paste("Incorrect counting of counts ncs1=",ncs1, "and ncounts=",ncounts))}



#  ------------------------------

mod.file = paste0("models/iCAR_",nspecies,"_species_cov.stan")

# parms = c("sdnoise",
#           "sdobs",
#           "sdbeta",
#           "alpha",
#           "sdalpha",
#           "BETA",
#           "ALPHA",
#           "beta",
#           "log_lik")

parms = c("log_lik",
          "sdnoise",
          paste0("alpha",1:nspecies),
          paste0("beta",1:nspecies),
          "obs",
          "ALPHA",
          "BETA",
          "A_cov",
          "B_cov")

## compile model
model = stan_model(file=mod.file)

print(paste(group_name,nspecies,"cov"))

## run sampler on model, data
stanfit <- sampling(model,
                          data=stan_data,
                          verbose=TRUE, refresh=50,
                          chains=3, iter=1100,
                          warmup=800,
                          cores = 3,
                          pars = parms,
                          control = list(adapt_delta = 0.9,
                                         max_treedepth = 15))


save(list = c("stanfit","stan_data","vintj","route_map","sel_bbs","species_list","species_df"),
     file = paste0("output/",group_name,"_",nspecies,"_cov_adj_output.RData"))


launch_shinystan(stanfit) 



  
#   #  ------------------------------
#   
  mod.file = paste0("models/iCAR_",nspecies,"_species_naive.stan")
  
  parms = c("log_lik",
            "sdnoise",
            paste0("alpha",1:nspecies),
            paste0("beta",1:nspecies),
            "obs",
            "ALPHA",
            "BETA")

  ## compile model
  model = stan_model(file=mod.file)


stan_data_naive <- stan_data


stan_data_naive[["age"]] <- NULL
stan_data_naive[["nages"]] <- NULL
stan_data_naive[["age_basispred"]] <- NULL
stan_data_naive[["nknots_age"]] <- NULL

print(paste(group_name,nspecies,"naive"))
  ## run sampler on model, data
  stanfit_naive <- sampling(model,
                           data=stan_data,
                           verbose=TRUE, refresh=50,
                           chains=4, iter=1100,
                           warmup=800,
                           cores = 4,
                           pars = parms,
                           control = list(adapt_delta = 0.9,
                                          max_treedepth = 15))


  save(list = c("stanfit_naive","stan_data_naive","vintj","route_map","sel_bbs","species_list","species_df"),
       file = paste0("output/",group_name,"_",nspecies,"_naive_output.RData"))

# 
#   # prior predictive check --------------------------------------------------
#   

#   
#   launch_shinystan(stanfit_naive)





####### see scipt "Results_plotting_initial.R"




# Plotting ----------------------------------------------------------------
# 
# 
# load(paste0("output/",group_name,"_",nspecies,"_age_adj_output.RData"))
# load(paste0("output/",group_name,"_",nspecies,"_naive_output.RData"))
# 
# 
# 
# library(tidybayes)
# library(scales)
# 
# nobs_scale = 1
# nobs = select(.data = sel_bbs,
#                  observer_age) %>% 
#   group_by(observer_age) %>% 
#   slice_sample(prop = nobs_scale)
# 
# age_samples = gather_draws(stanfit,age_pred[a])
# 
# age_effect = age_samples %>% group_by(a) %>% 
#   summarise(age_effect = exp(mean(.value)),
#             lci = exp(quantile(.value,0.025)),
#             uci = exp(quantile(.value,0.975)),
#             sd = (sd(.value)),
#             prec = 1/var(.value),
#             .groups = "keep")
# 
# 
# age_effect$age = age_effect$a + min(sel_bbs$observer_age)
# age_plot = ggplot(data = age_effect,aes(x = age,y = age_effect))+
#   geom_line(colour = viridis_pal(option = "plasma")(1))+
#   geom_ribbon(aes(x = age,ymin = lci,ymax = uci),alpha = 0.3,fill = viridis_pal(option = "plasma")(1))+
#   geom_dotplot(data = nobs,mapping = ggplot2::aes(x = observer_age),drop = TRUE,binaxis = "x", 
#                stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = FALSE,
#                fill = grDevices::grey(0.6),colour = grDevices::grey(0.6),alpha = 0.2,dotsize = 0.3)+
#   coord_cartesian(ylim = c(NA,2))+
#   labs(x = "Observer Age",
#        y = "Relative proportion of individuals counted",
#        title = "Multiplicative effect of observer age on BBS observations",
#        subtitle = paste("Species with",group_name,"frequency auditory cues"),
#        caption = paste("Species included:",paste(species_list,collapse = ", ")))
# print(age_plot)
# 
# 
# pdf(paste0("Figures/age_effect_",group_name,".pdf"),
#     width = 11,
#     height = 8.5)
# print(age_plot)
# dev.off()
# 
# 
# 
# 
# # plot the hyperparameters from the two models
# BETA_samples = gather_draws(stanfit,BETA[s])
# BETA_a = BETA_samples %>% group_by(s) %>%
#   summarise(Trend = 100*(exp(mean(.value))-1),
#             lci_Trend = 100*(exp(quantile(.value,0.025))-1),
#             uci_Trend = 100*(exp(quantile(.value,0.975))-1),
#             sd_Trend = 100*((sd(.value))-1),
#             .groups = "keep")
# BETA_a$species = species_list[BETA_a$s]
# BETA_a$model = "Age_adjusted"
# 
# 
# BETA_samples_n = gather_draws(stanfit_naive,BETA[s])
# BETA_n = BETA_samples_n %>% group_by(s) %>%
#   summarise(Trend = 100*(exp(mean(.value))-1),
#             lci_Trend = 100*(exp(quantile(.value,0.025))-1),
#             uci_Trend = 100*(exp(quantile(.value,0.975))-1),
#             sd_Trend = 100*((sd(.value))-1),
#             .groups = "keep")
# BETA_n$species = species_list[BETA_n$s]
# BETA_n$model = "Naive"
# 
# BETAs = bind_rows(BETA_a,BETA_n)
# 
# B_plot = ggplot(data = BETAs,aes(x = species,y = Trend,colour = model,group = model))+
#   geom_pointrange(aes(ymin = lci_Trend, ymax = uci_Trend),position = position_dodge(width = 0.25))+
#   coord_flip()
# 
# pdf(paste0("figures/","comparing_hyperparamter_trends.pdf"))
# print(B_plot)
# dev.off()
# 
# # 
# 
# 
# # compare the trends at each route ----------------------------------------
# 
# 
# laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system
# 
# locat = system.file("maps",
#                     package = "bbsBayes")
# map.file = "BBS_USGS_strata"
# 
# strata_map = read_sf(dsn = locat,
#                      layer = map.file)
# strata_map = st_transform(strata_map,crs = laea)
# 
# strata_map = strata_map %>% filter(COUNTRY == "CA")
# 
# 
# ####
# # add trend and abundance ----------------------------------------
# 
# beta_samples = gather_draws(stanfit,beta[r,s])
# alpha_samples = gather_draws(stanfit,alpha[r,s])
# 
# slopes = beta_samples %>% group_by(r,s) %>% 
#   summarise(b = mean(.value),
#             lci = quantile(.value,0.025),
#             uci = quantile(.value,0.975),
#             sd = sd(.value),
#             .groups = "keep")
# 
# 
# ints = alpha_samples %>% group_by(r,s) %>% 
#   summarise(int = mean(exp(.value)),
#             int_lci = quantile(exp(.value),0.025),
#             int_uci = quantile(exp(.value),0.975),
#             int_sd = sd(exp(.value)),
#             .groups = "keep")
# 
# slopes = full_join(slopes,ints)
# slopes$species = species_list[slopes$s]
# slopes$model = "Age_adjusted"
# 
# 
# beta_samples_n = gather_draws(stanfit_null,beta[r,s])
# alpha_samples_n = gather_draws(stanfit_null,alpha[r,s])
# 
# slopes_n = beta_samples_n %>% group_by(r,s) %>% 
#   summarise(b = mean(.value),
#             lci = quantile(.value,0.025),
#             uci = quantile(.value,0.975),
#             sd = sd(.value),
#             .groups = "keep")
# ints_n = alpha_samples_n %>% group_by(r,s) %>% 
#   summarise(int = mean(exp(.value)),
#             int_lci = quantile(exp(.value),0.025),
#             int_uci = quantile(exp(.value),0.975),
#             int_sd = sd(exp(.value)),
#             .groups = "keep")
# 
# slopes_n = full_join(slopes_n,ints_n)
# 
# slopes_n$species = species_list[slopes_n$s]
# slopes_n$model = "naive"
# 
# 
# 
# 
# 
# # lists of the species and model specific maps ----------------------------
# 
# 
# maps_age = vector(mode = "list",length = length(species_list))
# maps_n = vector(mode = "list",length = length(species_list))
# 
# 
# for(sps in 1:max(slopes$s)){
#   spsn = species_list[sps]
#   tmp <- slopes %>% filter(species == spsn)
#   maps_age[[sps]] = inner_join(route_map,tmp,by = c("routeF" = "r"))
#   
#   tmp <- slopes_n %>% filter(species == spsn)
#   maps_n[[sps]] = inner_join(route_map,tmp,by = c("routeF" = "r"))
#   
# }
# 
# # add mapping of trends ---------------------------------------------------
# 
# # plot the route-level trends maps compare adjacent maps
# 
# 
# for(sps in 1:max(slopes$s)){
#   spsn = species_list[sps]
#   
# breaks <- c(-0.07, -0.04, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.04, 0.07)
# labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
# labls = paste0(labls, " slope")
# 
# route_map_age <- maps_age[[sps]]
# route_map_naive <- maps_n[[sps]]
# 
# route_map_age$Tplot <- cut(route_map_age$b,breaks = c(-Inf, breaks, Inf),labels = labls)
# route_map_naive$Tplot <- cut(route_map_naive$b,breaks = c(-Inf, breaks, Inf),labels = labls)
# map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
#                  "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
# names(map_palette) <- labls
# 
# 
# tmap_age = ggplot(route_map_age)+
#   geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
#   geom_sf(aes(colour = Tplot,size = int))+
#   scale_colour_manual(values = map_palette, aesthetics = c("colour"),
#                       guide = guide_legend(reverse=TRUE),
#                       name = paste0("slope\n",1999,"-",2019))+
#   labs(title = paste(spsn,"Age-adjusted trends by route (size = abundance)"))
# 
# tmap_naive = ggplot(route_map_naive)+
#   geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
#   geom_sf(aes(colour = Tplot,size = int))+
#   scale_colour_manual(values = map_palette, aesthetics = c("colour"),
#                       guide = guide_legend(reverse=TRUE),
#                       name = paste0("slope\n",1999,"-",2019))+
#   labs(title = paste(spsn,"Naive trends by route (size = abundance)"))
# 
# 
# pdf(file = paste0("figures/",spsn,"trend_maps_age_adjusted.pdf"),
#     width = 11,
#     height = 8.5)
# 
# print(tmap_age)
# print(tmap_naive)
# 
# dev.off()
# 
# 
# }
# 
# # 
# # write.csv(route_map_out,
# #           file = paste0("output/",species," ",firstYear," ",lastYear,"_Canadian_trends_and_intercepts.csv"))
# # 
# 




# plot the comparisons between models -------------------------------------




