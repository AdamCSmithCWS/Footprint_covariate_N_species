
# Plotting ----------------------------------------------------------------



library(patchwork)
library(tidyverse)
library(rstan)
library(tidybayes)
library(scales)
library(sp)
library(sf)
library(loo)



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



load(paste0("output/",group_name,"_",nspecies,"_hier_cov_adj_output.RData"))
load(paste0("output/",group_name,"_",nspecies,"_naive_output.RData"))



A_cov_samples = gather_draws(stanfit_hier,A_cov)
B_cov_samples = gather_draws(stanfit_hier,B_cov)

a_cov_samples = gather_draws(stanfit_hier,a_cov[sp]) %>% 
  left_join(.,species_df,by = "sp") %>% 
  group_by(english,sp) %>% 
  summarise(cov_abund = (mean(.value)),
            cov_abund_lci = (quantile(.value,0.025)),
            cov_abund_uci = (quantile(.value,0.975)),
            .groups = "keep")

b_cov_samples = gather_draws(stanfit_hier,b_cov[sp]) %>% 
  left_join(.,species_df,by = "sp") %>% 
  group_by(english,sp) %>% 
  summarise(cov_effect_slope = (mean(.value)),
            cov_effect_slope_lci = (quantile(.value,0.05)),
            cov_effect_slope_uci = (quantile(.value,0.95)),
            .groups = "keep")

b_cov_h = B_cov_samples %>% 
  summarise(cov_effect_slope = (mean(.value)),
            cov_effect_slope_lci = (quantile(.value,0.05)),
            cov_effect_slope_uci = (quantile(.value,0.95))) %>% 
  mutate(english = "Mean all species")
b_covs = bind_rows(b_cov_h,b_cov_samples)
b_covs$english <- factor(b_covs$english,levels = c("Mean all species",species_list),ordered = TRUE)

slope_effect_plot = ggplot(data = b_covs,aes(x = english,y = cov_effect_slope))+
  geom_pointrange(aes(ymin = cov_effect_slope_lci,ymax = cov_effect_slope_uci),alpha = 0.2)+
  geom_abline(slope = 0, intercept = 0)+
  labs(x = "",
       y = "Effect of Increasing Human Footprint on Trend",
       title = "Covariate Effects on Route-level Population trends",
       subtitle = paste("Species of",group_name))+
  theme_classic()+
coord_flip()
print(slope_effect_plot)

pdf(paste0("Figures/Slope_effect_",group_name,".pdf"),
    width = 11,
    height = 8.5)
print(slope_effect_plot)
dev.off()


# model fit comparison ----------------------------------------------------

cov_loo <- try(loo(stanfit_hier),silent = TRUE)


naive_loo <- try(loo(stanfit_naive),silent = TRUE)

if(class(age_loo) != "try-error"){
age_loo_out[[group_name]] <- age_loo
naive_loo_out[[group_name]] <- naive_loo
}
## looic values are interpreted similarly to AIC values (smaller = better)
## it's an estimate of the error in predicting from a new dataset (smaller = less predicted error)
## looic also comes with an estimate of Standard Error
## generally, the age models have lower looic values, but the difference is very small
## relative to the Standard Error of the looic summaries
## therefore, it will be useful to explore the looic values more closely

# extract the point wise looic, merge with data ---------------------------
## every data-point has an associated looic value, so we can explore how each data point
## contributes to the estimated prediction error


# Which observations have the bad and very bad pareto-k-diagnostics --------
## which years, observer ages, average counts, species, regions...???
## high pareto-k-diagnostics indicate that the looic values are not as reliable as others
## the only proper way to fix this is to run a full cross-validation, but our lives are too
## short for that. So we'll have to ignore the diagnostic warnings overall, but it will 
## still be sueful to explore the conditions under which these diagnostics
## suggest we should be cautious. 


# explore how looic varies among species ----------------------------------
# are some species better predicted than others


# explore how the looic varies with year ----------------------------------
## the most important issue with fit are the first and last years
## the fit at the start and end of the time-series will determine
## how well the trend is estimated. The sume of the looic values in the fist and last years
## should give a particularly useful assessment of hte difference in fit of the two models.


# next generate plots of trends -------------------------------------------

# plot the hyperparameters from the two models
BETA_samples = gather_draws(stanfit_hier,BETA[s])
BETA_a = BETA_samples %>% group_by(s) %>% 
  summarise(Trend = 100*(exp(mean(.value))-1),
            lci_Trend = 100*(exp(quantile(.value,0.025))-1),
            uci_Trend = 100*(exp(quantile(.value,0.975))-1),
            .groups = "keep")
BETA_a$species = species_list[BETA_a$s] 
BETA_a$model = "Covariate_adjusted" 


BETA_samples_n = gather_draws(stanfit_naive,BETA[s])
BETA_n = BETA_samples_n %>% group_by(s) %>% 
  summarise(Trend = 100*(exp(mean(.value))-1),
            lci_Trend = 100*(exp(quantile(.value,0.025))-1),
            uci_Trend = 100*(exp(quantile(.value,0.975))-1),
            .groups = "keep")
BETA_n$species = species_list[BETA_n$s] 
BETA_n$model = "Naive" 



BETAs = bind_rows(BETA_a,BETA_n)

B_plot = ggplot(data = BETAs,aes(x = species,y = Trend,colour = model,group = model))+
  geom_pointrange(aes(ymin = lci_Trend, ymax = uci_Trend),position = position_dodge(width = 0.25))+
  coord_flip()

pdf(paste0("figures/",group_name,"_",nspecies,"_comparing_hyperparamter_trends.pdf"))
print(B_plot)
dev.off()


# compare the trends at each route ----------------------------------------


laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system

locat = system.file("maps",
                    package = "bbsBayes")
map.file = "BBS_USGS_strata"

strata_map = read_sf(dsn = locat,
                     layer = map.file)
strata_map = st_transform(strata_map,crs = laea)

strata_map = strata_map %>% filter(COUNTRY == "CA")


####
# add trend and abundance ----------------------------------------


sumr <- summary(stanfit_hier)$summary
slopes <- NULL
sumr_n <- summary(stanfit_naive)$summary
slopes_n <- NULL
for(i in 1:nspecies){
  sps = species_list[i]
  nr = stan_data[[paste0("nroutes",i)]]
sltmp <- as.data.frame(sumr[paste0("beta",i,"[",1:nr,"]"),])
names(sltmp) <- c("b","se_mean","sd","lci","lqrt","med","uqrt","uci","b_n_eff","b_Rhat")
sltmp <- sltmp[,c("b","lci","uci","b_n_eff","b_Rhat")]
for(j in 1:3){
  sltmp[,j] <- (exp(sltmp[,j])-1)*100
}
sltmp$r <- 1:nr

alphas <- as.data.frame(sumr[paste0("alpha",i,"[",1:nr,"]"),])
names(alphas) <- c("int","se_mean","int_sd","int_lci","lqrt","med","uqrt","int_uci","int_n_eff_b","int_Rhat_b")
alphas <- alphas[,c("int","int_lci","int_uci","int_n_eff_b","int_Rhat_b")]
for(j in 1:3){
  alphas[,j] <- exp(alphas[,j])
}
alphas$r <- 1:nr

sltmp <- full_join(sltmp,alphas)
sltmp$species <- sps

slopes <- bind_rows(slopes,sltmp)

rm("sltmp")

# same for naive models ---------------------------------------------------



sltmp <- as.data.frame(sumr_n[paste0("beta",i,"[",1:nr,"]"),])
names(sltmp) <- c("b","se_mean","sd","lci","lqrt","med","uqrt","uci","b_n_eff","b_Rhat")
sltmp <- sltmp[,c("b","lci","uci","b_n_eff","b_Rhat")]
for(j in 1:3){
  sltmp[,j] <- (exp(sltmp[,j])-1)*100
}
sltmp$r <- 1:nr

alphas <- as.data.frame(sumr_n[paste0("alpha",i,"[",1:nr,"]"),])
names(alphas) <- c("int","se_mean","int_sd","int_lci","lqrt","med","uqrt","int_uci","int_n_eff_b","int_Rhat_b")
alphas <- alphas[,c("int","int_lci","int_uci","int_n_eff_b","int_Rhat_b")]
for(j in 1:3){
  alphas[,j] <- exp(alphas[,j])
}
alphas$r <- 1:nr

sltmp <- full_join(sltmp,alphas)
sltmp$species <- sps

slopes_n <- bind_rows(slopes_n,sltmp)


}

# lists of the species and model specific maps ----------------------------


maps_age = vector(mode = "list",length = length(species_list))
maps_n = vector(mode = "list",length = length(species_list))
maps_dif = vector(mode = "list",length = length(species_list))


for(sps in 1:nspecies){
  spsn = species_list[sps]
  
  load(paste0("data/speciesRouteData/",group_name,"_",spsn,"_",nspecies,"_route_data.RData"))
  
  
 tmp <- slopes %>% filter(species == spsn)
  maps_age[[sps]] = inner_join(route_map,tmp,by = c("routeF" = "r"))
  
  tcomb_age <- tmp %>% select(b,lci,uci,r,species) %>% 
    rename(b_age = b,
           lci_age = lci,
           uci_age = uci)
  
  tmp <- slopes_n %>% filter(species == spsn)
  maps_n[[sps]] = inner_join(route_map,tmp,by = c("routeF" = "r"))
  
  tcomb <- tmp %>% select(b,lci,uci,r,species,int) %>% 
    rename(b_naive = b,
           lci_naive = lci,
           uci_naive = uci) %>% 
    full_join(.,tcomb_age) %>% 
    mutate(dif = b_naive - b_age)
  maps_dif[[sps]] = inner_join(route_map,tcomb,by = c("routeF" = "r"))
  
  
}


# Explore the differences in trends and the observers age --------
## are the route-level differences in trends strongly related to the age of the observer
## on a given route? they should be if the model is doing what we think it is, but you'll 
## want to confirm that. Try plotting the difference against the average observer age on the route







# add mapping of trends ---------------------------------------------------

# plot the route-level trends maps compare adjacent maps


for(sps in 1:nspecies){
  spsn = species_list[sps]
  
  breaks <- c(-0.07, -0.04, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.04, 0.07)
  breaks <- round((exp(breaks)-1)*100,1)
  labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
  labls = paste0(labls, " %/year")
  
  route_map_age <- maps_age[[sps]]
  route_map_naive <- maps_n[[sps]]

  route_map_age$Tplot <- cut(route_map_age$b,breaks = c(-Inf, breaks, Inf),labels = labls)
  route_map_naive$Tplot <- cut(route_map_naive$b,breaks = c(-Inf, breaks, Inf),labels = labls)
  map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                   "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
  names(map_palette) <- labls
  
  
  tmap_age = ggplot(route_map_age)+
    geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
    geom_sf(aes(colour = Tplot,size = int))+
    scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                        guide = guide_legend(reverse=TRUE),
                        name = paste0("slope\n",firstYear,"-",2019))+
    labs(title = paste(spsn,"Trends after removing covariate by route (size = abundance)"))
  
  tmap_naive = ggplot(route_map_naive)+
    geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
    geom_sf(aes(colour = Tplot,size = int))+
    scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                        guide = guide_legend(reverse=TRUE),
                        name = paste0("slope\n",firstYear,"-",2019))+
    labs(title = paste(spsn,"Naive trends by route (size = abundance)"))
  
 
   breaks <- c(-0.02, -0.01, -0.005,-0.0025,-0.001,0.001, 0.0025, 0.005, 0.01, 0.02)
  breaks <- round((exp(breaks)-1)*100,1)
  labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
  labls = paste0(labls, " %/year")
  
  route_map_dif <- maps_dif[[sps]]
  route_map_dif$Tplot <- cut(route_map_dif$dif,breaks = c(-Inf, breaks, Inf),labels = labls)
  map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                   "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
  names(map_palette) <- labls
  
  tmap_dif = ggplot(route_map_dif)+
    geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
    geom_sf(aes(colour = Tplot,size = int))+
    scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                        guide = guide_legend(reverse=TRUE),
                        name = paste0("slope\n",firstYear,"-",2019))+
    labs(title = paste(spsn,"Effect of covariate on trend (, size = abundance)"))
  
  
  
  
  pdf(file = paste0("figures/",spsn,"_",group_name,"_",nspecies,"_","trend_maps_age_adjusted.pdf"),
      width = 11,
      height = 8.5)
  
  print(tmap_age)
  print(tmap_naive)
  print(tmap_dif)
  dev.off()
  
  
}




}

save(list = c("age_loo_out","naive_loo_out"),
     file = "loo_output.RData")
















