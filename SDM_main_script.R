#############################################################################
# running analysis at a resolution of 250 meters
#0. settings wd, packages:
#setwd("Z:/jseuren/Stageonderzoek")
setwd("/vol/milkunstage/jseuren/Stageonderzoek")
source("Scripts/hulpscript_250m.R")
#.libPaths("/vol/milkunM/DepLibrary/Cluster_information/packages")
# packages:
library(tidyverse)
library(dplyr)
library(sf)
library(wdpar)
library(raster)
library(biomod2)
library(terra) 
library(tidyterra)
library(ggplot2)
library(snowfall)
library(maps)
library(caret)
library(foreach)
library(confintr)
library(leaflet)
#1a. defining parameters: building the model
switch_env <- T  # do we use continuous variables? 
switch_cat <- F  # do we use categorical variables?
pseudo_strategie <- 'random' # Method for pseudo-absence generation
pseudos_rep <- 5 # PA repetitions
pseudos_n <- 20000 # PAs per rep
model_algoritmes <- c('GLM', 'GAM', 'RF', 'GBM', 'CTA')# Algorithms used
eval_metrics <- c('ROC', 'TSS')# defining evaluation metrics
cv_strategie <- 'kfold' # strategy for cross-validation
cv_rep <- 8 # number of cv runs
cv_k <- 4 # number of subsets in which original set of data is split
cv_perc <- 0.70 # proportion data used for calibrating model
var_import <- 3 # nr of mutations to variables to assess their importance
# settings for analyses:
split_evalbox <- 'algo' # based on what models are split in evaluation plot
split_varbox <- c('algo', 'expl.var', 'PA')# bases on which runs are split in variable importance plot
ens_split_varbox <- c('algo', 'expl.var', 'merged.by.PA') # same but for ensemble model
fixed_var_rc <- 'mean' # level at which other variables are kept while generating response curves
folder_modelprestatie <- "./output/modelprestatie/250m" # folder for model evalution figures
folder_figuren <- "./output/figuren" # folder for general figures.
folder_emodelprestatie <- "./output/ens_mod_pres/250m" # same for ensemble models
folder_ensemble_figuren <- "./output/ensemble_figuren" # same for ensemble models
model_id <- "250_meter" # id-tag of modelling run
###############################################################################################
#1c.parameters: ensemble models
em_choice <- 'all' # all models can theoretically be used for the ensemble
em_by <- 'all' # based on which ensemble models are made
em_models <- 'EMwmean' # criteria for weighting model
em_selectie <- 'ROC' # metric used to select models
em_thresh <- 0.8 # selection threshold
em_eval_metrics <- c('ROC', 'TSS') # metrics to evaluate ensemble models
em_var_import <- 3 # nr of mutations to variables to assess their importance.
em_conf_sig <- 0.05 # significance level for CI
em_score_decay <- 'proportional' # how weighting based on performance scores chages relative to performance scores
################################################################################################
# 1d. Settings for projections
projectie_naam <- 'projectie_NL_250' # name for projections
modelkeuze <- 'all' # which models included for projections
binary_metric <- 'all' # metric to use for binarisation
filter_metric <- 'all' # metric to use for filtering
# 1e. plots
binair_methode <- 'TSS' # for binarisation predictions.
# load in relevant files for plots:
NL_shape <- st_read("./Input_data/env_data/nl_1km.shp")
NL_shape <- st_transform(NL_shape, 4326)
NL_shape # shapefile of the netherlands
Nederland <- map_data("world", region = "Netherlands")
volgorde_soorten <- c("adder","gladde slang", "hazelworm", "levendbarende hagedis",
                      "muurhagedis", "ringslang", "zandhagedis") # order of species for plots
gemeente_selectie <- c("Ede", "Epe", "Deurne", "Huizen", "Horst aan de Maas", "Roermond","Ommen", "Westerveld", "Zutphen", "Wormerland") # selection of case study areas
# loading in occurrence and environmental data:
# 2. Input data: load in NDFF (occurrence) data en env variables
NDFF_data <- inlezen_NDFF()
NDFF_per_soort <- doSplit(NDFF_data, "soort")
# environmental data:
Env_data <- inlezen_env_data(switch_env, switch_cat)
correlatie_matrix <- correlatie_check() # check for correlation
#switch_cat <- T
#Env_data_num <- inlezen_kaarten_lr_num(switch_env, switch_cat) # deze gebruiken voor uitknippen

# now to load in shapefiles of case study areas and use it to cut them out of initial data:
Env_data_num <- Env_data
shape_case_gemeentes <- st_read("./Input_data/shapefiles/case_gemeentes.shp")
shape_case_gemeentes <- st_transform(shape_case_gemeentes, crs = st_crs(Env_data)) 
shape_Ede <- st_read("./Input_data/shapefiles/ede_shapefile.shp")
shape_Epe <- st_read("./Input_data/shapefiles/epe_shapefile.shp")
shape_Ommen <- st_read("./Input_data/shapefiles/ommen_shapefile.shp")
shape_Deurne <- st_read("./Input_data/shapefiles/deurne_shapefile.shp")
shape_Roermond <- st_read("./Input_data/shapefiles/roermond_shapefile.shp")
shape_Westerveld <- st_read("./Input_data/shapefiles/westerveld_shapefile.shp")
shape_Huizen <- st_read("./Input_data/shapefiles/huizen_shapefile.shp")
shape_Zutphen <- st_read("./Input_data/shapefiles/zutphen_shapefile.shp")
shape_Wormerland <- st_read("./Input_data/shapefiles/wormerland_shapefile.shp")
shape_Horst <- st_read("./Input_data/shapefiles/horst_shapefile.shp")

#Env_data_case <- maskeren_case_gemeentes()
#Env_data_case <- stack(Env_data_case)
#writeRaster(Env_data_case, filename = "./Input_data/env_data/stack.tif", options = "INTERLEAVE = BAND", overwrite = TRUE)
Env_data_case <- stack("./Input_data/env_data/stack.tif")
#plot(Env_data_case) . 
# extact environmental data per municipality:
Env_data_basis <- as.list(Env_data)
Env_data_Deurne <- stack(lapply(Env_data_basis, uitknippen_deurne))
Env_data_Ede <- stack(lapply(Env_data_basis, uitknippen_ede))
Env_data_Epe <- stack(lapply(Env_data_basis, uitknippen_epe))
Env_data_Horst <- stack(lapply(Env_data_basis, uitknippen_horst))
Env_data_Huizen <- stack(lapply(Env_data_basis, uitknippen_huizen))
Env_data_Ommen <- stack(lapply(Env_data_basis, uitknippen_ommen))
Env_data_Roermond <- stack(lapply(Env_data_basis, uitknippen_roermond))
Env_data_Westerveld <- stack(lapply(Env_data_basis, uitknippen_westerveld))
Env_data_Wormerland <- stack(lapply(Env_data_basis, uitknippen_wormerland))
Env_data_Zutphen <- stack(lapply(Env_data_basis, uitknippen_zutphen))
####################################################################################
# Similar operation for occurrence data: per species
NDFF_ede <- filter(NDFF_data, gemeente == "Ede")
NDFF_ede_per_soort <- doSplit(NDFF_ede, "soort") 
NDFF_epe <- filter(NDFF_data, gemeente == "Epe")
NDFF_epe_per_soort <- doSplit(NDFF_epe, "soort") 
NDFF_ommen <- filter(NDFF_data, gemeente == "Ommen")
NDFF_ommen_per_soort <- doSplit(NDFF_ommen, "soort")
NDFF_deurne <- filter(NDFF_data, gemeente == "Deurne")
NDFF_deurne_per_soort <- doSplit(NDFF_deurne, "soort")
NDFF_roermond <- filter(NDFF_data, gemeente == "Roermond")
NDFF_roermond_per_soort <- doSplit(NDFF_roermond, "soort")
NDFF_westerveld <- filter(NDFF_data, gemeente == "Westerveld")
NDFF_westerveld_per_soort <- doSplit(NDFF_westerveld, "soort")
NDFF_huizen <- filter(NDFF_data, gemeente == "Huizen")
NDFF_huizen_per_soort <- doSplit(NDFF_huizen, "soort")
NDFF_zutphen <- filter(NDFF_data, gemeente == "Zutphen")
NDFF_zutphen_per_soort <- doSplit(NDFF_zutphen, "soort")
NDFF_wormerland <- filter(NDFF_data, gemeente == "Wormerland")
NDFF_wormerland_per_soort <- doSplit(NDFF_wormerland, "soort")
NDFF_horst <- filter(NDFF_data, gemeente == "Horst aan de Maas")
NDFF_horst_per_soort <- doSplit(NDFF_horst, "soort")
# remaining data for entirety of NL, minus case study areas:
NDFF_case <- filteren_NDFF()
NDFF_case_per_soort <- doSplit(NDFF_case, "soort") # per species
#####################################################################.
# now to generate and save PAs:
Data_heel_NL <- lapply(doSplit(NDFF_data, 'soort'), data_formatteren_NL_case) 
# PAs can be extracted from this formatted data, after which they are saved:
modelData_f = Data_heel_NL 
  PA_adder_f <- cbind(modelData_f$adder@PA.table, modelData_f$adder@coord)
  PA_gladde_slang_f <- cbind(modelData_f$`gladde slang`@PA.table, modelData_f$`gladde slang`@coord)
  PA_hazelworm_f <- cbind(modelData_f$hazelworm@PA.table, modelData_f$hazelworm@coord)
  PA_levendbarende_hagedis_f <- cbind(modelData_f$`levendbarende hagedis`@PA.table, modelData_f$`levendbarende hagedis`@coord)
  PA_muurhagedis_f <- cbind(modelData_f$muurhagedis@PA.table, modelData_f$muurhagedis@coord)
  PA_ringslang_f <- cbind(modelData_f$ringslang@PA.table, modelData_f$ringslang@coord)
  PA_zandhagedis_f <- cbind(modelData_f$zandhagedis@PA.table, modelData_f$zandhagedis@coord)
  PA_f <- list(adder = PA_adder_f, gladde_slang = PA_gladde_slang_f, hazelworm = PA_hazelworm_f, levendbarende_hagedis = PA_levendbarende_hagedis_f,
               muurhagedis = PA_muurhagedis_f, ringslang = PA_ringslang_f, zandhagedis = PA_zandhagedis_f)
PA_per_soort <- PA_f # per species
lapply(names(PA_per_soort), opslaan_PA_data) # save
# load in csv to test:
PA_per_soort <- ophalen_PAs()

# now to actually produce the SDMs per species.
SDMs_250m <- lapply(doSplit(NDFF_case, "soort"), SDM_pipeline_250)
SDMs_250m <- ophalen_modellen_250()
# Analyse the models:
evaluatie_250 <- lapply(SDMs_250m, evaluatie_250) # their performance
# save the performance
lapply(names(evaluatie_250), opslaan_evaluatie)
# feature importances:
var_belang_250 <- lapply(SDMs_250m, var_belang_250) 
# save to csv:
lapply(names(var_belang_250), opslaan_belang)
# visualise performacne per algorithm:
box_prestatie_250 <- lapply(SDMs_250m, prestatie_box_250) 
#save as pngs:
lapply(names(box_prestatie_250), opslaan_box_eval_mean)
# feature importance, per algorithm:
var_box_250 <- lapply(SDMs_250m, var_imp_box_250)
#save as pngs:
lapply(names(var_box_250), opslaan_var_imp_plot)
# Finally compose and save response curves:
respons_curves_250 <- lapply(SDMs_250m, resp_curves_250)
lapply(names(respons_curves_250), opslaan_resp_curves)
#######################################################################################
# compute ensemble models per species:
ensemble_NL_250m <- lapply(SDMs_250m, ensemble_pipeline_250)
ensemble_NL_250m <- ophalen_ensemble_modellen_250()
# evaluation and analysis:
ensemble_evaluatie_data <- lapply(ensemble_NL_250m,ensemble_model_evaluatie_250) 

# save to csvs
lapply(names(ensemble_evaluatie_data), opslaan_ensemble_evaluatie) 

# feature importance:
ens_var_belang_lijst <- lapply(ensemble_NL_250m, ens_var_belang_250m)

# save
lapply(names(ens_var_belang_lijst), ens_opslaan_belang)

# visualisation steps, incl. saving figures;
ens_box_prestatie_lijst <- lapply(ensemble_NL_250m, ens_prestatie_box_250)

lapply(names(ens_box_prestatie_lijst), opslaan_ens_box_eval_mean)

ens_var_imp_plot_lijst <- lapply(ensemble_NL_250m, ens_var_imp_box_250)
lapply(names(ens_var_imp_plot_lijst), opslaan_ens_var_imp_plot)

# finally response curves:
ensemble_respons_curves_lijst <- lapply(ensemble_NL_250m, ens_resp_curves_250)
# ensemble_rc_smooth <- lapply(ensemble_NL_250m, ens_resp_curves_250_smooth)
# save
lapply(names(ensemble_respons_curves_lijst), opslaan_ens_resp_curves)
###################################################################################
# Now to project the ensemble models, first for the entire coutnry, afterwards per case study municipality:
Env_data <- inlezen_env_data(switch_env, switch_cat) # data NL

projecties_NL_250 <- lapply(ensemble_NL_250m, ensemble_projecteren_250) # projections per species
projecties_NL_250 <- ophalen_250_projecties_NL()
projecties_NL_250 <- ophalen_250_projecties_NL()
# visualisation:
kaarten_NL_250 <- lapply(projecties_NL_250, projectie_kaart_250)
hw_projectie <- BIOMOD_EnsembleForecasting(bm.em = ensemble_NL_250m$hazelworm,
                                           proj.name = projectie_naam,
                                           new.env = Env_data,           
                                           models.chosen = modelkeuze,
                                           metric.binary =
                                             binary_metric,
                                           metric.filter =
                                             filter_metric)
plot_hw <- plot(hw_projectie, plot.output = 'facet')
lh_projectie <- BIOMOD_EnsembleForecasting(bm.em = ensemble_NL_250m$levendbarende_hagedis,
proj.name = projectie_naam,
                                           new.env = Env_data,        n
                                           models.chosen = modelkeuze,
                                           metric.binary 
                                           binary_metric,
                                           metric.filter =
                                             filter_metric)
plot_lh <- plot(lh_projectie, plot.output = 'facet')
plot_adder <- plot(projecties_NL_250$adder, plot.output = 'facet')
plot_gs <- plot(projecties_NL_250$gladde_slang, plot.output = 'facet')
plot_mh <- plot(projecties_NL_250$muurhagedis, plot.output = 'facet')
plot_rs <- plot(projecties_NL_250$ringslang, plot.output = 'facet')
plot_zh <- plot(projecties_NL_250$zandhagedis, plot.output = 'facet')
kaarten_NL_250 <- list(adder = plot_adder, gladde_slang = plot_gs, hazelworm = plot_hw,
                      levendbarende_hagedis = plot_lh, muurhagedis = plot_mh, ringslang = plot_rs,
                      zandhagedis = plot_zh)
projecties_NL_250 <- list(adder = projecties_NL_250$adder, gladde_slang = projecties_NL_250$gladde_slang,
                      hazelworm = hw_projectie, levendbarende_hagedis = lh_projectie,
                      muurhagedis = projecties_NL_250$muurhagedis, ringslang = projecties_NL_250$ringslang,
                      zandhagedis = projecties_NL_250$zandhagedis)
# save these maps:
lapply(names(kaarten_NL_250), kaarten_opslaan_250)

# binary range per species:
range_NL_250 <- lapply(projecties_NL_250, range_binair_250)
range_plots_NL_250 <- lapply(range_NL_250, plot_range_250)
# save range plots: UNFORTUNATELY DOES NOT WORK
lapply(names(range_plots_NL_250), opslaan_range_plot_250) 
# interactive map with all projections:
geschiktheid_NL_per_soort()
#####################################################################################
# Now to project, visualise and save for every case study area.
# Deurne
projectie_Deurne <- lapply(ensemble_NL_250m, projectie_Deurne)
projectie_Deurne <- ophalen_projecties_Deurne()
kaarten_Deurne <- lapply(projectie_Deurne, proj_kaart_Deurne)
lapply(names(kaarten_Deurne), kaarten_Deurne_opslaan)
ranges_Deurne <- lapply(projectie_Deurne, range_binair_Deurne)
r_plots_Deurne <- lapply(ranges_Deurne, range_plot_Deurne)
geschiktheid_Deurne_per_soort()
# Ede
projectie_Ede <- lapply(ensemble_NL_250m, projectie_Ede)
projectie_Ede <- ophalen_projecties_Ede()
kaarten_Ede <- lapply(projectie_Ede, proj_kaart_Ede)
lapply(names(kaarten_Ede), kaarten_Ede_opslaan)
ranges_Ede <- lapply(projectie_Ede, range_binair_Ede)
r_plots_Ede <- lapply(ranges_Ede, range_plot_Ede)
geschiktheid_Ede_per_soort()

# Epe
projectie_Epe <- lapply(ensemble_NL_250m, projectie_Epe)
projectie_Epe <- ophalen_projecties_Epe()
kaarten_Epe <- lapply(projectie_Epe, proj_kaart_Epe)
lapply(names(kaarten_Epe), kaarten_Epe_opslaan)
ranges_Epe <- lapply(projectie_Epe, range_binair_Epe)
r_plots_Epe <- lapply(ranges_Epe, range_plot_Epe)
geschiktheid_Epe_per_soort()

# Horst
projectie_Horst <- lapply(ensemble_NL_250m, projectie_Horst)
projectie_Horst <- ophalen_projecties_Horst()
kaarten_Horst <- lapply(projectie_Horst, proj_kaart_Horst)
lapply(names(kaarten_Horst), kaarten_Horst_opslaan)
ranges_Horst <- lapply(projectie_Horst, range_binair_Horst)
r_plots_Horst <- lapply(ranges_Horst, range_plot_Horst)
geschiktheid_Horst_per_soort()

# Huizen
projectie_Huizen <- lapply(ensemble_NL_250m, projectie_Huizen)
projectie_Huizen <- ophalen_projecties_Huizen()
kaarten_Huizen <- lapply(projectie_Huizen, proj_kaart_Huizen)
lapply(names(kaarten_Huizen), kaarten_Huizen_opslaan)
ranges_Huizen <- lapply(projectie_Huizen, range_binair_Huizen)
r_plots_Huizen <- lapply(ranges_Huizen, range_plot_Huizen)
geschiktheid_Huizen_per_soort()

# Ommen
projectie_Ommen <- lapply(ensemble_NL_250m, projectie_Ommen)
projectie_Ommen <- ophalen_projecties_Ommen()
kaarten_Ommen <- lapply(projectie_Ommen, proj_kaart_Ommen)
lapply(names(kaarten_Ommen), kaarten_Ommen_opslaan)
# binaire range ophalen:
ranges_Ommen <- lapply(projectie_Ommen, range_binair_Ommen)
# binaire kaarten produceren:
r_plots_Ommen <- lapply(ranges_Ommen, range_plot_Ommen)
geschiktheid_Ommen_per_soort()

# Roermond
projectie_Roermond <- lapply(ensemble_NL_250m, projectie_Roermond)
projectie_Roermond <- ophalen_projecties_Roermond()
kaarten_Roermond <- lapply(projectie_Roermond, proj_kaart_Roermond)
lapply(names(kaarten_Roermond), kaarten_Roermond_opslaan)
ranges_Roermond <- lapply(projectie_Roermond, range_binair_Roermond)
r_plots_Roermond <- lapply(ranges_Roermond, range_plot_Roermond)
geschiktheid_Roermond_per_soort()

# Westerveld
projectie_Westerveld <- lapply(ensemble_NL_250m, projectie_Westerveld)
projectie_Westerveld <- ophalen_projecties_Westerveld()
kaarten_Westerveld <- lapply(projectie_Westerveld, proj_kaart_Westerveld)
lapply(names(kaarten_Westerveld), kaarten_Westerveld_opslaan)
ranges_Westerveld <- lapply(projectie_Westerveld, range_binair_Westerveld)
r_plots_Westerveld<- lapply(ranges_Westerveld, range_plot_Westerveld)
geschiktheid_Westerveld_per_soort()

# Wormerland
projectie_Wormerland <- lapply(ensemble_NL_250m, projectie_Wormerland)
projectie_Wormerland <- ophalen_projecties_Wormerland()
kaarten_Wormerland <- lapply(projectie_Wormerland, proj_kaart_Wormerland)
lapply(names(kaarten_Wormerland), kaarten_Wormerland_opslaan)
ranges_Wormerland <- lapply(projectie_Wormerland, range_binair_Wormerland)
r_plots_Wormerland <- lapply(ranges_Wormerland, range_plot_Wormerland)
geschiktheid_Wormerland_per_soort()

# Zutphen
projectie_Zutphen <- lapply(ensemble_NL_250m, projectie_Zutphen)
projectie_Zutphen <- ophalen_projecties_Zutphen()
kaarten_Zutphen <- lapply(projectie_Zutphen, proj_kaart_Zutphen)
lapply(names(kaarten_Zutphen), kaarten_Zutphen_opslaan)
ranges_Zutphen <- lapply(projectie_Zutphen, range_binair_Zutphen)
r_plots_Zutphen <- lapply(ranges_Zutphen, range_plot_Zutphen)
geschiktheid_Zutphen_per_soort()
###############################################################################
# Now to compare the predictions to the real occurrences to assess their performance.
# For this, we must rasterise predictions, and also make a raster combining PAs and occurrences.
# PAs to raster:
PA_ede_raster <- lapply(PA_per_soort, uitknippen_PA_ede)
PA_epe_raster <- lapply(PA_per_soort, uitknippen_PA_epe)
PA_deurne_raster <- lapply(PA_per_soort, uitknippen_PA_deurne)
PA_roermond_raster <- lapply(PA_per_soort, uitknippen_PA_roermond)
PA_ommen_raster <- lapply(PA_per_soort, uitknippen_PA_ommen)
PA_huizen_raster <- lapply(PA_per_soort, uitknippen_PA_huizen)
PA_westerveld_raster <- lapply(PA_per_soort, uitknippen_PA_westerveld)
PA_wormerland_raster <- lapply(PA_per_soort, uitknippen_PA_wormerland)
PA_horst_raster <- lapply(PA_per_soort, uitknippen_PA_horst)
PA_zutphen_raster <- lapply(PA_per_soort, uitknippen_PA_zutphen)

# occurrences to raster:
waarnemingen_rasters_Ede <- lapply(NDFF_ede_per_soort, raster_maken_waarnemingen_ede)
waarnemingen_rasters_Epe <- lapply(NDFF_epe_per_soort, raster_maken_waarnemingen_epe)
waarnemingen_rasters_Ommen <- lapply(NDFF_ommen_per_soort, raster_maken_waarnemingen_ommen)
waarnemingen_rasters_Deurne <- lapply(NDFF_deurne_per_soort, raster_maken_waarnemingen_deurne)
waarnemingen_rasters_Roermond <- lapply(NDFF_roermond_per_soort, raster_maken_waarnemingen_roermond)
waarnemingen_rasters_Huizen <- lapply(NDFF_huizen_per_soort, raster_maken_waarnemingen_huizen)
waarnemingen_rasters_Westerveld <- lapply(NDFF_westerveld_per_soort, raster_maken_waarnemingen_westerveld)
waarnemingen_rasters_Wormerland <- lapply(NDFF_wormerland_per_soort, raster_maken_waarnemingen_wormerland)
waarnemingen_rasters_Horst <- lapply(NDFF_horst_per_soort, raster_maken_waarnemingen_horst)
waarnemingen_rasters_Zutphen <- lapply(NDFF_zutphen_per_soort, raster_maken_waarnemingen_zutphen)
##################################################################################################
# Merge rasters for PAs and occurrences and save them to dataframe:
PA_W_ede <- merge_ede()
echte_data_Ede <- df_data_ede()
PA_W_epe <- merge_epe()
echte_data_Epe <- df_data_epe()
PA_W_deurne <- merge_deurne()
echte_data_Deurne <- df_data_deurne()
PA_W_roermond <- merge_roermond()
echte_data_Roermond <- df_data_roermond()
PA_W_huizen <- merge_huizen()
echte_data_Huizen <- df_data_huizen()
PA_W_ommen <- merge_ommen()
echte_data_Ommen <- df_data_ommen()
PA_W_westerveld <- merge_westerveld()
echte_data_Westerveld <- df_data_westerveld()
PA_W_wormerland <- merge_wormerland()
echte_data_Wormerland <- df_data_wormerland()
PA_W_horst <- merge_horst()
echte_data_Horst <- df_data_horst()
PA_W_zutphen <- merge_zutphen()
echte_data_Zutphen <- df_data_zutphen()
####################################################################################
# Now to rasterise the predictions by the ensemble SDMs.
raster_voorspelling_Ede <- lapply(range_NL_250, raster_projectie_ede)
raster_voorspelling_Epe <- lapply(range_NL_250, raster_projectie_epe)
raster_voorspelling_Deurne <- lapply(range_NL_250, raster_projectie_deurne)
raster_voorspelling_Horst <- lapply(range_NL_250, raster_projectie_horst)
raster_voorspelling_Huizen <- lapply(range_NL_250, raster_projectie_huizen)
raster_voorspelling_Ommen <- lapply(range_NL_250, raster_projectie_ommen)
raster_voorspelling_Roermond <- lapply(range_NL_250, raster_projectie_roermond)
raster_voorspelling_Westerveld <- lapply(range_NL_250, raster_projectie_westerveld)
raster_voorspelling_Wormerland <- lapply(range_NL_250, raster_projectie_wormerland)
raster_voorspelling_Zutphen <- lapply(range_NL_250, raster_projectie_zutphen)
# Also extract these predictions and save in dataframe
voorspelling_Ede <- df_voorspeld_ede()
voorspelling_Epe <- df_voorspeld_epe()
voorspelling_Deurne <- df_voorspeld_deurne()
voorspelling_Horst <- df_voorspeld_horst()
voorspelling_Huizen <- df_voorspeld_huizen()
voorspelling_Ommen <- df_voorspeld_ommen()
voorspelling_Roermond <- df_voorspeld_roermond()
voorspelling_Westerveld <- df_voorspeld_westerveld()
voorspelling_Wormerland <- df_voorspeld_wormerland()
voorspelling_Zutphen <- df_voorspeld_zutphen()
#####################################################################################
# Now create a confusion matrix comparing the real data and predictions for every municipality:
# Ede:
confusion_matrices_Ede <- cm_ede()
# Epe:
confusion_matrices_Epe <- cm_epe()
# Deurne:
confusion_matrices_Deurne <- cm_deurne()
# Horst:
confusion_matrices_Horst <- cm_horst()
# Huizen:
confusion_matrices_Huizen <- cm_huizen()
# Ommen: 
confusion_matrices_Ommen <- cm_ommen()
# Roermond:
confusion_matrices_Roermond <- cm_roermond()
# Westerveld:
confusion_matrices_Westerveld <- cm_westerveld()
# Wormerland:
confusion_matrices_Wormerland <- cm_wormerland()
# Zutphen:
confusion_matrices_Zutphen <- cm_zutphen()
#############################################################################
# Extract the data from these matrices:
# Ede:
df_cm_ede <- lapply(confusion_matrices_Ede, ophalen_info_cm_ede)                        
data_cm_ede <- bind_rows(df_cm_ede, .id = "soort")
data_cm_ede <- data_cm_ede %>%
  mutate(gemeente = c("Ede", "Ede", "Ede", "Ede", "Ede", "Ede"),
         .before = soort)
# Epe:
df_cm_epe <- lapply(confusion_matrices_Epe, ophalen_info_cm_epe)                        
data_cm_epe <- bind_rows(df_cm_epe, .id = "soort")
data_cm_epe <- data_cm_epe %>%
  mutate(gemeente = c("Epe", "Epe", "Epe", "Epe", "Epe", "Epe"),
         .before = soort)
# Deurne:
df_cm_deurne <- lapply(confusion_matrices_Deurne, ophalen_info_cm_deurne)                        
data_cm_deurne <- bind_rows(df_cm_deurne, .id = "soort")
data_cm_deurne <- data_cm_deurne %>%
  mutate(gemeente = c("Deurne", "Deurne", "Deurne"),
         .before = soort)
# Horst:
df_cm_horst <- lapply(confusion_matrices_Horst, ophalen_info_cm_horst)                        
data_cm_horst <- bind_rows(df_cm_horst, .id = "soort")
data_cm_horst <- data_cm_horst %>%
  mutate(gemeente = c("Horst", "Horst"),
         .before = soort)
# Huizen:
df_cm_huizen <- lapply(confusion_matrices_Huizen, ophalen_info_cm_huizen)                        
data_cm_huizen <- bind_rows(df_cm_huizen, .id = "soort")
data_cm_huizen <- data_cm_huizen %>%
  mutate(gemeente = c("Huizen", "Huizen", "Huizen"),
         .before = soort)
# Ommen:
df_cm_ommen <- lapply(confusion_matrices_Ommen, ophalen_info_cm_ommen)                        
data_cm_ommen <- bind_rows(df_cm_ommen, .id = "soort")
data_cm_ommen <- data_cm_ommen %>%
  mutate(gemeente = c("Ommen", "Ommen", "Ommen", "Ommen", "Ommen"),
         .before = soort)
# Roermond:
df_cm_roermond <- lapply(confusion_matrices_Roermond, ophalen_info_cm_roermond)                        
data_cm_roermond <- bind_rows(df_cm_roermond, .id = "soort")
data_cm_roermond <- data_cm_roermond %>%
  mutate(gemeente = c("Roermond", "Roermond", "Roermond", "Roermond", "Roermond"),
         .before = soort)
# Westerveld:
df_cm_westerveld <- lapply(confusion_matrices_Westerveld, ophalen_info_cm_westerveld)                        
data_cm_westerveld <- bind_rows(df_cm_westerveld, .id = "soort")
data_cm_westerveld <- data_cm_westerveld %>%
  mutate(gemeente = c("Westerveld", "Westerveld", "Westerveld", "Westerveld", "Westerveld", "Westerveld", "Westerveld"),
         .before = soort)
# Wormerland:
df_cm_wormerland <- lapply(confusion_matrices_Wormerland, ophalen_info_cm_wormerland)                        
data_cm_wormerland <- bind_rows(df_cm_wormerland, .id = "soort")
data_cm_wormerland <- data_cm_wormerland %>%
  mutate(gemeente = "Wormerland",
         .before = soort)
# Zutphen:
df_cm_zutphen <- lapply(confusion_matrices_Zutphen, ophalen_info_cm_zutphen)                        
data_cm_zutphen <- bind_rows(df_cm_zutphen, .id = "soort")
data_cm_zutphen <- data_cm_zutphen %>%
  mutate(gemeente = c("Zutphen", "Zutphen", "Zutphen"),
         .before = soort)
# Combine in one big dataframe:
data_cm_gecombineerd <- bind_rows(data_cm_deurne, data_cm_ede, data_cm_epe, data_cm_horst, data_cm_huizen, data_cm_ommen, data_cm_roermond,
data_cm_westerveld, data_cm_wormerland, data_cm_zutphen)
# Add column with calculated TSS:
data_cm_gecombineerd <- data_cm_gecombineerd %>%
  mutate(TSS = data_cm_gecombineerd$Sensitivity + data_cm_gecombineerd$Specificity - 1,
         .after = Specificity)
# save for further analysis:
write.csv(data_cm_gecombineerd, file = "./output/modelprestatie/250m/data_confusion_matrix_250.csv",
          row.names = F)