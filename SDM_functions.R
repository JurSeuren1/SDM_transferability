# In this script the functions used in the main script are defined
# working directory:.
setwd("/vol/milkunstage/jseuren/Stageonderzoek")
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
# generic function to split data
doSplit <- function(x,y=cols){
  if(length(y)==1){
    x_split <- split(x,list(x[,y]),drop = T)
  }else{
    x_split <- split(x,as.list(x[,y]),drop = T)
  }
}
# function to load in occurrence data:
inlezen_NDFF <- function(loc="./Input_data/NDFF_soortdata.csv")
{
  df <- read.csv("./Input_data/NDFF_data_incl_gemeente.csv", sep = ";") 
  colnames(df) <- c("soort", "aanwezigheid", "X", "Y", "gemeente") 
  df[df == "Adder"] <- "adder"
  df[df == "Gladde slang"] <- "gladde slang"
  df[df == "Hazelworm"] <- "hazelworm"
  df[df == "Levendbarende hagedis"] <- "levendbarende hagedis"
  df[df == "Muurhagedis"] <- "muurhagedis"
  df[df == "Ringslang"] <- "ringslang"
  df[df == "Zandhagedis"] <- "zandhagedis"
  return(df)
}

# loading in environmental data
inlezen_env_data <- function(s_klim_f = switch_env, s_cat_f = switch_cat,
                             NL_f = st_read("./Input_data/env_data/nl_1km.shp")){
  leeg_f <- raster(xmn = 3.083333, xmx = 7.225, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  NL_f <- st_transform(NL_f, 4326)
  if(s_klim_f){ # continuous variables:
    env_num_stack_f <- stack(c(gem_temp = raster("./Input_data/env_data/env_bio_1.tif"), 
                               temp_sg_var = raster("./Input_data/env_data/env_bio_4.tif"),
                               neerslag_jaarlijks = raster("./Input_data/env_data/env_bio_12.tif"),
                               neerslag_sg_var = raster("./Input_data/env_data/env_bio_15.tif")))
    env_num_crop_f <- crop(env_num_stack_f, NL_f)# voor NL
    env_num_crop_f <- resample(env_num_crop_f, leeg_f, method = "bilinear")# voor NL
    # also load in proportional categorical data: this conversion works for now!
    lu_f <- rast("./Input_data/env_data/land_use_proportions.tif")
    lu_s_f <- stack(lu_f)
    lu_NL <- crop(lu_s_f, NL_f)
    lu_s_NL <- resample(lu_NL, leeg_f, method = 'bilinear')
    st_f <- rast("./Input_data/env_data/LBK_general.tif")
    st_s_f <- stack(st_f)
    st_NL <- crop(st_s_f, NL_f)
    st_s_NL <- resample(st_NL, leeg_f, method = 'bilinear')
    cat_prop_stack_f <- stack(c(land_use = lu_s_NL,
                                soi_type = st_s_NL))
    # crop to NL:
    cat_prop_crop_f <- crop(cat_prop_stack_f, NL_f)# voor NL
    cat_prop_crop_f <- resample(cat_prop_crop_f, leeg_f, method = "bilinear")# voor NL
    crs(cat_prop_crop_f) <- crs(env_num_crop_f)
    env_stack_f <- stack(env_num_crop_f, cat_prop_crop_f)
  }
  if(s_cat_f){
    Cat_f <- raster::stack(c(LGN = as.factor(raster("./Input_data/env_data/env_LGN_hr.tif")), # categorische variabelen
                             bodemsoort = as.factor(raster("./Input_data/env_data/env_bodemsoort_hr2.tif")))) 
    LGN_coarse_f <- projectRaster(from = Cat_f$LGN, to = env_crop_f, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                                  method = 'ngb')
    LGN_c_factor_f <- as.factor(LGN_coarse_f)
    bodem_coarse_f <- projectRaster(from = Cat_f$bodemsoort, to = env_crop_f, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                                    method = 'ngb')
    bodem_c_factor_f <- as.factor(bodem_coarse_f)
    Cat_coarse_f <- stack(LGN_c_factor_f, bodem_c_factor_f)
    atr_LGN <- data.frame(
      ID = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),
      landgebruik = c("nvt","gras","landbouw","kassen","boomgaard","loofbos","naaldbos","zoet water","zout water",
                      "bebouwing","bos in bebouwd gebied","gras in bebouwd gebied","verhard","overige moerasvegetatie","open zand","duinen",
                      "heide","matig  vergraste heide","sterk vergraste heide","hoogveen","bos in hoogveengebied",
                      "rietvegetatie","bos in moerasgebied","gras (kustgebied)","struikvegetatie in hoogveengebied","struikvegetatie in moerasgebied",
                      "overige struikvegetatie"
      ))
    levels(Cat_coarse_f[[1]]) <- atr_LGN
    atr_bodem <- data.frame(ID = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),
                            bodemsoort = c("geen data", "kalkrijke zandgronden","water","kalkarme zandgronden","dijk","zeeklei", "moerige gronden","veengronden","moerasgrond","kalkrijke kleigronden","leemgronden","kalkarme kleigronden","kalkarme zavelgronden","kalkrijke zavelgronden", "lössgronden","gors","droge dalbodem"))
    levels(Cat_coarse_f[[2]]) <- atr_bodem
  }
  
  df_out <- if(s_cat_f & s_klim_f) {
    raster::stack(c(env_stack_f, Cat_coarse_f))
  } else{
    if(s_klim_f) {
      env_stack_f
    } else { Cat_coarse_f
    }
  }
  return(df_out)
}
# correlation check:
correlatie_check <- function(variabelen_f = Env_data){
  cor_raster_f <- layerStats(variabelen_f, 'pearson', na.rm=T)
  corr_matrix_f <- cor_raster_f$'pearson correlation coefficient'
  print(corr_matrix_f)
  return(corr_matrix_f)
}
# loading in environmental data numerically, works best for rasterisation:
inlezen_kaarten_lr_num <- function(s_klim_f = switch_env, s_cat_f = switch_cat,
                                   NL_f = st_read("./Input_data/env_data/nl_1km.shp")){
  leeg_f <- raster(xmn = 3.083333, xmx = 7.225, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  if(s_klim_f){
    env_stack_f <- stack(c(gem_temp = raster("./Input_data/env_data/bio_1.tif"), # continue variabelen
                           temp_sg_var = raster("./Input_data/env_data/bio_4.tif"),
                           neerslag_jaarlijks = raster("./Input_data/env_data/bio_12.tif"),
                           neerslag_sg_var = raster("./Input_data/env_data/bio_15.tif")))
    env_crop_f <- crop(env_stack_f, NL_f)# voor NL
    env_crop_f <- resample(env_crop_f, leeg_f, method = "bilinear")
  }
  if(s_cat_f){
    Cat_f <- raster::stack(c(LGN = raster("Input_data/env_data/land_use_proportions.tif"), # categorische variabelen
                             bodemsoort = raster("Input_data/env_data/LBK_general.tif")))
    Cat_f <- resample(Cat_f, leeg_f, method = "ngb")
  }
  df_out <- if(s_cat_f & s_klim_f) {
    raster::stack(c(env_crop_f, Cat_f))
  } else{
    if(s_klim_f) {
      env_crop_f
    } else { Cat_f
    }
  }
  return(df_out)
}

# masking case study areas
maskeren_case_gemeentes <- function(env_data_f = Env_data, shape_f = shape_case_gemeentes) {
  data_zonder_case_f <- mask(env_data_f, shape_f, inverse = TRUE)
  return(data_zonder_case_f)
}
# extract data per municipality:
uitknippen_deurne <- function(env_f = Env_data_basis, shape_f = shape_Deurne)
{
  data_f <- crop(env_f, shape_f)
  return(data_f)
}
uitknippen_ede <- function(env_f = Env_data_basis, shape_f = shape_Ede)
{
  data_f <- crop(env_f, shape_f)
  return(data_f)
}
uitknippen_epe <- function(env_f = Env_data_basis, shape_f = shape_Epe)
{
  data_f <- crop(env_f, shape_f)
  return(data_f)
}
uitknippen_horst <- function(env_f = Env_data_basis, shape_f = shape_Horst)
{
  data_f <- crop(env_f, shape_f)
  return(data_f)
}
uitknippen_huizen <- function(env_f = Env_data_basis, shape_f = shape_Huizen)
{
  data_f <- crop(env_f, shape_f)
  return(data_f)
}
uitknippen_ommen <- function(env_f = Env_data_basis, shape_f = shape_Ommen)
{
  data_f <- crop(env_f, shape_f)
  return(data_f)
}
uitknippen_roermond <- function(env_f = Env_data_basis, shape_f = shape_Roermond)
{
  data_f <- crop(env_f, shape_f)
  return(data_f)
}
uitknippen_westerveld <- function(env_f = Env_data_basis, shape_f = shape_Westerveld)
{
  data_f <- crop(env_f, shape_f)
  return(data_f)
}
uitknippen_wormerland <- function(env_f = Env_data_basis, shape_f = shape_Wormerland)
{
  data_f <- crop(env_f, shape_f)
  return(data_f)
}
uitknippen_zutphen <- function(env_f = Env_data_basis, shape_f = shape_Zutphen)
{
  data_f <- crop(env_f, shape_f)
  return(data_f)
}
# filtering occurrence records for whole country
filteren_NDFF <- function(ndff_f = NDFF_data, selectie_f = gemeente_selectie)
{
  ndff_filter_f <- ndff_f %>%
    filter(!(gemeente %in% selectie_f))         
  return(ndff_filter_f)
}
# Formatting data and extracting PAs
data_formatteren_NL_case <- function(ndff_f = NDFF_case, Env_data_f = Env_data, pseudo_strategie_f = pseudo_strategie,
                                     pseudos_rep_f = 1, pseudos_n_f = pseudos_n)
{
  modelData_f <- BIOMOD_FormatingData(resp.var = ndff_f$aanwezigheid, # aanwezigheidskolom
                                      expl.var = Env_data_f,              # omgevingsvariabelen
                                      resp.xy = ndff_f[, c("X", "Y")],    # coördinaten waarnemingen
                                      resp.name = unique(ndff_f$soort),   # soortnamen
                                      PA.nb.rep = pseudos_rep_f,          # ingestelde aantal PA datasets
                                      PA.nb.absences = pseudos_n_f,       # ingestelde aantal PAs per dataset
                                      PA.strategy = pseudo_strategie_f,
                                      filter.raster = TRUE)
  return(modelData_f)
}
# PA extraction:
extract_PAs <- function(modelData_f = Data_heel_NL) 
{
  PA_adder_f <- cbind(modelData_f$adder@PA.table, modelData_f$adder@coord)
  PA_gladde_slang_f <- cbind(modelData_f$`gladde slang`@PA.table, modelData_f$`gladde slang`@coord)
  PA_hazelworm_f <- cbind(modelData_f$hazelworm@PA.table, modelData_f$hazelworm@coord)
  PA_levendbarende_hagedis_f <- cbind(modelData_f$`levendbarende hagedis`@PA.table, modelData_f$`levendbarende hagedis`@coord)
  PA_muurhagedis_f <- cbind(modelData_f$muurhagedis@PA.table, modelData_f$muurhagedis@coord)
  PA_ringslang_f <- cbind(modelData_f$ringslang@PA.table, modelData_f$ringslang@coord)
  PA_zandhagedis_f <- cbind(modelData_f$zandhagedis@PA.table, modelData_f$zandhagedis@coord)
  PA_f <- list(adder = PA_adder_f, gladde_slang = PA_gladde_slang_f, hazelworm = PA_hazelworm_f, levendbarende_hagedis = PA_levendbarende_hagedis_f,
               muurhagedis = PA_muurhagedis_f, ringslang = PA_ringslang_f, zandhagedis = PA_zandhagedis_f)
  return(PA_f)
}
# save PAs:
opslaan_PA_data <- function(x){
  PA_f <- PA_per_soort[[x]]
  write.csv(PA_f, file = paste0("./output/pseudo_absences_250/Nederland/", x , ".csv"), row.names = F)
} 
ophalen_PAs <- function() {
  adder_f <- read.csv("./output/pseudo_absences_250/Nederland/adder.csv")
  gladde_slang_f <- read.csv("./output/pseudo_absences_250/Nederland/gladde_slang.csv")
  hazelworm_f <- read.csv("./output/pseudo_absences_250/Nederland/hazelworm.csv")
  levendbarende_hagedis_f <- read.csv("./output/pseudo_absences_250/Nederland/levendbarende_hagedis.csv")
  muurhagedis_f <- read.csv("./output/pseudo_absences_250/Nederland/muurhagedis.csv")
  ringslang_f <- read.csv("./output/pseudo_absences_250/Nederland/ringslang.csv")
  zandhagedis_f <- read.csv("./output/pseudo_absences_250/Nederland/zandhagedis.csv")
  PAs_f <- list(adder = adder_f, gladde_slang = gladde_slang_f, hazelworm = hazelworm_f,
                levendbarende_hagedis = levendbarende_hagedis_f, muurhagedis = muurhagedis_f,
                ringslang = ringslang_f, zandhagedis = zandhagedis_f)
  return(PAs_f)
}
# Function for individual SDMs:
SDM_pipeline_250 <- function(ndff_f = NDFF_case, Env_data_f = Env_data_case, pseudo_strategie_f = pseudo_strategie,pseudos_rep_f = pseudos_rep, pseudos_n_f = pseudos_n,model_algoritmes_f = model_algoritmes, eval_metrics_f = eval_metrics, cpus_f = cpus,cv_strat_f = cv_strategie,cv_rep_f = cv_rep, cv_k_f = cv_k, cv_perc_f = cv_perc,
                            var_import_f = var_import, model_id_f = model_id)
{ # eerst de input-data op de goede manier formatteren:
  modelData_f <- BIOMOD_FormatingData(resp.var = ndff_f$aanwezigheid, # aanwezigheidskolom
                                      expl.var = Env_data_f,              # omgevingsvariabelen
                                      resp.xy = ndff_f[, c("X", "Y")],    # coördinaten waarnemingen
                                      resp.name = unique(ndff_f$soort),   # soortnamen
                                      PA.nb.rep = pseudos_rep_f,          # ingestelde aantal PA datasets
                                      PA.nb.absences = pseudos_n_f,       # ingestelde aantal PAs per dataset
                                      PA.strategy = pseudo_strategie_f, # strategie PAs
                                      filter.raster = TRUE)   # ingestelde strategie voor genereren PAs
  # vervolgens kunnen met deze data modellen worden opgesteld: 
  #modOptions_f <- BIOMOD_ModelingOptions()
  # Los/enkel model: inclusief splitten data voor cross-validation:
  Model_soort_f <- BIOMOD_Modeling(bm.format = modelData_f, # hiervoor geformatteerde data
                                   modeling.id = model_id_f,  # eventueel mogelijkheid om een simulatie set te runnen
                                   models = model_algoritmes_f,# gewenste algoritmen, zoals gedefinieerd in hoofdscript
                                   CV.strategy = pseudo_strategie_f,   # strategie voor kruisvalidatie, ingesteld in hoofdscript
                                   CV.nb.rep = cv_rep_f,       # aantal herhalingen kruisvalidatie, ingesteld in hoofdscript 
                                   CV.k = cv_k_f,            # hoeveel delen uit oorspronkelijke dataset, ingesteld in hoofdscript (onderbouwen in verslag!)
                                   CV.perc = cv_perc,      # proportie achtergehouden voor kruisvalidatie
                                   metric.eval = eval_metrics, # welke parameters om model te beoordelen, ingesteld in hoofdscript
                                   var.import = var_import_f)  # aantal aanpassingen variabelen om belang te toetsen, ingesteld in hoofdscript
  return(Model_soort_f)
}
# retrieve models after they were saved:
ophalen_modellen_250 <- function(){
  adder_f = load("./adder/adder.250_meter.models.out")
  gladde_slang_f = load("./gladde.slang/gladde.slang.250_meter.models.out")
  hazelworm_f = load("./hazelworm/hazelworm.250_meter.models.out")
  levendbarende_hagedis_f = load("./levendbarende.hagedis/levendbarende.hagedis.250_meter.models.out")
  muurhagedis_f = load("./muurhagedis/muurhagedis.250_meter.models.out")
  ringslang_f = load("./ringslang/ringslang.250_meter.models.out")
  zandhagedis_f = load("./zandhagedis/zandhagedis.250_meter.models.out")
  modellen_f <- list(adder = adder.250_meter.models.out, gladde_slang = gladde.slang.250_meter.models.out, 
                     hazelworm = hazelworm.250_meter.models.out, levendbarende_hagedis = levendbarende.hagedis.250_meter.models.out,
                     muurhagedis = muurhagedis.250_meter.models.out, ringslang = ringslang.250_meter.models.out,
                     zandhagedis = zandhagedis.250_meter.models.out)
  return(modellen_f)
}
# model evaluation:
evaluatie_250 <- function(model_f = SDMs_250m, folder_modelprestatie_f = folder_modelprestatie) {
  model_eval_f <- data.frame(get_evaluations(model_f), check.rows = TRUE)
  return(model_eval_f)
}
# saving
opslaan_evaluatie <- function(x, data_f = evaluatie_250){
  data_f <- evaluatie_250[[x]]
  write.csv(data_f, file = paste0("./output/modelprestatie/250m/", x , ".csv"))
}
# feature importance:
var_belang_250 <- function(model_f = SDMs_250m) {
  var_belang_f <- data.frame(get_variables_importance(model_f), check.rows = TRUE)
  return(var_belang_f)
}
# and save:
opslaan_belang <- function(x, belang_f = var_belang_250){
  belang_f <- var_belang_250[[x]]
  write.csv(belang_f, file = paste0("./output/belang_variabelen/250m/", x , ".csv"))
}

# visualisation and saving this of the previous elements:
prestatie_box_250 <- function(Model_soort_f = SDMs_250m, split_evalbox_f = split_evalbox,
                              folder_figuren_f = folder_figuren) {
  plot_evalmean_f <- bm_PlotEvalMean(bm.out = Model_soort_f,
                                     group.by = split_evalbox_f)
  return(plot_evalmean_f)
}

opslaan_box_eval_mean <- function(x, eval_box_f = box_prestatie_250){
  box_ev_mean_f <- eval_box_f[[x]]$plot
  ggsave(box_ev_mean_f, device = png, file = paste0("./output/figuren/gem_prestatie/250m/", x , ".png"))
}


var_imp_box_250 <- function(Model_soort_f = SDMs_250m, split_varbox_f = split_varbox,
                            folder_figuren_f = folder_figuren) {
  plot_var_imp_f <- bm_PlotVarImpBoxplot(bm.out = Model_soort_f, 
                                         group.by = split_varbox_f)
  return(plot_var_imp_f)
}

opslaan_var_imp_plot <- function(x, var_imp_f = var_box_250){
  box_var_imp_f <- var_imp_f[[x]]$plot
  ggsave(box_var_imp_f, device = png, file = paste0("./output/figuren/var_belang/250m/", x , ".png"))
}
# response curves:
resp_curves_250 <-  function(Model_soort_f = SDMs_250m, fixed_var_rc_f = fixed_var_rc,
                             folder_figuren_f = folder_figuren, Env_data_f = Env_data_case) {
  models_f <- get_built_models(Model_soort_f, run = "RUN1")
  plot_rc_f <- bm_PlotResponseCurves(bm.out = Model_soort_f, 
                                     models.chosen = models_f,
                                     fixed.var = fixed_var_rc_f)
  #plot_rc_f <- plot_rc_f$plot +
  #             theme(legend.text = element_text(size = 0.01),
  #                   legend.key.spacing.x = unit(0.5, 'pt'),
  #                   legend.key.spacing.y = unit(0.5, 'pt'),
  #                   legend.key.height = unit(0.3, 'pt'),
  #                   legend.key.width = unit(1, 'pt'))
  # return(plot_rc_f)
}
# save:
opslaan_resp_curves <- function(x, resp_curves_f = respons_curves_250){
  respons_curves_f <- resp_curves_f[[x]]$plot
  ggsave(respons_curves_f, device = png, file = paste0("./output/figuren/resp_curves/250m/", x , ".png"))
}

# creating ensemble models:
ensemble_pipeline_250 <- function(soortmodellen_f = SDMs_250m, Env_data_f = Env_data_case, em_choice_f = em_choice, em_by_f = em_by,
                                  em_models_f = em_models, em_selectie_f = em_selectie,
                                  em_thresh_f = em_thresh, em_eval_metrics_f = em_eval_metrics,
                                  em_var_import_f = em_var_import,
                                  em_score_decay_f = em_score_decay)
{ 
  Ensemble_soort_f <- BIOMOD_EnsembleModeling(bm.mod = soortmodellen_f,
                                              models.chosen = em_choice_f,
                                              em.by = em_by_f,
                                              em.algo = em_models_f,
                                              metric.select = em_selectie_f,
                                              metric.select.thresh = em_thresh_f,
                                              metric.eval = em_eval_metrics_f,
                                              var.import = em_var_import_f,
                                              EMwmean.decay = em_score_decay_f,
                                              do.progress = TRUE)
  return(Ensemble_soort_f)
}
# retrieve models from pc
ophalen_ensemble_modellen_250 <- function(){
  adder_f = load("./adder/adder.250_meter.ensemble.models.out")
  gladde_slang_f = load("./gladde.slang/gladde.slang.250_meter.ensemble.models.out")
  hazelworm_f = load("./hazelworm/hazelworm.250_meter.ensemble.models.out")
  levendbarende_hagedis_f = load("./levendbarende.hagedis/levendbarende.hagedis.250_meter.ensemble.models.out")
  muurhagedis_f = load("./muurhagedis/muurhagedis.250_meter.ensemble.models.out")
  ringslang_f = load("./ringslang/ringslang.250_meter.ensemble.models.out")
  zandhagedis_f = load("./zandhagedis/zandhagedis.250_meter.ensemble.models.out")
  
  modellen_f <- list(adder = adder.250_meter.ensemble.models.out, gladde_slang = gladde.slang.250_meter.ensemble.models.out, 
                     hazelworm = hazelworm.250_meter.ensemble.models.out, levendbarende_hagedis = levendbarende.hagedis.250_meter.ensemble.models.out,
                     muurhagedis = muurhagedis.250_meter.ensemble.models.out, ringslang = ringslang.250_meter.ensemble.models.out,
                     zandhagedis = zandhagedis.250_meter.ensemble.models.out)
  return(modellen_f)
}
# evaluation and analysis, saving every element in similar fashion to individual models:
ensemble_model_evaluatie_250 <- function(model_f = ensemble_NL_250m, folder_modelprestatie_f = folder_emodelprestatie) {
  model_eval_f <- data.frame(get_evaluations(model_f), check.rows = TRUE)
  return(model_eval_f)
}

opslaan_ensemble_evaluatie <- function(x){
  edata_f <- ensemble_evaluatie_data[[x]]
  write.csv(edata_f, file = paste0("./output/ens_mod_pres/250m/", x , ".csv"))
}


ens_var_belang_250m <- function(emodel_f = ensemble_NL_250m) {
  ens_var_belang_f <- data.frame(get_variables_importance(emodel_f), check.rows = TRUE)
  return(ens_var_belang_f)
}

ens_opslaan_belang <- function(x){
  ens_belang_f <- ens_var_belang_lijst[[x]]
  write.csv(ens_belang_f, file = paste0("./output/ens_var_belang/250m/", x , ".csv"))
}  

ens_prestatie_box_250 <- function(ensModel_soort_f = ensemble_NL_250m, split_evalbox_f = split_evalbox,
                                  folder_figuren_f = folder_ensemble_figuren) {
  plot_ens_evalmean_f <- bm_PlotEvalMean(bm.out = ensModel_soort_f,
                                         group.by = split_evalbox_f)
  return(plot_ens_evalmean_f)
}

opslaan_ens_box_eval_mean <- function(x, ens_eval_box_f = ens_box_prestatie_lijst){
  ens_box_ev_mean_f <- ens_eval_box_f[[x]]$plot
  ggsave(ens_box_ev_mean_f, device = png, file = paste0("./output/ensemble_figuren/gem_prestatie/250m/", x , ".png"))
}

ens_var_imp_box_250 <- function(ensModel_soort_f = ensemble_NL_250m, ens_split_varbox_f = ens_split_varbox,
                                folder_figuren_f = folder_ensemble_figuren) {
  plot_ens_var_imp_f <- bm_PlotVarImpBoxplot(bm.out = ensModel_soort_f, 
                                             group.by = ens_split_varbox_f)
  return(plot_ens_var_imp_f)
}

opslaan_ens_var_imp_plot <- function(x, ens_var_imp_f = ens_var_imp_plot_lijst){
  ens_box_var_imp_f <- ens_var_imp_f[[x]]$plot
  ggsave(ens_box_var_imp_f, device = png, file = paste0("./output/ensemble_figuren/ensemble_var_belang/250m/", x , ".png"))
}
# response curves:
ens_resp_curves_250 <-  function(ensModel_soort_f = ensemble_NL_250m, fixed_var_rc_f = fixed_var_rc,
                                 folder_ens_figuren_f = folder_ensemble_figuren) {
  plot_ens_rc_f <- bm_PlotResponseCurves(bm.out = ensModel_soort_f, 
                                         models.chosen = 'all',
                                         fixed.var = fixed_var_rc_f)
  return(plot_ens_rc_f)
}

# save:
opslaan_ens_resp_curves <- function(x, curve_lijst_f = ensemble_respons_curves_lijst){
  ens_respons_curves_f <- curve_lijst_f[[x]]$plot
  ggsave(ens_respons_curves_f, device = "png", file = paste0("./output/ensemble_figuren/resp_curves/250m/", x , ".png"))
}
#######################################################################################
# projections, first for whole NL, then case study areas:
ensemble_projecteren_250 <- function(ensModel_soort_f = ensemble_NL_250m, projectie_naam_f =projectie_naam, env_data_f = Env_data, modelkeuze_f = modelkeuze, binary_metric_f = binary_metric, filter_metric_f = filter_metric) {
  ensemble_projectie <- BIOMOD_EnsembleForecasting(bm.em = ensModel_soort_f,
                                                   proj.name = projectie_naam_f,
                                                   new.env = env_data_f,           # huidige data om kaartje NL te krijgen
                                                   models.chosen = modelkeuze_f,
                                                   metric.binary = binary_metric_f,
                                                   metric.filter = filter_metric_f)
  return(ensemble_projectie)
}
# retrieve them from pc:
ophalen_250_projecties_NL <- function(){
  adder_f = load("./adder/proj_projectie_NL_250/adder.projectie_NL_250.ensemble.projection.out")
  gladde_slang_f = load("./gladde.slang/proj_projectie_NL_250/gladde.slang.projectie_NL_250.ensemble.projection.out")
  hazelworm_f = load("./hazelworm/hazelworm.projectie_NL_250.ensemble.projection.out")
  levendbarende_hagedis_f = load("./levendbarende.hagedis/proj_projectie_NL_250/levendbarende.hagedis.projectie_NL_250.ensemble.projection.out")
  muurhagedis_f = load("./muurhagedis/proj_projectie_NL_250/muurhagedis.projectie_NL250.ensemble.projection.out")
  ringslang_f = load("./ringslang/proj_projectie_NL_250/ringslang.projectie_NL_250.ensemble.projection.out")
  zandhagedis_f = load("./zandhagedis/proj_projectie_NL_250/zandhagedis.projectie_NL_250.ensemble.projection.out")
  
  projecties_f <- list(adder = adder.projectie_NL_250.ensemble.projection.out, gladde_slang = gladde.slang.projectie_NL_250.ensemble.projection.out, 
                       hazelworm = hazelworm.projectie_NL_250.ensemble.projection.out, levendbarende_hagedis = levendbarende.hagedis.projectie_NL_250.ensemble.projection.out,
                       muurhagedis = muurhagedis.projectie_NL_250.ensemble.projection.out, ringslang = ringslang.projectie_NL_250.ensemble.projection.out,
                       zandhagedis = zandhagedis.projectie_NL_250.ensemble.projection.out)
  return(projecties_f)
}
# visualisation of projections:
projectie_kaart_250 <- function(proj_f = projecties_NL_250, x){
  projectie_plot <- plot(proj_f,
                         plot.output = 'facet')
  return(projectie_plot)
}
# save as png:
kaarten_opslaan_250 <- function(x, kaarten_soorten_f = kaarten_NL_250, folder_proj_f = 
                                 folder_projecties){
  kaart_soort_f <- kaarten_soorten_f[[x]]
  ggsave(kaart_soort_f, device = "png", file = paste0("./output/projecties/Nederland_250/", x ,".png"),
         height=6, width=8)
}
# save binary maps
range_binair_250 <- function(proj_f = projecties_NL_250, binary_metric_f = binair_methode) {
  range_f <- get_predictions(proj_f,
                             metric.binary = binary_metric_f,
                             model.as.col = TRUE)
  return(range_f)
}
# plot binairy range
plot_range_250 <- function(range_f = range_NL_250){
  range_plots_f <- plot(range_f,
                        main = "geschiktheid habitat")
  return(range_plots_f)
}

# interactive map:
geschiktheid_NL_per_soort <- function(projectie_f = projecties_250, shape_f = NL_shape){
  ens_projectie_adder <- unwrap(projectie_f$adder@proj.out@val)
  ens_projectie_gladde_slang <- unwrap(projectie_f$gladde_slang@proj.out@val)
  ens_projectie_hazelworm <- unwrap(projectie_f$hazelworm@proj.out@val)
  ens_projectie_levendbarende_hagedis <- unwrap(projectie_f$levendbarende_hagedis@proj.out@val)
  ens_projectie_muurhagedis <- unwrap(projectie_f$muurhagedis@proj.out@val)
  ens_projectie_ringslang <- unwrap(projectie_f$ringslang@proj.out@val)
  ens_projectie_zandhagedis <- unwrap(projectie_f$zandhagedis@proj.out@val)
  leaflet() %>%
    addProviderTiles("OpenStreetMap",
                     group = "OpenStreetMap") %>%
    addPolygons(data = shape_f$geometry,
                group = "landsgrens") %>%
    addRasterImage(ens_projectie_adder,
                   group = "adder") %>%
    addRasterImage(ens_projectie_gladde_slang,
                   group = "gladde_slang") %>%
    addRasterImage(ens_projectie_hazelworm,
                   group = "hazelworm")%>%
    addRasterImage(ens_projectie_levendbarende_hagedis,
                   group = "levendbarende_hagedis")%>%
    addRasterImage(ens_projectie_muurhagedis,
                   group = "muurhagedis")%>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "ringslang") %>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "zandhagedis") %>%
    addLayersControl(baseGroups = c('OpenStreetMap'), overlayGroups = c("adder", "gladde_slang",
                                                                        "hazelworm", "levendbarende_hagedis", "muurhagedis", "ringslang", "zandhagedis"))
}
####################################################################################
# Now the same steps per municipality:
# Deurne:
projectie_Deurne <- function(ensModel_soort_f = ensemble_NL_250m, projectie_naam_f =
                               'Deurne_250', modelkeuze_f = modelkeuze, binary_metric_f = 
                               binary_metric, filter_metric_f = filter_metric) {
  ensemble_projectie <- BIOMOD_EnsembleForecasting(bm.em = ensModel_soort_f,
                                                   proj.name = projectie_naam_f,
                                                   new.env = Env_data_Deurne, # data in gemeente Deurne
                                                   models.chosen = modelkeuze_f,
                                                   metric.binary = binary_metric_f,
                                                   metric.filter = filter_metric_f)
  return(ensemble_projectie)
}

proj_kaart_Deurne <- function(proj_f = projectie_Deurne, x, ede_f = shape_Deurne){
  projectie_plot <- plot(proj_f,
                         plot.output = 'facet')
  return(projectie_plot)
}
ophalen_projecties_Deurne <- function(){
  adder_f = load("./adder/proj_Deurne_250/adder.Deurne_250.ensemble.projection.out")
  gladde_slang_f = load("./gladde.slang/proj_Deurne_250/gladde.slang.Deurne_250.ensemble.projection.out")
  hazelworm_f = load("./hazelworm/proj_Deurne_250/hazelworm.Deurne_250.ensemble.projection.out")
  levendbarende_hagedis_f = load("./levendbarende.hagedis/proj_Deurne_250/levendbarende.hagedis.Deurne_250.ensemble.projection.out")
  muurhagedis_f = load("./muurhagedis/proj_Deurne_250/muurhagedis.Deurne_250.ensemble.projection.out")
  ringslang_f = load("./ringslang/proj_Deurne_250/ringslang.Deurne_250.ensemble.projection.out")
  zandhagedis_f = load("./zandhagedis/proj_Deurne_250/zandhagedis.Deurne_250.ensemble.projection.out")
  
  projecties_f <- list(adder = adder.Deurne_250.ensemble.projection.out, gladde_slang = gladde.slang.Deurne_250.ensemble.projection.out, 
                       hazelworm = hazelworm.Deurne_250.ensemble.projection.out, levendbarende_hagedis = levendbarende.hagedis.Deurne_250.ensemble.projection.out,
                       muurhagedis = muurhagedis.Deurne_250.ensemble.projection.out, ringslang = ringslang.Deurne_250.ensemble.projection.out,
                       zandhagedis = zandhagedis.Deurne_250.ensemble.projection.out)
  return(projecties_f)
} 

kaarten_Deurne_opslaan <- function(x, kaarten_soorten_f = kaarten_Deurne){
  kaart_soort_f <- kaarten_soorten_f[[x]]
  ggsave(kaart_soort_f, device = "png", file = paste0("./output/projecties/Deurne_250/", x ,".png"),
         height=6, width=8)
}

range_binair_Deurne <- function(proj_f = projectie_Deurne, binary_metric_f = binair_methode) {
  range_f <- get_predictions(proj_f,
                             metric.binary = binary_metric_f,
                             model.as.col = TRUE)
  return(range_f)
}
:
range_plot_Deurne <- function(range_f = ranges_Deurne){
  range_plots_f <- plot(range_f)
  return(range_plots_f)
}

geschiktheid_Deurne_per_soort <- function(projectie_f = projectie_Deurne, shape_Deurne_f = shape_Deurne){
  ens_projectie_adder <- unwrap(projectie_f$adder@proj.out@val)
  ens_projectie_gladde_slang <- unwrap(projectie_f$gladde_slang@proj.out@val)
  ens_projectie_hazelworm <- unwrap(projectie_f$hazelworm@proj.out@val)
  ens_projectie_levendbarende_hagedis <- unwrap(projectie_f$levendbarende_hagedis@proj.out@val)
  ens_projectie_muurhagedis <- unwrap(projectie_f$muurhagedis@proj.out@val)
  ens_projectie_ringslang <- unwrap(projectie_f$ringslang@proj.out@val)
  ens_projectie_zandhagedis <- unwrap(projectie_f$zandhagedis@proj.out@val)
  leaflet() %>%
    addProviderTiles("OpenStreetMap",
                     group = "OpenStreetMap") %>%
    addPolygons(data = shape_Deurne_f$geometry,
                group = "gemeentegrens") %>%
    addRasterImage(ens_projectie_adder,
                   group = "adder") %>%
    addRasterImage(ens_projectie_gladde_slang,
                   group = "gladde_slang") %>%
    addRasterImage(ens_projectie_hazelworm,
                   group = "hazelworm")%>%
    addRasterImage(ens_projectie_levendbarende_hagedis,
                   group = "levendbarende_hagedis")%>%
    addRasterImage(ens_projectie_muurhagedis,
                   group = "muurhagedis")%>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "ringslang") %>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "zandhagedis") %>%
    addLayersControl(baseGroups = c('OpenStreetMap'), overlayGroups = c("adder", "gladde_slang",
                                                                        "hazelworm", "levendbarende_hagedis", "muurhagedis", "ringslang", "zandhagedis"))
}

# Ede:
projectie_Ede <- function(ensModel_soort_f = ensemble_NL_250m, projectie_naam_f =
                            'Ede_250', modelkeuze_f = modelkeuze, binary_metric_f = 
                            binary_metric, filter_metric_f = filter_metric) {
  ensemble_projectie <- BIOMOD_EnsembleForecasting(bm.em = ensModel_soort_f,
                                                   proj.name = projectie_naam_f,
                                                   new.env = Env_data_Ede, # data in gemeente Ede
                                                   models.chosen = modelkeuze_f,
                                                   metric.binary = binary_metric_f,
                                                   metric.filter = filter_metric_f)
  return(ensemble_projectie)
}

proj_kaart_Ede <- function(proj_f = projectie_Ede, x, ede_f = shape_Ede){
  projectie_plot <- plot(proj_f,
                         plot.output = 'facet')
  return(projectie_plot)
}
ophalen_projecties_Ede <- function(){
  adder_f = load("./adder/proj_Ede_250/adder.Ede_250.ensemble.projection.out")
  gladde_slang_f = load("./gladde.slang/proj_Ede_250/gladde.slang.Ede_250.ensemble.projection.out")
  hazelworm_f = load("./hazelworm/proj_Ede_250/hazelworm.Ede_250.ensemble.projection.out")
  levendbarende_hagedis_f = load("./levendbarende.hagedis/proj_Ede_250/levendbarende.hagedis.Ede_250.ensemble.projection.out")
  muurhagedis_f = load("./muurhagedis/proj_Ede_250/muurhagedis.Ede_250.ensemble.projection.out")
  ringslang_f = load("./ringslang/proj_Ede_250/ringslang.Ede_250.ensemble.projection.out")
  zandhagedis_f = load("./zandhagedis/proj_Ede_250/zandhagedis.Ede_250.ensemble.projection.out")
  
  projecties_f <- list(adder = adder.Ede_250.ensemble.projection.out, gladde_slang = gladde.slang.Ede_250.ensemble.projection.out, 
                       hazelworm = hazelworm.Ede_250.ensemble.projection.out, levendbarende_hagedis = levendbarende.hagedis.Ede_250.ensemble.projection.out,
                       muurhagedis = muurhagedis.Ede_250.ensemble.projection.out, ringslang = ringslang.Ede_250.ensemble.projection.out,
                       zandhagedis = zandhagedis.Ede_250.ensemble.projection.out)
  return(projecties_f)
} 
kaarten_Ede_opslaan <- function(x, kaarten_soorten_f = kaarten_Ede){
  kaart_soort_f <- kaarten_soorten_f[[x]]
  ggsave(kaart_soort_f, device = "png", file = paste0("./output/projecties/Ede_250/", x ,".png"),
         height=6, width=8)
}

range_binair_Ede <- function(proj_f = projectie_Ede, binary_metric_f = binair_methode) {
  range_f <- get_predictions(proj_f,
                             metric.binary = binary_metric_f,
                             model.as.col = TRUE)
  return(range_f)
}

range_plot_Ede <- function(range_f = ranges_Ede){
  range_plots_f <- plot(range_f)
  return(range_plots_f)
}

geschiktheid_Ede_per_soort <- function(projectie_f = projectie_Ede, shape_Ede_f = shape_Ede){
  ens_projectie_adder <- unwrap(projectie_f$adder@proj.out@val)
  ens_projectie_gladde_slang <- unwrap(projectie_f$gladde_slang@proj.out@val)
  ens_projectie_hazelworm <- unwrap(projectie_f$hazelworm@proj.out@val)
  ens_projectie_levendbarende_hagedis <- unwrap(projectie_f$levendbarende_hagedis@proj.out@val)
  ens_projectie_muurhagedis <- unwrap(projectie_f$muurhagedis@proj.out@val)
  ens_projectie_ringslang <- unwrap(projectie_f$ringslang@proj.out@val)
  ens_projectie_zandhagedis <- unwrap(projectie_f$zandhagedis@proj.out@val)
  leaflet() %>%
    addProviderTiles("OpenStreetMap",
                     group = "OpenStreetMap") %>%
    addPolygons(data = shape_Ede_f$geometry,
                group = "gemeentegrens") %>%
    addRasterImage(ens_projectie_adder,
                   group = "adder") %>%
    addRasterImage(ens_projectie_gladde_slang,
                   group = "gladde_slang") %>%
    addRasterImage(ens_projectie_hazelworm,
                   group = "hazelworm")%>%
    addRasterImage(ens_projectie_levendbarende_hagedis,
                   group = "levendbarende_hagedis")%>%
    addRasterImage(ens_projectie_muurhagedis,
                   group = "muurhagedis")%>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "ringslang") %>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "zandhagedis") %>%
    addLayersControl(baseGroups = c('OpenStreetMap'), overlayGroups = c("adder", "gladde_slang",
                                                                        "hazelworm", "levendbarende_hagedis", "muurhagedis", "ringslang", "zandhagedis"))
}
# Epe:
projectie_Epe <- function(ensModel_soort_f = ensemble_NL_250m, projectie_naam_f =
                            'Epe_250', modelkeuze_f = modelkeuze, binary_metric_f = 
                            binary_metric, filter_metric_f = filter_metric) {
  ensemble_projectie <- BIOMOD_EnsembleForecasting(bm.em = ensModel_soort_f,
                                                   proj.name = projectie_naam_f,
                                                   new.env = Env_data_Epe, # data in gemeente Ede
                                                   models.chosen = modelkeuze_f,
                                                   metric.binary = binary_metric_f,
                                                   metric.filter = filter_metric_f)
  return(ensemble_projectie)
}

proj_kaart_Epe <- function(proj_f = projectie_Epe, x, ede_f = shape_Epe){
  projectie_plot <- plot(proj_f,
                         plot.output = 'facet')
  return(projectie_plot)
}
ophalen_projecties_Epe <- function(){
  adder_f = load("./adder/proj_Epe_250/adder.Epe_250.ensemble.projection.out")
  gladde_slang_f = load("./gladde.slang/proj_Epe_250/gladde.slang.Epe_250.ensemble.projection.out")
  hazelworm_f = load("./hazelworm/proj_Epe_250/hazelworm.Epe_250.ensemble.projection.out")
  levendbarende_hagedis_f = load("./levendbarende.hagedis/proj_Epe_250/levendbarende.hagedis.Epe_250.ensemble.projection.out")
  muurhagedis_f = load("./muurhagedis/proj_Epe_250/muurhagedis.Epe_250.ensemble.projection.out")
  ringslang_f = load("./ringslang/proj_Epe_250/ringslang.Epe_250.ensemble.projection.out")
  zandhagedis_f = load("./zandhagedis/proj_Epe_250/zandhagedis.Epe_250.ensemble.projection.out")
  
  projecties_f <- list(adder = adder.Epe_250.ensemble.projection.out, gladde_slang = gladde.slang.Epe_250.ensemble.projection.out, 
                       hazelworm = hazelworm.Epe_250.ensemble.projection.out, levendbarende_hagedis = levendbarende.hagedis.Epe_250.ensemble.projection.out,
                       muurhagedis = muurhagedis.Epe_250.ensemble.projection.out, ringslang = ringslang.Epe_250.ensemble.projection.out,
                       zandhagedis = zandhagedis.Epe_250.ensemble.projection.out)
  return(projecties_f)
} 

kaarten_Epe_opslaan <- function(x, kaarten_soorten_f = kaarten_Epe){
  kaart_soort_f <- kaarten_soorten_f[[x]]
  ggsave(kaart_soort_f, device = "png", file = paste0("./output/projecties/Epe_250/", x ,".png"),
         height=6, width=8)
}

range_binair_Epe <- function(proj_f = projectie_Epe, binary_metric_f = binair_methode) {
  range_f <- get_predictions(proj_f,
                             metric.binary = binary_metric_f,
                             model.as.col = TRUE)
  return(range_f)
}

range_plot_Epe <- function(range_f = ranges_Epe){
  range_plots_f <- plot(range_f)
  return(range_plots_f)
}

geschiktheid_Epe_per_soort <- function(projectie_f = projectie_Epe, shape_Epe_f = shape_Epe){
  ens_projectie_adder <- unwrap(projectie_f$adder@proj.out@val)
  ens_projectie_gladde_slang <- unwrap(projectie_f$gladde_slang@proj.out@val)
  ens_projectie_hazelworm <- unwrap(projectie_f$hazelworm@proj.out@val)
  ens_projectie_levendbarende_hagedis <- unwrap(projectie_f$levendbarende_hagedis@proj.out@val)
  ens_projectie_muurhagedis <- unwrap(projectie_f$muurhagedis@proj.out@val)
  ens_projectie_ringslang <- unwrap(projectie_f$ringslang@proj.out@val)
  ens_projectie_zandhagedis <- unwrap(projectie_f$zandhagedis@proj.out@val)
  leaflet() %>%
    addProviderTiles("OpenStreetMap",
                     group = "OpenStreetMap") %>%
    addPolygons(data = shape_Epe_f$geometry,
                group = "gemeentegrens") %>%
    addRasterImage(ens_projectie_adder,
                   group = "adder") %>%
    addRasterImage(ens_projectie_gladde_slang,
                   group = "gladde_slang") %>%
    addRasterImage(ens_projectie_hazelworm,
                   group = "hazelworm")%>%
    addRasterImage(ens_projectie_levendbarende_hagedis,
                   group = "levendbarende_hagedis")%>%
    addRasterImage(ens_projectie_muurhagedis,
                   group = "muurhagedis")%>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "ringslang") %>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "zandhagedis") %>%
    addLayersControl(baseGroups = c('OpenStreetMap'), overlayGroups = c("adder", "gladde_slang",
                                                                        "hazelworm", "levendbarende_hagedis", "muurhagedis", "ringslang", "zandhagedis"))
}

# Horst:
projectie_Horst <- function(ensModel_soort_f = ensemble_NL_250m, projectie_naam_f =
                              'Horst_250', modelkeuze_f = modelkeuze, binary_metric_f = 
                              binary_metric, filter_metric_f = filter_metric) {
  ensemble_projectie <- BIOMOD_EnsembleForecasting(bm.em = ensModel_soort_f,
                                                   proj.name = projectie_naam_f,
                                                   new.env = Env_data_Horst, # data in gemeente Ede
                                                   models.chosen = modelkeuze_f,
                                                   metric.binary = binary_metric_f,
                                                   metric.filter = filter_metric_f)
  return(ensemble_projectie)
}

proj_kaart_Horst <- function(proj_f = projectie_Horst, x, ede_f = shape_Horst){
  projectie_plot <- plot(proj_f,
                         plot.output = 'facet')
  return(projectie_plot)
}
ophalen_projecties_Horst <- function(){
  adder_f = load("./adder/proj_Horst_250/adder.Horst_250.ensemble.projection.out")
  gladde_slang_f = load("./gladde.slang/proj_Horst_250/gladde.slang.Horst_250.ensemble.projection.out")
  hazelworm_f = load("./hazelworm/proj_Horst_250/hazelworm.Horst_250.ensemble.projection.out")
  levendbarende_hagedis_f = load("./levendbarende.hagedis/proj_Horst_250/levendbarende.hagedis.Horst_250.ensemble.projection.out")
  muurhagedis_f = load("./muurhagedis/proj_Horst_250/muurhagedis.Horst_250.ensemble.projection.out")
  ringslang_f = load("./ringslang/proj_Horst_250/ringslang.Horst_250.ensemble.projection.out")
  zandhagedis_f = load("./zandhagedis/proj_Horst_250/zandhagedis.Horst_250.ensemble.projection.out")
  
  projecties_f <- list(adder = adder.Horst_250.ensemble.projection.out, gladde_slang = gladde.slang.Horst_250.ensemble.projection.out, 
                       hazelworm = hazelworm.Horst_250.ensemble.projection.out, levendbarende_hagedis = levendbarende.hagedis.Horst_250.ensemble.projection.out,
                       muurhagedis = muurhagedis.Horst_250.ensemble.projection.out, ringslang = ringslang.Horst_250.ensemble.projection.out,
                       zandhagedis = zandhagedis.Horst_250.ensemble.projection.out)
  return(projecties_f)
} 

kaarten_Horst_opslaan <- function(x, kaarten_soorten_f = kaarten_Horst){
  kaart_soort_f <- kaarten_soorten_f[[x]]
  ggsave(kaart_soort_f, device = "png", file = paste0("./output/projecties/Horst_250/", x ,".png"),
         height=6, width=8)
}

range_binair_Horst <- function(proj_f = projectie_Horst, binary_metric_f = binair_methode) {
  range_f <- get_predictions(proj_f,
                             metric.binary = binary_metric_f,
                             model.as.col = TRUE)
  return(range_f)
}

range_plot_Horst <- function(range_f = ranges_Horst){
  range_plots_f <- plot(range_f)
  return(range_plots_f)
}


geschiktheid_Horst_per_soort <- function(projectie_f = projectie_Horst, shape_Horst_f = shape_Horst){
  ens_projectie_adder <- unwrap(projectie_f$adder@proj.out@val)
  ens_projectie_gladde_slang <- unwrap(projectie_f$gladde_slang@proj.out@val)
  ens_projectie_hazelworm <- unwrap(projectie_f$hazelworm@proj.out@val)
  ens_projectie_levendbarende_hagedis <- unwrap(projectie_f$levendbarende_hagedis@proj.out@val)
  ens_projectie_muurhagedis <- unwrap(projectie_f$muurhagedis@proj.out@val)
  ens_projectie_ringslang <- unwrap(projectie_f$ringslang@proj.out@val)
  ens_projectie_zandhagedis <- unwrap(projectie_f$zandhagedis@proj.out@val)
  leaflet() %>%
    addProviderTiles("OpenStreetMap",
                     group = "OpenStreetMap") %>%
    addPolygons(data = shape_Horst_f$geometry,
                group = "gemeentegrens") %>%
    addRasterImage(ens_projectie_adder,
                   group = "adder") %>%
    addRasterImage(ens_projectie_gladde_slang,
                   group = "gladde_slang") %>%
    addRasterImage(ens_projectie_hazelworm,
                   group = "hazelworm")%>%
    addRasterImage(ens_projectie_levendbarende_hagedis,
                   group = "levendbarende_hagedis")%>%
    addRasterImage(ens_projectie_muurhagedis,
                   group = "muurhagedis")%>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "ringslang") %>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "zandhagedis") %>%
    addLayersControl(baseGroups = c('OpenStreetMap'), overlayGroups = c("adder", "gladde_slang",
                                                                        "hazelworm", "levendbarende_hagedis", "muurhagedis", "ringslang", "zandhagedis"))
}

# Huizen: 
projectie_Huizen <- function(ensModel_soort_f = ensemble_NL_250m, projectie_naam_f =
                               'Huizen_250', modelkeuze_f = modelkeuze, binary_metric_f = 
                               binary_metric, filter_metric_f = filter_metric) {
  ensemble_projectie <- BIOMOD_EnsembleForecasting(bm.em = ensModel_soort_f,
                                                   proj.name = projectie_naam_f,
                                                   new.env = Env_data_Huizen, # data in gemeente Ede
                                                   models.chosen = modelkeuze_f,
                                                   metric.binary = binary_metric_f,
                                                   metric.filter = filter_metric_f)
  return(ensemble_projectie)
}

proj_kaart_Huizen <- function(proj_f = projectie_Huizen, x, ede_f = shape_Huizen){
  projectie_plot <- plot(proj_f,
                         plot.output = 'facet')
  return(projectie_plot)
}
ophalen_projecties_Huizen <- function(){
  adder_f = load("./adder/proj_Huizen_250/adder.Huizen_250.ensemble.projection.out")
  gladde_slang_f = load("./gladde.slang/proj_Huizen_250/gladde.slang.Huizen_250.ensemble.projection.out")
  hazelworm_f = load("./hazelworm/proj_Huizen_250/hazelworm.Huizen_250.ensemble.projection.out")
  levendbarende_hagedis_f = load("./levendbarende.hagedis/proj_Huizen_250/levendbarende.hagedis.Huizen_250.ensemble.projection.out")
  muurhagedis_f = load("./muurhagedis/proj_Huizen_250/muurhagedis.Huizen_250.ensemble.projection.out")
  ringslang_f = load("./ringslang/proj_Huizen_250/ringslang.Huizen_250.ensemble.projection.out")
  zandhagedis_f = load("./zandhagedis/proj_Huizen_250/zandhagedis.Huizen_250.ensemble.projection.out")
  
  projecties_f <- list(adder = adder.Huizen_250.ensemble.projection.out, gladde_slang = gladde.slang.Huizen_250.ensemble.projection.out, 
                       hazelworm = hazelworm.Huizen_250.ensemble.projection.out, levendbarende_hagedis = levendbarende.hagedis.Huizen_250.ensemble.projection.out,
                       muurhagedis = muurhagedis_250.Huizen.ensemble.projection.out, ringslang = ringslang.Huizen_250.ensemble.projection.out,
                       zandhagedis = zandhagedis.Huizen_250.ensemble.projection.out)
  return(projecties_f)
} 

kaarten_Huizen_opslaan <- function(x, kaarten_soorten_f = kaarten_Huizen){
  kaart_soort_f <- kaarten_soorten_f[[x]]
  ggsave(kaart_soort_f, device = "png", file = paste0("./output/projecties/Huizen_250/", x ,".png"),
         height=6, width=8)
}

range_binair_Huizen <- function(proj_f = projectie_Huizen, binary_metric_f = binair_methode) {
  range_f <- get_predictions(proj_f,
                             metric.binary = binary_metric_f,
                             model.as.col = TRUE)
  return(range_f)
}

range_plot_Huizen <- function(range_f = ranges_Huizen){
  range_plots_f <- plot(range_f)
  return(range_plots_f)
}

geschiktheid_Huizen_per_soort <- function(projectie_f = projectie_Huizen, shape_Huizen_f = shape_Huizen){
  ens_projectie_adder <- unwrap(projectie_f$adder@proj.out@val)
  ens_projectie_gladde_slang <- unwrap(projectie_f$gladde_slang@proj.out@val)
  ens_projectie_hazelworm <- unwrap(projectie_f$hazelworm@proj.out@val)
  ens_projectie_levendbarende_hagedis <- unwrap(projectie_f$levendbarende_hagedis@proj.out@val)
  ens_projectie_muurhagedis <- unwrap(projectie_f$muurhagedis@proj.out@val)
  ens_projectie_ringslang <- unwrap(projectie_f$ringslang@proj.out@val)
  ens_projectie_zandhagedis <- unwrap(projectie_f$zandhagedis@proj.out@val)
  leaflet() %>%
    addProviderTiles("OpenStreetMap",
                     group = "OpenStreetMap") %>%
    addPolygons(data = shape_Huizen_f$geometry,
                group = "gemeentegrens") %>%
    addRasterImage(ens_projectie_adder,
                   group = "adder") %>%
    addRasterImage(ens_projectie_gladde_slang,
                   group = "gladde_slang") %>%
    addRasterImage(ens_projectie_hazelworm,
                   group = "hazelworm")%>%
    addRasterImage(ens_projectie_levendbarende_hagedis,
                   group = "levendbarende_hagedis")%>%
    addRasterImage(ens_projectie_muurhagedis,
                   group = "muurhagedis")%>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "ringslang") %>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "zandhagedis") %>%
    addLayersControl(baseGroups = c('OpenStreetMap'), overlayGroups = c("adder", "gladde_slang",
                                                                        "hazelworm", "levendbarende_hagedis", "muurhagedis", "ringslang", "zandhagedis"))
}

# Ommen:
projectie_Ommen <- function(ensModel_soort_f = ensemble_NL_250m, projectie_naam_f =
                              'Ommen_250', modelkeuze_f = modelkeuze, binary_metric_f = 
                              binary_metric, filter_metric_f = filter_metric) {
  ensemble_projectie <- BIOMOD_EnsembleForecasting(bm.em = ensModel_soort_f,
                                                   proj.name = projectie_naam_f,
                                                   new.env = Env_data_Ommen, # data in gemeente Ede
                                                   models.chosen = modelkeuze_f,
                                                   metric.binary = binary_metric_f,
                                                   metric.filter = filter_metric_f)
  return(ensemble_projectie)
}

proj_kaart_Ommen <- function(proj_f = projectie_Ommen, x, ede_f = shape_Ommen){
  projectie_plot <- plot(proj_f,
                         plot.output = 'facet')
  return(projectie_plot)
}
ophalen_projecties_Ommen <- function(){
  adder_f = load("./adder/proj_Ommen_250/adder.Ommen_250.ensemble.projection.out")
  gladde_slang_f = load("./gladde.slang/proj_Ommen_250/gladde.slang.Ommen_250.ensemble.projection.out")
  hazelworm_f = load("./hazelworm/proj_Ommen_250/hazelworm.Ommen_250.ensemble.projection.out")
  levendbarende_hagedis_f = load("./levendbarende.hagedis/proj_Ommen_250/levendbarende.hagedis.Ommen_250.ensemble.projection.out")
  muurhagedis_f = load("./muurhagedis/proj_Ommen_250/muurhagedis.Ommen_250.ensemble.projection.out")
  ringslang_f = load("./ringslang/proj_Ommen_250/ringslang.Ommen_250.ensemble.projection.out")
  zandhagedis_f = load("./zandhagedis/proj_Ommen_250/zandhagedis.Ommen_250.ensemble.projection.out")
  
  projecties_f <- list(adder = adder.Ommen_250.ensemble.projection.out, gladde_slang = gladde.slang.Ommen_250.ensemble.projection.out, 
                       hazelworm = hazelworm.Ommen_250.ensemble.projection.out, levendbarende_hagedis = levendbarende.hagedis.Ommen_250.ensemble.projection.out,
                       muurhagedis = muurhagedis.Ommen_250.ensemble.projection.out, ringslang = ringslang.Ommen_250.ensemble.projection.out,
                       zandhagedis = zandhagedis.Ommen_250.ensemble.projection.out)
  return(projecties_f)
} 

kaarten_Ommen_opslaan <- function(x, kaarten_soorten_f = kaarten_Ommen){
  kaart_soort_f <- kaarten_soorten_f[[x]]
  ggsave(kaart_soort_f, device = "png", file = paste0("./output/projecties/Ommen_250/", x ,".png"),
         height=6, width=8)
}

range_binair_Ommen <- function(proj_f = projectie_Ommen, binary_metric_f = binair_methode) {
  range_f <- get_predictions(proj_f,
                             metric.binary = binary_metric_f,
                             model.as.col = TRUE)
  return(range_f)
}

range_plot_Ommen <- function(range_f = ranges_Ommen){
  range_plots_f <- plot(range_f)
  return(range_plots_f)
}

geschiktheid_Ommen_per_soort <- function(projectie_f = projectie_Ommen, shape_Ommen_f = shape_Ommen){
  ens_projectie_adder <- unwrap(projectie_f$adder@proj.out@val)
  ens_projectie_gladde_slang <- unwrap(projectie_f$gladde_slang@proj.out@val)
  ens_projectie_hazelworm <- unwrap(projectie_f$hazelworm@proj.out@val)
  ens_projectie_levendbarende_hagedis <- unwrap(projectie_f$levendbarende_hagedis@proj.out@val)
  ens_projectie_muurhagedis <- unwrap(projectie_f$muurhagedis@proj.out@val)
  ens_projectie_ringslang <- unwrap(projectie_f$ringslang@proj.out@val)
  ens_projectie_zandhagedis <- unwrap(projectie_f$zandhagedis@proj.out@val)
  leaflet() %>%
    addProviderTiles("OpenStreetMap",
                     group = "OpenStreetMap") %>%
    addPolygons(data = shape_Ommen_f$geometry,
                group = "gemeentegrens") %>%
    addRasterImage(ens_projectie_adder,
                   group = "adder") %>%
    addRasterImage(ens_projectie_gladde_slang,
                   group = "gladde_slang") %>%
    addRasterImage(ens_projectie_hazelworm,
                   group = "hazelworm")%>%
    addRasterImage(ens_projectie_levendbarende_hagedis,
                   group = "levendbarende_hagedis")%>%
    addRasterImage(ens_projectie_muurhagedis,
                   group = "muurhagedis")%>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "ringslang") %>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "zandhagedis") %>%
    addLayersControl(baseGroups = c('OpenStreetMap'), overlayGroups = c("adder", "gladde_slang",
                                                                        "hazelworm", "levendbarende_hagedis", "muurhagedis", "ringslang", "zandhagedis"))
}

# Roermond:
projectie_Roermond <- function(ensModel_soort_f = ensemble_NL_250m, projectie_naam_f =
                                 'Roermond_250', modelkeuze_f = modelkeuze, binary_metric_f = 
                                 binary_metric, filter_metric_f = filter_metric) {
  ensemble_projectie <- BIOMOD_EnsembleForecasting(bm.em = ensModel_soort_f,
                                                   proj.name = projectie_naam_f,
                                                   new.env = Env_data_Roermond, # data in gemeente Ede
                                                   models.chosen = modelkeuze_f,
                                                   metric.binary = binary_metric_f,
                                                   metric.filter = filter_metric_f)
  return(ensemble_projectie)
}

proj_kaart_Roermond <- function(proj_f = projectie_Roermond, x, ede_f = shape_Roermond){
  projectie_plot <- plot(proj_f,
                         plot.output = 'facet')
  return(projectie_plot)
}
ophalen_projecties_Roermond <- function(){
  adder_f = load("./adder/proj_Roermond_250/adder.Roermond_250.ensemble.projection.out")
  gladde_slang_f = load("./gladde.slang/proj_Roermond_250/gladde.slang.Roermond_250.ensemble.projection.out")
  hazelworm_f = load("./hazelworm/proj_Roermond_250/hazelworm.Roermond_250.ensemble.projection.out")
  levendbarende_hagedis_f = load("./levendbarende.hagedis/proj_Roermond_250/levendbarende.hagedis.Roermond_250.ensemble.projection.out")
  muurhagedis_f = load("./muurhagedis/proj_Roermond_250/muurhagedis.Roermond_250.ensemble.projection.out")
  ringslang_f = load("./ringslang/proj_Roermond_250/ringslang.Roermond_250.ensemble.projection.out")
  zandhagedis_f = load("./zandhagedis/proj_Roermond_250/zandhagedis.Roermond_250.ensemble.projection.out")
  
  projecties_f <- list(adder = adder.Roermond_250.ensemble.projection.out, gladde_slang = gladde.slang.Roermond_250.ensemble.projection.out, 
                       hazelworm = hazelworm.Roermond_250.ensemble.projection.out, levendbarende_hagedis = levendbarende.hagedis.Roermond_250.ensemble.projection.out,
                       muurhagedis = muurhagedis.Roermond_250.ensemble.projection.out, ringslang = ringslang.Roermond_250.ensemble.projection.out,
                       zandhagedis = zandhagedis.Roermond_250.ensemble.projection.out)
  return(projecties_f)
} 
kaarten_Roermond_opslaan <- function(x, kaarten_soorten_f = kaarten_Roermond){
  kaart_soort_f <- kaarten_soorten_f[[x]]
  ggsave(kaart_soort_f, device = "png", file = paste0("./output/projecties/Roermond_250/", x ,".png"),
         height=6, width=8)
}
range_binair_Roermond <- function(proj_f = projectie_Roermond, binary_metric_f = binair_methode) {
  range_f <- get_predictions(proj_f,
                             metric.binary = binary_metric_f,
                             model.as.col = TRUE)
  return(range_f)
}

range_plot_Roermond <- function(range_f = ranges_Roermond){
  range_plots_f <- plot(range_f)
  return(range_plots_f)
}

geschiktheid_Roermond_per_soort <- function(projectie_f = projectie_Roermond, shape_Roermond_f = shape_Roermond){
  ens_projectie_adder <- unwrap(projectie_f$adder@proj.out@val)
  ens_projectie_gladde_slang <- unwrap(projectie_f$gladde_slang@proj.out@val)
  ens_projectie_hazelworm <- unwrap(projectie_f$hazelworm@proj.out@val)
  ens_projectie_levendbarende_hagedis <- unwrap(projectie_f$levendbarende_hagedis@proj.out@val)
  ens_projectie_muurhagedis <- unwrap(projectie_f$muurhagedis@proj.out@val)
  ens_projectie_ringslang <- unwrap(projectie_f$ringslang@proj.out@val)
  ens_projectie_zandhagedis <- unwrap(projectie_f$zandhagedis@proj.out@val)
  leaflet() %>%
    addProviderTiles("OpenStreetMap",
                     group = "OpenStreetMap") %>%
    addPolygons(data = shape_Roermond_f$geometry,
                group = "gemeentegrens") %>%
    addRasterImage(ens_projectie_adder,
                   group = "adder") %>%
    addRasterImage(ens_projectie_gladde_slang,
                   group = "gladde_slang") %>%
    addRasterImage(ens_projectie_hazelworm,
                   group = "hazelworm")%>%
    addRasterImage(ens_projectie_levendbarende_hagedis,
                   group = "levendbarende_hagedis")%>%
    addRasterImage(ens_projectie_muurhagedis,
                   group = "muurhagedis")%>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "ringslang") %>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "zandhagedis") %>%
    addLayersControl(baseGroups = c('OpenStreetMap'), overlayGroups = c("adder", "gladde_slang",
                                                                        "hazelworm", "levendbarende_hagedis", "muurhagedis", "ringslang", "zandhagedis"))
}

# Westerveld:
projectie_Westerveld <- function(ensModel_soort_f = ensemble_NL_250m, projectie_naam_f =
                                   'Westerveld_250', modelkeuze_f = modelkeuze, binary_metric_f = 
                                   binary_metric, filter_metric_f = filter_metric) {
  ensemble_projectie <- BIOMOD_EnsembleForecasting(bm.em = ensModel_soort_f,
                                                   proj.name = projectie_naam_f,
                                                   new.env = Env_data_Westerveld, # data in gemeente Ede
                                                   models.chosen = modelkeuze_f,
                                                   metric.binary = binary_metric_f,
                                                   metric.filter = filter_metric_f)
  return(ensemble_projectie)
}
proj_kaart_Westerveld <- function(proj_f = projectie_Westerveld, x, ede_f = shape_Westerveld){
  projectie_plot <- plot(proj_f,
                         plot.output = 'facet')
  return(projectie_plot)
}
ophalen_projecties_Westerveld <- function(){
  adder_f = load("./adder/proj_Westerveld_250/adder.Westerveld_250.ensemble.projection.out")
  gladde_slang_f = load("./gladde.slang/proj_Westerveld_250/gladde.slang.Westerveld_250.ensemble.projection.out")
  hazelworm_f = load("./hazelworm/proj_Westerveld_250/hazelworm.Westerveld_250.ensemble.projection.out")
  levendbarende_hagedis_f = load("./levendbarende.hagedis/proj_Westerveld_250/levendbarende.hagedis.Westerveld_250.ensemble.projection.out")
  muurhagedis_f = load("./muurhagedis/proj_Westerveld_250/muurhagedis.Westerveld_250.ensemble.projection.out")
  ringslang_f = load("./ringslang/proj_Westerveld_250/ringslang.Westerveld_250.ensemble.projection.out")
  zandhagedis_f = load("./zandhagedis/proj_Westerveld_250/zandhagedis.Westerveld_250.ensemble.projection.out")
  
  projecties_f <- list(adder = adder.Westerveld_250.ensemble.projection.out, gladde_slang = gladde.slang.Westerveld_250.ensemble.projection.out, 
                       hazelworm = hazelworm.Westerveld_250.ensemble.projection.out, levendbarende_hagedis = levendbarende.hagedis.Westerveld_250.ensemble.projection.out,
                       muurhagedis = muurhagedis.Westerveld_250.ensemble.projection.out, ringslang = ringslang.Westerveld_250.ensemble.projection.out,
                       zandhagedis = zandhagedis.Westerveld_250.ensemble.projection.out)
  return(projecties_f)
} 
kaarten_Westerveld_opslaan <- function(x, kaarten_soorten_f = kaarten_Westerveld){
  kaart_soort_f <- kaarten_soorten_f[[x]]
  ggsave(kaart_soort_f, device = "png", file = paste0("./output/projecties/Westerveld_250/", x ,".png"),
         height=6, width=8)
}

range_binair_Westerveld <- function(proj_f = projectie_Westerveld, binary_metric_f = binair_methode) {
  range_f <- get_predictions(proj_f,
                             metric.binary = binary_metric_f,
                             model.as.col = TRUE)
  return(range_f)
}

range_plot_Westerveld <- function(range_f = ranges_Westerveld){
  range_plots_f <- plot(range_f)
  return(range_plots_f)
}

geschiktheid_Westerveld_per_soort <- function(projectie_f = projectie_Westerveld, shape_Westerveld_f = shape_Westerveld){
  ens_projectie_adder <- unwrap(projectie_f$adder@proj.out@val)
  ens_projectie_gladde_slang <- unwrap(projectie_f$gladde_slang@proj.out@val)
  ens_projectie_hazelworm <- unwrap(projectie_f$hazelworm@proj.out@val)
  ens_projectie_levendbarende_hagedis <- unwrap(projectie_f$levendbarende_hagedis@proj.out@val)
  ens_projectie_muurhagedis <- unwrap(projectie_f$muurhagedis@proj.out@val)
  ens_projectie_ringslang <- unwrap(projectie_f$ringslang@proj.out@val)
  ens_projectie_zandhagedis <- unwrap(projectie_f$zandhagedis@proj.out@val)
  leaflet() %>%
    addProviderTiles("OpenStreetMap",
                     group = "OpenStreetMap") %>%
    addPolygons(data = shape_Westerveld_f$geometry,
                group = "gemeentegrens") %>%
    addRasterImage(ens_projectie_adder,
                   group = "adder") %>%
    addRasterImage(ens_projectie_gladde_slang,
                   group = "gladde_slang") %>%
    addRasterImage(ens_projectie_hazelworm,
                   group = "hazelworm")%>%
    addRasterImage(ens_projectie_levendbarende_hagedis,
                   group = "levendbarende_hagedis")%>%
    addRasterImage(ens_projectie_muurhagedis,
                   group = "muurhagedis")%>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "ringslang") %>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "zandhagedis") %>%
    addLayersControl(baseGroups = c('OpenStreetMap'), overlayGroups = c("adder", "gladde_slang",
                                                                        "hazelworm", "levendbarende_hagedis", "muurhagedis", "ringslang", "zandhagedis"))
}

# Wormerland:
projectie_Wormerland <- function(ensModel_soort_f = ensemble_NL_250m, projectie_naam_f =
                                   'Wormerland_250', modelkeuze_f = modelkeuze, binary_metric_f = 
                                   binary_metric, filter_metric_f = filter_metric) {
  ensemble_projectie <- BIOMOD_EnsembleForecasting(bm.em = ensModel_soort_f,
                                                   proj.name = projectie_naam_f,
                                                   new.env = Env_data_Wormerland, # data in gemeente Ede
                                                   models.chosen = modelkeuze_f,
                                                   metric.binary = binary_metric_f,
                                                   metric.filter = filter_metric_f)
  return(ensemble_projectie)
}
proj_kaart_Wormerland <- function(proj_f = projectie_Wormerland, x, ede_f = shape_Wormerland){
  projectie_plot <- plot(proj_f,
                         plot.output = 'facet')
  return(projectie_plot)
}
ophalen_projecties_Wormerland <- function(){
  adder_f = load("./adder/proj_Wormerland_250/adder.Wormerland_250.ensemble.projection.out")
  gladde_slang_f = load("./gladde.slang/proj_Wormerland_250/gladde.slang.Wormerland_250.ensemble.projection.out")
  hazelworm_f = load("./hazelworm/proj_Wormerland_250/hazelworm.Wormerland_250.ensemble.projection.out")
  levendbarende_hagedis_f = load("./levendbarende.hagedis/proj_Wormerland_250/levendbarende.hagedis.Wormerland_250.ensemble.projection.out")
  muurhagedis_f = load("./muurhagedis/proj_Wormerland_250/muurhagedis.Wormerland_250.ensemble.projection.out")
  ringslang_f = load("./ringslang/proj_Wormerland_250/ringslang.Wormerland_250.ensemble.projection.out")
  zandhagedis_f = load("./zandhagedis/proj_Wormerland_250/zandhagedis.Wormerland_250.ensemble.projection.out")
  
  projecties_f <- list(adder = adder.Wormerland_250.ensemble.projection.out, gladde_slang = gladde.slang.Wormerland_250.ensemble.projection.out, 
                       hazelworm = hazelworm.Wormerland_250.ensemble.projection.out, levendbarende_hagedis = levendbarende.hagedis_250.Wormerland.ensemble.projection.out,
                       muurhagedis = muurhagedis.Wormerland_250.ensemble.projection.out, ringslang = ringslang.Wormerland_250.ensemble.projection.out,
                       zandhagedis = zandhagedis.Wormerland_250.ensemble.projection.out)
  return(projecties_f)
} 
kaarten_Wormerland_opslaan <- function(x, kaarten_soorten_f = kaarten_Wormerland){
  kaart_soort_f <- kaarten_soorten_f[[x]]
  ggsave(kaart_soort_f, device = "png", file = paste0("./output/projecties/Wormerland_250/", x ,".png"),
         height=6, width=8)
}

range_binair_Wormerland <- function(proj_f = projectie_Wormerland, binary_metric_f = binair_methode) {
  range_f <- get_predictions(proj_f,
                             metric.binary = binary_metric_f,
                             model.as.col = TRUE)
  return(range_f)
}

range_plot_Wormerland <- function(range_f = ranges_Wormerland){
  range_plots_f <- plot(range_f)
  return(range_plots_f)
}

geschiktheid_Wormerland_per_soort <- function(projectie_f = projectie_Wormerland, shape_Wormerland_f = shape_Wormerland){
  ens_projectie_adder <- unwrap(projectie_f$adder@proj.out@val)
  ens_projectie_gladde_slang <- unwrap(projectie_f$gladde_slang@proj.out@val)
  ens_projectie_hazelworm <- unwrap(projectie_f$hazelworm@proj.out@val)
  ens_projectie_levendbarende_hagedis <- unwrap(projectie_f$levendbarende_hagedis@proj.out@val)
  ens_projectie_muurhagedis <- unwrap(projectie_f$muurhagedis@proj.out@val)
  ens_projectie_ringslang <- unwrap(projectie_f$ringslang@proj.out@val)
  ens_projectie_zandhagedis <- unwrap(projectie_f$zandhagedis@proj.out@val)
  leaflet() %>%
    addProviderTiles("OpenStreetMap",
                     group = "OpenStreetMap") %>%
    addPolygons(data = shape_Wormerland_f$geometry,
                group = "gemeentegrens") %>%
    addRasterImage(ens_projectie_adder,
                   group = "adder") %>%
    addRasterImage(ens_projectie_gladde_slang,
                   group = "gladde_slang") %>%
    addRasterImage(ens_projectie_hazelworm,
                   group = "hazelworm")%>%
    addRasterImage(ens_projectie_levendbarende_hagedis,
                   group = "levendbarende_hagedis")%>%
    addRasterImage(ens_projectie_muurhagedis,
                   group = "muurhagedis")%>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "ringslang") %>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "zandhagedis") %>%
    addLayersControl(baseGroups = c('OpenStreetMap'), overlayGroups = c("adder", "gladde_slang",
                                                                        "hazelworm", "levendbarende_hagedis", "muurhagedis", "ringslang", "zandhagedis"))
}

# Zutphen:
projectie_Zutphen <- function(ensModel_soort_f = ensemble_NL_250m, projectie_naam_f =
                                'Zutphen_250', modelkeuze_f = modelkeuze, binary_metric_f = 
                                binary_metric, filter_metric_f = filter_metric) {
  ensemble_projectie <- BIOMOD_EnsembleForecasting(bm.em = ensModel_soort_f,
                                                   proj.name = projectie_naam_f,
                                                   new.env = Env_data_Zutphen, # data in gemeente Ede
                                                   models.chosen = modelkeuze_f,
                                                   metric.binary = binary_metric_f,
                                                   metric.filter = filter_metric_f)
  return(ensemble_projectie)
}
proj_kaart_Zutphen <- function(proj_f = projectie_Zutphen, x, ede_f = shape_Zutphen){
  projectie_plot <- plot(proj_f,
                         plot.output = 'facet')
  return(projectie_plot)
}
ophalen_projecties_Zutphen <- function(){
  adder_f = load("./adder/proj_Zutphen_250/adder.Zutphen_250.ensemble.projection.out")
  gladde_slang_f = load("./gladde.slang/proj_Zutphen_250/gladde.slang.Zutphen_250.ensemble.projection.out")
  hazelworm_f = load("./hazelworm/proj_Zutphen_250/hazelworm.Zutphen_250.ensemble.projection.out")
  levendbarende_hagedis_f = load("./levendbarende.hagedis/proj_Zutphen_250/levendbarende.hagedis.Zutphen_250.ensemble.projection.out")
  muurhagedis_f = load("./muurhagedis/proj_Zutphen_250/muurhagedis.Zutphen_250.ensemble.projection.out")
  ringslang_f = load("./ringslang/proj_Zutphen_250/ringslang.Zutphen_250.ensemble.projection.out")
  zandhagedis_f = load("./zandhagedis/proj_Zutphen_250/zandhagedis.Zutphen_250.ensemble.projection.out")
  
  projecties_f <- list(adder = adder.Zutphen_250.ensemble.projection.out, gladde_slang = gladde.slang.Zutphen_250.ensemble.projection.out, 
                       hazelworm = hazelworm.Zutphen_250.ensemble.projection.out, levendbarende_hagedis = levendbarende.hagedis.Zutphen_250.ensemble.projection.out,
                       muurhagedis = muurhagedis.Zutphen_250.ensemble.projection.out, ringslang = ringslang.Zutphen_250.ensemble.projection.out,
                       zandhagedis = zandhagedis.Zutphen_250.ensemble.projection.out)
  return(projecties_f)
} 
kaarten_Zutphen_opslaan <- function(x, kaarten_soorten_f = kaarten_Zutphen){
  kaart_soort_f <- kaarten_soorten_f[[x]]
  ggsave(kaart_soort_f, device = "png", file = paste0("./output/projecties/Zutphen_250/", x ,".png"),
         height=6, width=8)
}
range_binair_Zutphen <- function(proj_f = projectie_Zutphen, binary_metric_f = binair_methode) {
  range_f <- get_predictions(proj_f,
                             metric.binary = binary_metric_f,
                             model.as.col = TRUE)
  return(range_f)
}
range_plot_Zutphen <- function(range_f = ranges_Zutphen){
  range_plots_f <- plot(range_f)
  return(range_plots_f)
}

geschiktheid_Zutphen_per_soort <- function(projectie_f = projectie_Zutphen, shape_Zutphen_f = shape_Zutphen){
  ens_projectie_adder <- unwrap(projectie_f$adder@proj.out@val)
  ens_projectie_gladde_slang <- unwrap(projectie_f$gladde_slang@proj.out@val)
  ens_projectie_hazelworm <- unwrap(projectie_f$hazelworm@proj.out@val)
  ens_projectie_levendbarende_hagedis <- unwrap(projectie_f$levendbarende_hagedis@proj.out@val)
  ens_projectie_muurhagedis <- unwrap(projectie_f$muurhagedis@proj.out@val)
  ens_projectie_ringslang <- unwrap(projectie_f$ringslang@proj.out@val)
  ens_projectie_zandhagedis <- unwrap(projectie_f$zandhagedis@proj.out@val)
  leaflet() %>%
    addProviderTiles("OpenStreetMap",
                     group = "OpenStreetMap") %>%
    addPolygons(data = shape_Zutphen_f$geometry,
                group = "gemeentegrens") %>%
    addRasterImage(ens_projectie_adder,
                   group = "adder") %>%
    addRasterImage(ens_projectie_gladde_slang,
                   group = "gladde_slang") %>%
    addRasterImage(ens_projectie_hazelworm,
                   group = "hazelworm")%>%
    addRasterImage(ens_projectie_levendbarende_hagedis,
                   group = "levendbarende_hagedis")%>%
    addRasterImage(ens_projectie_muurhagedis,
                   group = "muurhagedis")%>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "ringslang") %>%
    addRasterImage(ens_projectie_zandhagedis,
                   group = "zandhagedis") %>%
    addLayersControl(baseGroups = c('OpenStreetMap'), overlayGroups = c("adder", "gladde_slang",
                                                                        "hazelworm", "levendbarende_hagedis", "muurhagedis", "ringslang", "zandhagedis"))
}
####################################################################################################################
# functions related to testing the performance:
# Retrieving and rasterising PAs per municipalities
uitknippen_PA_ede <- function(pa_f = PA_per_soort,gebied_f = shape_Ede) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  pa_f[,1] <- 0
  raster_f <- rasterize(x = pa_f[2:3], y = leeg_f, field = pa_f[,1], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
uitknippen_PA_epe <- function(pa_f = PA_per_soort,gebied_f = shape_Epe) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  pa_f[,1] <- 0
  raster_f <- rasterize(x = pa_f[2:3], y = leeg_f, field = pa_f[,1], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
uitknippen_PA_deurne <- function(pa_f = PA_per_soort,gebied_f = shape_Deurne) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  pa_f[,1] <- 0
  raster_f <- rasterize(x = pa_f[2:3], y = leeg_f, field = pa_f[,1], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
uitknippen_PA_roermond <- function(pa_f = PA_per_soort,gebied_f = shape_Roermond) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  pa_f[,1] <- 0
  raster_f <- rasterize(x = pa_f[2:3], y = leeg_f, field = pa_f[,1], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
uitknippen_PA_ommen <- function(pa_f = PA_per_soort,gebied_f = shape_Ommen) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  pa_f[,1] <- 0
  raster_f <- rasterize(x = pa_f[2:3], y = leeg_f, field = pa_f[,1], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
uitknippen_PA_huizen <- function(pa_f = PA_per_soort,gebied_f = shape_Huizen) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  pa_f[,1] <- 0
  raster_f <- rasterize(x = pa_f[2:3], y = leeg_f, field = pa_f[,1], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
uitknippen_PA_westerveld <- function(pa_f = PA_per_soort,gebied_f = shape_Westerveld) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  pa_f[,1] <- 0
  raster_f <- rasterize(x = pa_f[2:3], y = leeg_f, field = pa_f[,1], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
uitknippen_PA_wormerland <- function(pa_f = PA_per_soort,gebied_f = shape_Wormerland) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  pa_f[,1] <- 0
  raster_f <- rasterize(x = pa_f[2:3], y = leeg_f, field = pa_f[,1], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
uitknippen_PA_horst <- function(pa_f = PA_per_soort,gebied_f = shape_Horst) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  pa_f[,1] <- 0
  raster_f <- rasterize(x = pa_f[2:3], y = leeg_f, field = pa_f[,1], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
uitknippen_PA_zutphen <- function(pa_f = PA_per_soort,gebied_f = shape_Zutphen) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  pa_f[,1] <- 0
  raster_f <- rasterize(x = pa_f[2:3], y = leeg_f, field = pa_f[,1], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
####################################################################################
# Rasterising occurrence data:
raster_maken_waarnemingen_ede <- function(w_f = NDFF_ede_per_soort, gebied_f = shape_Ede) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  raster_f <- rasterize(x = w_f[3:4], y = leeg_f, field = w_f[,2], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
raster_maken_waarnemingen_epe <- function(w_f = NDFF_epe_per_soort, gebied_f = shape_Epe) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  raster_f <- rasterize(x = w_f[3:4], y = leeg_f, field = w_f[,2], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
raster_maken_waarnemingen_deurne <- function(w_f = NDFF_deurne_per_soort, gebied_f = shape_Deurne) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  raster_f <- rasterize(x = w_f[3:4], y = leeg_f, field = w_f[,2], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
raster_maken_waarnemingen_roermond <- function(w_f = NDFF_roermond_per_soort, gebied_f = shape_Roermond) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  raster_f <- rasterize(x = w_f[3:4], y = leeg_f, field = w_f[,2], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
raster_maken_waarnemingen_ommen <- function(w_f = NDFF_ommen_per_soort, gebied_f = shape_Ommen) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  raster_f <- rasterize(x = w_f[3:4], y = leeg_f, field = w_f[,2], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
raster_maken_waarnemingen_huizen <- function(w_f = NDFF_huizen_per_soort, gebied_f = shape_Huizen) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  raster_f <- rasterize(x = w_f[3:4], y = leeg_f, field = w_f[,2], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
raster_maken_waarnemingen_westerveld <- function(w_f = NDFF_westerveld_per_soort, gebied_f = shape_Westerveld) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  raster_f <- rasterize(x = w_f[3:4], y = leeg_f, field = w_f[,2], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
raster_maken_waarnemingen_wormerland <- function(w_f = NDFF_wormerland_per_soort, gebied_f = shape_Wormerland) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  raster_f <- rasterize(x = w_f[3:4], y = leeg_f, field = w_f[,2], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
raster_maken_waarnemingen_horst <- function(w_f = NDFF_horst_per_soort, gebied_f = shape_Horst) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  raster_f <- rasterize(x = w_f[3:4], y = leeg_f, field = w_f[,2], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
raster_maken_waarnemingen_zutphen <- function(w_f = NDFF_zutphen_per_soort, gebied_f = shape_Zutphen) {
  leeg_f <- raster(xmn = 3.083333, xmx = 7.224998, ymn = 50.75, ymx = 53.74167,
                   resolution = c(0.004166665, 0.004166665), crs = "+proj=longlat +datum=WGS84 +no_defs")
  raster_f <- rasterize(x = w_f[3:4], y = leeg_f, field = w_f[,2], fun = mean)
  raster_f <- crop(raster_f, gebied_f)
  return(raster_f)
}
##################################################################################################################
# Rasterise predictions:
# Ede:
raster_projectie_ede <- function(voorspelling_f = range_NL_250, gebied_f = shape_Ede)
{
  raster_voorspelling_f <- raster(voorspelling_f)
  raster_gebied_f <- crop(raster_voorspelling_f, gebied_f)
  return(raster_gebied_f)
}
# Epe:
raster_projectie_epe <- function(voorspelling_f = range_NL_250, gebied_f = shape_Epe)
{
  raster_voorspelling_f <- raster(voorspelling_f)
  raster_gebied_f <- crop(raster_voorspelling_f, gebied_f)
  return(raster_gebied_f)
}
# Deurne:
raster_projectie_deurne <- function(voorspelling_f = range_NL_250, gebied_f = shape_Deurne)
{
  raster_voorspelling_f <- raster(voorspelling_f)
  raster_gebied_f <- crop(raster_voorspelling_f, gebied_f)
  return(raster_gebied_f)
}
# Horst:
raster_projectie_horst <- function(voorspelling_f = range_NL_250, gebied_f = shape_Horst)
{
  raster_voorspelling_f <- raster(voorspelling_f)
  raster_gebied_f <- crop(raster_voorspelling_f, gebied_f)
  return(raster_gebied_f)
}
# Huizen:
raster_projectie_huizen <- function(voorspelling_f = range_NL_250, gebied_f = shape_Huizen)
{
  raster_voorspelling_f <- raster(voorspelling_f)
  raster_gebied_f <- crop(raster_voorspelling_f, gebied_f)
  return(raster_gebied_f)
}
# Ommen:
raster_projectie_ommen <- function(voorspelling_f = range_NL_250, gebied_f = shape_Ommen)
{
  raster_voorspelling_f <- raster(voorspelling_f)
  raster_gebied_f <- crop(raster_voorspelling_f, gebied_f)
  return(raster_gebied_f)
}
# Roermond:
raster_projectie_roermond <- function(voorspelling_f = range_NL_250, gebied_f = shape_Roermond)
{
  raster_voorspelling_f <- raster(voorspelling_f)
  raster_gebied_f <- crop(raster_voorspelling_f, gebied_f)
  return(raster_gebied_f)
}
# Westerveld:
raster_projectie_westerveld <- function(voorspelling_f = range_NL_250, gebied_f = shape_Westerveld)
{
  raster_voorspelling_f <- raster(voorspelling_f)
  raster_gebied_f <- crop(raster_voorspelling_f, gebied_f)
  return(raster_gebied_f)
}
# Wormerland:
raster_projectie_wormerland <- function(voorspelling_f = range_NL_250, gebied_f = shape_Wormerland)
{
  raster_voorspelling_f <- raster(voorspelling_f)
  raster_gebied_f <- crop(raster_voorspelling_f, gebied_f)
  return(raster_gebied_f)
}
# Zutphen:
raster_projectie_zutphen <- function(voorspelling_f = range_NL_250, gebied_f = shape_Zutphen)
{
  raster_voorspelling_f <- raster(voorspelling_f)
  raster_gebied_f <- crop(raster_voorspelling_f, gebied_f)
  return(raster_gebied_f)
}
# Merge PA and occurrence rasters per species per municipality:
merge_ede <- function(pa_f = PA_ede_raster, w_f = waarnemingen_rasters_Ede)
{
  adder_ede = raster::merge(w_f$adder, pa_f$adder)
  gladde_slang_ede = raster::merge(w_f$`gladde slang`, pa_f$gladde_slang)
  hazelworm_ede = raster::merge(w_f$hazelworm, pa_f$hazelworm) 
  levendbarende_hagedis_ede = raster::merge(w_f$`levendbarende hagedis`, pa_f$levendbarende_hagedis)
  muurhagedis_ede = pa_f$muurhagedis
  ringslang_ede = raster::merge(w_f$ringslang, pa_f$ringslang)
  zandhagedis_ede = raster::merge(w_f$zandhagedis, pa_f$zandhagedis)
  list_ede_f <- list(adder = adder_ede, gladde_slang = gladde_slang_ede, hazelworm = hazelworm_ede,levendbarende_hagedis = levendbarende_hagedis_ede, muurhagedis = muurhagedis_ede,ringslang = ringslang_ede, zandhagedis = zandhagedis_ede)  
  return(list_ede_f)
}
merge_epe <- function(pa_f = PA_epe_raster, w_f = waarnemingen_rasters_Epe)
{
  adder_epe = raster::merge(w_f$adder, pa_f$adder)
  gladde_slang_epe = raster::merge(w_f$`gladde slang`, pa_f$gladde_slang)
  hazelworm_epe = raster::merge(w_f$hazelworm, pa_f$hazelworm) 
  levendbarende_hagedis_epe = raster::merge(w_f$`levendbarende hagedis`, pa_f$levendbarende_hagedis)
  muurhagedis_epe = pa_f$muurhagedis
  ringslang_epe = raster::merge(w_f$ringslang, pa_f$ringslang)
  zandhagedis_epe = raster::merge(w_f$zandhagedis, pa_f$zandhagedis)
  list_epe_f <- list(adder = adder_epe, gladde_slang = gladde_slang_epe, hazelworm = hazelworm_epe,
                     levendbarende_hagedis = levendbarende_hagedis_epe, muurhagedis = muurhagedis_epe,
                     ringslang = ringslang_epe, zandhagedis = zandhagedis_epe)  
  return(list_epe_f)
}
merge_deurne <- function(pa_f = PA_deurne_raster, w_f = waarnemingen_rasters_Deurne)
{
  adder_deurne = pa_f$adder
  gladde_slang_deurne = raster::merge(w_f$`gladde slang`, pa_f$gladde_slang)
  hazelworm_deurne = pa_f$hazelworm
  levendbarende_hagedis_deurne = raster::merge(w_f$`levendbarende hagedis`, pa_f$levendbarende_hagedis)
  muurhagedis_deurne = pa_f$muurhagedis
  ringslang_deurne = pa_f$ringslang
  zandhagedis_deurne = raster::merge(w_f$zandhagedis, pa_f$zandhagedis)
  list_deurne_f <- list(adder = adder_deurne, gladde_slang = gladde_slang_deurne, hazelworm = hazelworm_deurne,
                        levendbarende_hagedis = levendbarende_hagedis_deurne, muurhagedis = muurhagedis_deurne,
                        ringslang = ringslang_deurne, zandhagedis = zandhagedis_deurne)  
  return(list_deurne_f)
}
merge_roermond <- function(pa_f = PA_roermond_raster, w_f = waarnemingen_rasters_Roermond)
{
  adder_roermond = raster::merge(w_f$adder, pa_f$adder)
  gladde_slang_roermond = pa_f$gladde_slang
  hazelworm_roermond = raster::merge(w_f$hazelworm, pa_f$hazelworm) 
  levendbarende_hagedis_roermond = raster::merge(w_f$`levendbarende hagedis`, pa_f$levendbarende_hagedis)
  muurhagedis_roermond = pa_f$muurhagedis
  ringslang_roermond = raster::merge(w_f$ringslang, pa_f$ringslang)
  zandhagedis_roermond = raster::merge(w_f$zandhagedis, pa_f$zandhagedis)
  list_roermond_f <- list(adder = adder_roermond, gladde_slang = gladde_slang_roermond, hazelworm = hazelworm_roermond,
                          levendbarende_hagedis = levendbarende_hagedis_roermond, muurhagedis = muurhagedis_roermond,
                          ringslang = ringslang_roermond, zandhagedis = zandhagedis_roermond)  
  return(list_roermond_f)
}
merge_huizen <- function(pa_f = PA_huizen_raster, w_f = waarnemingen_rasters_Huizen)
{
  adder_huizen = pa_f$adder
  gladde_slang_huizen = pa_f$gladde_slang
  hazelworm_huizen = pa_f$hazelworm
  levendbarende_hagedis_huizen = raster::merge(w_f$`levendbarende hagedis`, pa_f$levendbarende_hagedis)
  muurhagedis_huizen = pa_f$muurhagedis
  ringslang_huizen = raster::merge(w_f$ringslang, pa_f$ringslang)
  zandhagedis_huizen = raster::merge(w_f$zandhagedis, pa_f$zandhagedis)
  list_huizen_f <- list(adder = adder_huizen, gladde_slang = gladde_slang_huizen, hazelworm = hazelworm_huizen,
                        levendbarende_hagedis = levendbarende_hagedis_huizen, muurhagedis = muurhagedis_huizen,
                        ringslang = ringslang_huizen, zandhagedis = zandhagedis_huizen)  
  return(list_huizen_f)
}
merge_ommen <- function(pa_f = PA_ommen_raster, w_f = waarnemingen_rasters_Ommen)
{
  adder_ommen = raster::merge(w_f$adder, pa_f$adder)
  gladde_slang_ommen = pa_f$gladde_slang
  hazelworm_ommen = raster::merge(w_f$hazelworm, pa_f$hazelworm) 
  levendbarende_hagedis_ommen = raster::merge(w_f$`levendbarende hagedis`, pa_f$levendbarende_hagedis)
  muurhagedis_ommen = pa_f$muurhagedis
  ringslang_ommen = raster::merge(w_f$ringslang, pa_f$ringslang)
  zandhagedis_ommen = raster::merge(w_f$zandhagedis, pa_f$zandhagedis)
  list_ommen_f <- list(adder = adder_ommen, gladde_slang = gladde_slang_ommen, hazelworm = hazelworm_ommen,
                       levendbarende_hagedis = levendbarende_hagedis_ommen, muurhagedis = muurhagedis_ommen,
                       ringslang = ringslang_ommen, zandhagedis = zandhagedis_ommen)  
  return(list_ommen_f)
}
merge_westerveld <- function(pa_f = PA_westerveld_raster, w_f = waarnemingen_rasters_Westerveld)
{
  adder_PA_westerveld = raster::merge(w_f$adder, pa_f$adder)
  gladde_slang_PA_westerveld = raster::merge(w_f$`gladde slang`, pa_f$gladde_slang)
  hazelworm_PA_westerveld = raster::merge(w_f$hazelworm, pa_f$hazelworm) 
  levendbarende_hagedis_PA_westerveld = raster::merge(w_f$`levendbarende hagedis`, pa_f$levendbarende_hagedis)
  muurhagedis_PA_westerveld = merge(w_f$muurhagedis, pa_f$muurhagedis)
  ringslang_PA_westerveld = raster::merge(w_f$ringslang, pa_f$ringslang)
  zandhagedis_PA_westerveld = raster::merge(w_f$zandhagedis, pa_f$zandhagedis)
  list_PA_westerveld_f <- list(adder = adder_PA_westerveld, gladde_slang = gladde_slang_PA_westerveld, hazelworm = hazelworm_PA_westerveld,
                               levendbarende_hagedis = levendbarende_hagedis_PA_westerveld, muurhagedis = muurhagedis_PA_westerveld,
                               ringslang = ringslang_PA_westerveld, zandhagedis = zandhagedis_PA_westerveld)  
  return(list_PA_westerveld_f)
}
merge_wormerland <- function(pa_f = PA_wormerland_raster, w_f = waarnemingen_rasters_Wormerland)
{
  adder_wormerland = pa_f$adder
  gladde_slang_wormerland = pa_f$gladde_slang
  hazelworm_wormerland = pa_f$hazelworm
  levendbarende_hagedis_wormerland = pa_f$levendbarende_hagedis
  muurhagedis_wormerland = pa_f$muurhagedis
  ringslang_wormerland = raster::merge(w_f$ringslang, pa_f$ringslang)
  zandhagedis_wormerland = pa_f$zandhagedis
  list_wormerland_f <- list(adder = adder_wormerland, gladde_slang = gladde_slang_wormerland, hazelworm = hazelworm_wormerland,
                            levendbarende_hagedis = levendbarende_hagedis_wormerland, muurhagedis = muurhagedis_wormerland,
                            ringslang = ringslang_wormerland, zandhagedis = zandhagedis_wormerland)  
  return(list_wormerland_f)
}
merge_horst <- function(pa_f = PA_horst_raster, w_f = waarnemingen_rasters_Horst)
{
  adder_horst = pa_f$adder
  gladde_slang_horst = raster::merge(w_f$`gladde slang`, pa_f$gladde_slang)
  hazelworm_horst = pa_f$hazelworm
  levendbarende_hagedis_horst = raster::merge(w_f$`levendbarende hagedis`, pa_f$levendbarende_hagedis)
  muurhagedis_horst = pa_f$muurhagedis
  ringslang_horst = pa_f$ringslang
  zandhagedis_horst = pa_f$zandhagedis
  list_horst_f <- list(adder = adder_horst, gladde_slang = gladde_slang_horst, hazelworm = hazelworm_horst,
                       levendbarende_hagedis = levendbarende_hagedis_horst, muurhagedis = muurhagedis_horst,
                       ringslang = ringslang_horst, zandhagedis = zandhagedis_horst)  
  return(list_horst_f)
}
merge_zutphen <- function(pa_f = PA_zutphen_raster, w_f = waarnemingen_rasters_Zutphen)
{
  adder_zutphen = pa_f$adder
  gladde_slang_zutphen = pa_f$gladde_slang
  hazelworm_zutphen = raster::merge(w_f$hazelworm, pa_f$hazelworm) 
  levendbarende_hagedis_zutphen = raster::merge(w_f$`levendbarende hagedis`, pa_f$levendbarende_hagedis)
  muurhagedis_zutphen = pa_f$muurhagedis
  ringslang_zutphen = raster::merge(w_f$ringslang, pa_f$ringslang)
  zandhagedis_zutphen = pa_f$zandhagedis
  list_zutphen_f <- list(adder = adder_zutphen, gladde_slang = gladde_slang_zutphen, hazelworm = hazelworm_zutphen,
                         levendbarende_hagedis = levendbarende_hagedis_zutphen, muurhagedis = muurhagedis_zutphen,
                         ringslang = ringslang_zutphen, zandhagedis = zandhagedis_zutphen)  
  return(list_zutphen_f)
}
#################################################################################################################
#Extract data from the rasters and save them per municipality and species:
# Ede:
df_data_ede <- function(raster_f = PA_W_ede)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(adder = adder, gladde_slang = gladde_slang, hazelworm = hazelworm,
               levendbarende_hagedis = levendbarende_hagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(df_f)
}
# Epe:
df_data_epe <- function(raster_f = PA_W_epe)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(adder = adder, gladde_slang = gladde_slang, hazelworm = hazelworm,
               levendbarende_hagedis = levendbarende_hagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(df_f)
}
# Deurne:
df_data_deurne <- function(raster_f = PA_W_deurne)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(gladde_slang = gladde_slang, levendbarende_hagedis = levendbarende_hagedis,
               zandhagedis = zandhagedis)
  return(df_f)
}
# Horst:
df_data_horst <- function(raster_f = PA_W_horst)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(gladde_slang = gladde_slang, levendbarende_hagedis = levendbarende_hagedis)
  return(df_f)
}
# Huizen:
df_data_huizen <- function(raster_f = PA_W_huizen)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(levendbarende_hagedis = levendbarende_hagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(df_f)
}
# Ommen:
df_data_ommen <- function(raster_f = PA_W_ommen)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(adder = adder, hazelworm = hazelworm,
               levendbarende_hagedis = levendbarende_hagedis, 
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(df_f)
}
# Westerveld:
df_data_westerveld <- function(raster_f = PA_W_westerveld)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(adder = adder, gladde_slang = gladde_slang, hazelworm = hazelworm,
               levendbarende_hagedis = levendbarende_hagedis, muurhagedis = muurhagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(df_f)
}
# Wormerland:
df_data_wormerland <- function(raster_f = PA_W_wormerland)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(ringslang = ringslang)
  return(df_f)
}
# Zutphen:
df_data_zutphen <- function(raster_f = PA_W_zutphen)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(hazelworm = hazelworm, levendbarende_hagedis = levendbarende_hagedis, 
               ringslang = ringslang)
  return(df_f)
}
# Roermond:
df_data_roermond <- function(raster_f = PA_W_roermond)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(adder = adder, hazelworm = hazelworm,
               levendbarende_hagedis = levendbarende_hagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(df_f)
}
# do the same for the projections raster:
# Ede:
df_voorspeld_ede <- function(raster_f = raster_voorspelling_Ede)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(adder = adder, gladde_slang = gladde_slang, hazelworm = hazelworm,
               levendbarende_hagedis = levendbarende_hagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(df_f)
}
# Epe:
df_voorspeld_epe <- function(raster_f = raster_voorspelling_Epe)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(adder = adder, gladde_slang = gladde_slang, hazelworm = hazelworm,
               levendbarende_hagedis = levendbarende_hagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(df_f)
}
# Deurne:
df_voorspeld_deurne <- function(raster_f = raster_voorspelling_Deurne)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(gladde_slang = gladde_slang, levendbarende_hagedis = levendbarende_hagedis, 
               zandhagedis = zandhagedis)
  return(df_f)
}
# Horst:
df_voorspeld_horst <- function(raster_f = raster_voorspelling_Horst)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(gladde_slang = gladde_slang,
               levendbarende_hagedis = levendbarende_hagedis)
  return(df_f)
}
# Huizen:
df_voorspeld_huizen <- function(raster_f = raster_voorspelling_Huizen)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(levendbarende_hagedis = levendbarende_hagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(df_f)
}
# Ommen:
df_voorspeld_ommen <- function(raster_f = raster_voorspelling_Ommen)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(adder = adder, hazelworm = hazelworm,
               levendbarende_hagedis = levendbarende_hagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(df_f)
}
# Roermond:
df_voorspeld_roermond <- function(raster_f = raster_voorspelling_Roermond)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(adder = adder, hazelworm = hazelworm,
               levendbarende_hagedis = levendbarende_hagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(df_f)
}
# Westerveld:
df_voorspeld_westerveld <- function(raster_f = raster_voorspelling_Westerveld)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(adder = adder, gladde_slang = gladde_slang, hazelworm = hazelworm,
               levendbarende_hagedis = levendbarende_hagedis, muurhagedis = muurhagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(df_f)
}
# Wormerland:
df_voorspeld_wormerland <- function(raster_f = raster_voorspelling_Wormerland)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(ringslang = ringslang)
  return(df_f)
}
# Zutphen:
df_voorspeld_zutphen <- function(raster_f = raster_voorspelling_Zutphen)
{ 
  adder <- as.data.frame(values(raster_f$adder))
  colnames(adder) = "aanwezigheid"
  adder$aanwezigheid <- as.factor(adder$aanwezigheid)
  
  gladde_slang <- as.data.frame(values(raster_f$gladde_slang))
  colnames(gladde_slang) = "aanwezigheid"
  gladde_slang$aanwezigheid = as.factor(gladde_slang$aanwezigheid)
  
  hazelworm <- as.data.frame(values(raster_f$hazelworm))
  colnames(hazelworm) = "aanwezigheid"
  hazelworm$aanwezigheid = as.factor(hazelworm$aanwezigheid)
  
  levendbarende_hagedis <- as.data.frame(values(raster_f$levendbarende_hagedis))
  colnames(levendbarende_hagedis) = "aanwezigheid"
  levendbarende_hagedis$aanwezigheid = as.factor(levendbarende_hagedis$aanwezigheid)
  
  muurhagedis <- as.data.frame(values(raster_f$muurhagedis))
  colnames(muurhagedis) = "aanwezigheid"
  muurhagedis$aanwezigheid = as.factor(muurhagedis$aanwezigheid)
  
  ringslang <- as.data.frame(values(raster_f$ringslang))
  colnames(ringslang) = "aanwezigheid"
  ringslang$aanwezigheid = as.factor(ringslang$aanwezigheid)
  
  zandhagedis <- as.data.frame(values(raster_f$zandhagedis))
  colnames(zandhagedis) = "aanwezigheid"
  zandhagedis$aanwezigheid = as.factor(zandhagedis$aanwezigheid)
  df_f <- list(hazelworm = hazelworm, levendbarende_hagedis = levendbarende_hagedis,
               ringslang = ringslang)
  return(df_f)
}
################################################################################
# Generate confusion matrices to compare predictions against real data:
# Ede:
cm_ede <- function(echt_f = echte_data_Ede, voorspeld_f = voorspelling_Ede)
{
  adder <- confusionMatrix(data = voorspeld_f$adder$aanwezigheid, 
                           reference = echt_f$adder$aanwezigheid,
                           positive = "1")
  gladde_slang <- confusionMatrix(data = voorspeld_f$gladde_slang$aanwezigheid, 
                                  reference = echt_f$gladde_slang$aanwezigheid,
                                  positive = "1")
  hazelworm <- confusionMatrix(data = voorspeld_f$hazelworm$aanwezigheid, 
                               reference = echt_f$hazelworm$aanwezigheid,
                               positive = "1")
  levendbarende_hagedis <- confusionMatrix(data = voorspeld_f$levendbarende_hagedis$aanwezigheid, 
                                           reference = echt_f$levendbarende_hagedis$aanwezigheid,
                                           positive = "1")
  ringslang <- confusionMatrix(data = voorspeld_f$ringslang$aanwezigheid, 
                               reference = echt_f$ringslang$aanwezigheid,
                               positive = "1")
  zandhagedis <- confusionMatrix(data = voorspeld_f$zandhagedis$aanwezigheid, 
                                 reference = echt_f$zandhagedis$aanwezigheid,
                                 positive = "1")
  cm_f <- list(adder = adder, gladde_slang = gladde_slang, hazelworm = hazelworm,
               levendbarende_hagedis = levendbarende_hagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(cm_f)
}
# Epe:
cm_epe <- function(echt_f = echte_data_Epe, voorspeld_f = voorspelling_Epe)
{
  adder <- confusionMatrix(data = voorspeld_f$adder$aanwezigheid, 
                           reference = echt_f$adder$aanwezigheid,
                           positive = "1")
  gladde_slang <- confusionMatrix(data = voorspeld_f$gladde_slang$aanwezigheid, 
                                  reference = echt_f$gladde_slang$aanwezigheid,
                                  positive = "1")
  hazelworm <- confusionMatrix(data = voorspeld_f$hazelworm$aanwezigheid, 
                               reference = echt_f$hazelworm$aanwezigheid,
                               positive = "1")
  levendbarende_hagedis <- confusionMatrix(data = voorspeld_f$levendbarende_hagedis$aanwezigheid, 
                                           reference = echt_f$levendbarende_hagedis$aanwezigheid,
                                           positive = "1")
  ringslang <- confusionMatrix(data = voorspeld_f$ringslang$aanwezigheid, 
                               reference = echt_f$ringslang$aanwezigheid,
                               positive = "1")
  zandhagedis <- confusionMatrix(data = voorspeld_f$zandhagedis$aanwezigheid, 
                                 reference = echt_f$zandhagedis$aanwezigheid,
                                 positive = "1")
  cm_f <- list(adder = adder, gladde_slang = gladde_slang, hazelworm = hazelworm,
               levendbarende_hagedis = levendbarende_hagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(cm_f)
}
# Deurne:
cm_deurne <- function(echt_f = echte_data_Deurne, voorspeld_f = voorspelling_Deurne)
{
  gladde_slang <- confusionMatrix(data = voorspeld_f$gladde_slang$aanwezigheid, 
                                  reference = echt_f$gladde_slang$aanwezigheid,
                                  positive = "1")
  levendbarende_hagedis <- confusionMatrix(data = voorspeld_f$levendbarende_hagedis$aanwezigheid, 
                                           reference = echt_f$levendbarende_hagedis$aanwezigheid,
                                           positive = "1")
  zandhagedis <- confusionMatrix(data = voorspeld_f$zandhagedis$aanwezigheid, 
                                 reference = echt_f$zandhagedis$aanwezigheid,
                                 positive = "1")
  cm_f <- list(gladde_slang = gladde_slang,
               levendbarende_hagedis = levendbarende_hagedis,
               zandhagedis = zandhagedis)
  return(cm_f)
}
# Horst:
cm_horst <- function(echt_f = echte_data_Horst, voorspeld_f = voorspelling_Horst)
{
  gladde_slang <- confusionMatrix(data = voorspeld_f$gladde_slang$aanwezigheid, 
                                  reference = echt_f$gladde_slang$aanwezigheid,
                                  positive = "1")
  levendbarende_hagedis <- confusionMatrix(data = voorspeld_f$levendbarende_hagedis$aanwezigheid, 
                                           reference = echt_f$levendbarende_hagedis$aanwezigheid,
                                           positive = "1")
  cm_f <- list(gladde_slang = gladde_slang,
               levendbarende_hagedis = levendbarende_hagedis)
  return(cm_f)
}
# Huizen: 
cm_huizen <- function(echt_f = echte_data_Huizen, voorspeld_f = voorspelling_Huizen)
{
  levendbarende_hagedis <- confusionMatrix(data = voorspeld_f$levendbarende_hagedis$aanwezigheid, 
                                           reference = echt_f$levendbarende_hagedis$aanwezigheid,
                                           positive = "1")
  ringslang <- confusionMatrix(data = voorspeld_f$ringslang$aanwezigheid, 
                               reference = echt_f$ringslang$aanwezigheid,
                               positive = "1")
  zandhagedis <- confusionMatrix(data = voorspeld_f$zandhagedis$aanwezigheid, 
                                 reference = echt_f$zandhagedis$aanwezigheid,
                                 positive = "1")
  cm_f <- list(levendbarende_hagedis = levendbarende_hagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(cm_f)
}
# Ommen:
cm_ommen <- function(echt_f = echte_data_Ommen, voorspeld_f = voorspelling_Ommen)
{
  adder <- confusionMatrix(data = voorspeld_f$adder$aanwezigheid, 
                           reference = echt_f$adder$aanwezigheid,
                           positive = "1")
  hazelworm <- confusionMatrix(data = voorspeld_f$hazelworm$aanwezigheid, 
                               reference = echt_f$hazelworm$aanwezigheid,
                               positive = "1")
  levendbarende_hagedis <- confusionMatrix(data = voorspeld_f$levendbarende_hagedis$aanwezigheid, 
                                           reference = echt_f$levendbarende_hagedis$aanwezigheid,
                                           positive = "1")
  ringslang <- confusionMatrix(data = voorspeld_f$ringslang$aanwezigheid, 
                               reference = echt_f$ringslang$aanwezigheid,
                               positive = "1")
  zandhagedis <- confusionMatrix(data = voorspeld_f$zandhagedis$aanwezigheid, 
                                 reference = echt_f$zandhagedis$aanwezigheid,
                                 positive = "1")
  cm_f <- list(adder = adder, hazelworm = hazelworm,
               levendbarende_hagedis = levendbarende_hagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(cm_f)
}
# Roermond:
cm_roermond <- function(echt_f = echte_data_Roermond, voorspeld_f = voorspelling_Roermond)
{
  adder <- confusionMatrix(data = voorspeld_f$adder$aanwezigheid, 
                           reference = echt_f$adder$aanwezigheid,
                           positive = "1")
  hazelworm <- confusionMatrix(data = voorspeld_f$hazelworm$aanwezigheid, 
                               reference = echt_f$hazelworm$aanwezigheid,
                               positive = "1")
  levendbarende_hagedis <- confusionMatrix(data = voorspeld_f$levendbarende_hagedis$aanwezigheid, 
                                           reference = echt_f$levendbarende_hagedis$aanwezigheid,
                                           positive = "1")
  ringslang <- confusionMatrix(data = voorspeld_f$ringslang$aanwezigheid, 
                               reference = echt_f$ringslang$aanwezigheid,
                               positive = "1")
  zandhagedis <- confusionMatrix(data = voorspeld_f$zandhagedis$aanwezigheid, 
                                 reference = echt_f$zandhagedis$aanwezigheid,
                                 positive = "1")
  cm_f <- list(adder = adder, hazelworm = hazelworm,
               levendbarende_hagedis = levendbarende_hagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(cm_f)
}
# Westerveld:
cm_westerveld <- function(echt_f = echte_data_Westerveld, voorspeld_f = voorspelling_Westerveld)
{
  adder <- confusionMatrix(data = voorspeld_f$adder$aanwezigheid, 
                           reference = echt_f$adder$aanwezigheid,
                           positive = "1")
  gladde_slang <- confusionMatrix(data = voorspeld_f$gladde_slang$aanwezigheid, 
                                  reference = echt_f$gladde_slang$aanwezigheid,
                                  positive = "1")
  hazelworm <- confusionMatrix(data = voorspeld_f$hazelworm$aanwezigheid, 
                               reference = echt_f$hazelworm$aanwezigheid,
                               positive = "1")
  levendbarende_hagedis <- confusionMatrix(data = voorspeld_f$levendbarende_hagedis$aanwezigheid, 
                                           reference = echt_f$levendbarende_hagedis$aanwezigheid,
                                           positive = "1")
  muurhagedis <- confusionMatrix(data = voorspeld_f$muurhagedis$aanwezigheid, 
                                 reference = echt_f$muurhagedis$aanwezigheid,
                                 positive = "1")
  ringslang <- confusionMatrix(data = voorspeld_f$ringslang$aanwezigheid, 
                               reference = echt_f$ringslang$aanwezigheid,
                               positive = "1")
  zandhagedis <- confusionMatrix(data = voorspeld_f$zandhagedis$aanwezigheid, 
                                 reference = echt_f$zandhagedis$aanwezigheid,
                                 positive = "1")
  cm_f <- list(adder = adder, gladde_slang = gladde_slang, hazelworm = hazelworm,
               levendbarende_hagedis = levendbarende_hagedis, muurhagedis = muurhagedis,
               ringslang = ringslang, zandhagedis = zandhagedis)
  return(cm_f)
}
# Wormerland:
cm_wormerland <- function(echt_f = echte_data_Wormerland, voorspeld_f = voorspelling_Wormerland)
{
  ringslang <- confusionMatrix(data = voorspeld_f$ringslang$aanwezigheid, 
                               reference = echt_f$ringslang$aanwezigheid,
                               positive = "1")
  cm_f <- list(ringslang = ringslang)
  return(cm_f)
}
# Zutphen:
cm_zutphen <- function(echt_f = echte_data_Zutphen, voorspeld_f = voorspelling_Zutphen)
{
  hazelworm <- confusionMatrix(data = voorspeld_f$hazelworm$aanwezigheid, 
                               reference = echt_f$hazelworm$aanwezigheid,
                               positive = "1")
  levendbarende_hagedis <- confusionMatrix(data = voorspeld_f$levendbarende_hagedis$aanwezigheid, 
                                           reference = echt_f$levendbarende_hagedis$aanwezigheid,
                                           positive = "1")
  ringslang <- confusionMatrix(data = voorspeld_f$ringslang$aanwezigheid, 
                               reference = echt_f$ringslang$aanwezigheid,
                               positive = "1")
  cm_f <- list(hazelworm = hazelworm, levendbarende_hagedis = levendbarende_hagedis, 
               ringslang = ringslang)
  return(cm_f)
}
# Extract the data from the confusion matrices and add them to a df
# Ede:
ophalen_info_cm_ede <- function(cm_f = confusion_matrices_Ede) {
  df_f1 <- data.frame(t(cm_f$overall))
  df_f2 <- data.frame(t(cm_f$byClass))
  df_f <- cbind(df_f1, df_f2)
  return(df_f)
}
# Epe:
ophalen_info_cm_epe <- function(cm_f = confusion_matrices_Epe) {
  df_f1 <- data.frame(t(cm_f$overall))
  df_f2 <- data.frame(t(cm_f$byClass))
  df_f <- cbind(df_f1, df_f2)
  return(df_f)
}
# Deurne:
ophalen_info_cm_deurne <- function(cm_f = confusion_matrices_Deurne) {
  df_f1 <- data.frame(t(cm_f$overall))
  df_f2 <- data.frame(t(cm_f$byClass))
  df_f <- cbind(df_f1, df_f2)
  return(df_f)
}
# Horst:
ophalen_info_cm_horst <- function(cm_f = confusion_matrices_Horst) {
  df_f1 <- data.frame(t(cm_f$overall))
  df_f2 <- data.frame(t(cm_f$byClass))
  df_f <- cbind(df_f1, df_f2)
  return(df_f)
}
# Huizen:
ophalen_info_cm_huizen <- function(cm_f = confusion_matrices_Huizen) {
  df_f1 <- data.frame(t(cm_f$overall))
  df_f2 <- data.frame(t(cm_f$byClass))
  df_f <- cbind(df_f1, df_f2)
  return(df_f)
}
# Ommen:
ophalen_info_cm_ommen <- function(cm_f = confusion_matrices_Ommen) {
  df_f1 <- data.frame(t(cm_f$overall))
  df_f2 <- data.frame(t(cm_f$byClass))
  df_f <- cbind(df_f1, df_f2)
  return(df_f)
}
# Roermond:
ophalen_info_cm_roermond <- function(cm_f = confusion_matrices_Roermond) {
  df_f1 <- data.frame(t(cm_f$overall))
  df_f2 <- data.frame(t(cm_f$byClass))
  df_f <- cbind(df_f1, df_f2)
  return(df_f)
}
# Westerveld:
ophalen_info_cm_westerveld <- function(cm_f = confusion_matrices_Westerveld) {
  df_f1 <- data.frame(t(cm_f$overall))
  df_f2 <- data.frame(t(cm_f$byClass))
  df_f <- cbind(df_f1, df_f2)
  return(df_f)
}
# Wormerland:
ophalen_info_cm_wormerland <- function(cm_f = confusion_matrices_Wormerland) {
  df_f1 <- data.frame(t(cm_f$overall))
  df_f2 <- data.frame(t(cm_f$byClass))
  df_f <- cbind(df_f1, df_f2)
  return(df_f)
}
# Zutphen:
ophalen_info_cm_zutphen <- function(cm_f = confusion_matrices_Zutphen) {
  df_f1 <- data.frame(t(cm_f$overall))
  df_f2 <- data.frame(t(cm_f$byClass))
  df_f <- cbind(df_f1, df_f2)
  return(df_f)
}