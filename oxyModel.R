rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## packages
# get glmtools by first installing devtools, e.g. install.packages('devtools), then
# devtools::install_github('USGS-R/glmtools', ref='ggplot_overhaul')
library(devtools)
library(glmtools) 
library(dplyr)
library(ggplot2)
library(lubridate)
library(pracma)
library(readr)
library(LakeMetabolizer)
library(adagio)
library(zoo)
library(GenSA)

#devtools::install_github('LynetteGao/Limno_DataScience')
# library(simpleAnoxia)
library(odemR)


## load example data
lks <- list.dirs(path = 'inst/extdata/', full.names = TRUE, recursive = F)

ii = lks
  print(paste0('Running ',ii))
  data <- read.csv(paste0(ii,'/', list.files(ii, pattern = 'temperatures.csv', include.dirs = T)))
  meteo <- read.csv(paste0(ii,'/', list.files(ii, pattern = 'NLDAS', include.dirs = T)))
  
  chidx <- match(as.POSIXct(data$date),as.POSIXct(meteo$time))
  wind <- meteo$WindSpeed[chidx]
  
  if (length( list.files(ii, pattern = 'wq_data', include.dirs = T)) > 0){
  wq_data<- paste0(ii,'/', list.files(ii, pattern = 'wq_data', include.dirs = T))
 
  for(jj in wq_data){
      raw_obs <- read.csv(jj)

      if ('ActivityDepthHeighMeasure.MeasureValue' %in% colnames(raw_obs)){
        wq <- raw_obs %>%
          dplyr::filter(CharacteristicName== "Dissolved oxygen (DO)") %>%
          # dplyr::select(c('ActivityStartDate', 'ActivityStartTime.Time', 'ActivityDepthHeightMeasure.MeasureValue', 'ResultMeasureValue'))
          dplyr::select(c('ActivityStartDate', 'ActivityDepthHeighMeasure.MeasureValue', 'ResultMeasureValue'))
       wq <- rename(wq, 'ActivityDepthHeightMeasure.MeasureValue' = 'ActivityDepthHeighMeasure.MeasureValue')
      } else {
        wq <- raw_obs %>%
          dplyr::filter(CharacteristicName== "Dissolved oxygen (DO)") %>%
          # dplyr::select(c('ActivityStartDate', 'ActivityStartTime.Time', 'ActivityDepthHeightMeasure.MeasureValue', 'ResultMeasureValue'))
          dplyr::select(c('ActivityStartDate', 'ActivityDepthHeightMeasure.MeasureValue', 'ResultMeasureValue'))
        
      }
      
      
      obs<- wq
      # obs<-rbind(obs,more_obs)
  }
  obs$ActivityStartDate<-as.POSIXct(obs$ActivityStartDate)
  }

  if (is.factor(obs$ActivityDepthHeightMeasure.MeasureValue)){
    obs$ActivityDepthHeightMeasure.MeasureValue <-  as.numeric(as.character(obs$ActivityDepthHeightMeasure.MeasureValue))
  }
  if (is.factor(obs$ResultMeasureValue)){
    obs$ResultMeasureValue <-  as.numeric(as.character(obs$ResultMeasureValue))
  }
  
  # outlier detection
  outlier_values <- boxplot.stats(obs$ResultMeasureValue)$out 
  uvx <- match(outlier_values, obs$ResultMeasureValue)
  obs$ResultMeasureValue[uvx] <- NA
  
  eg_nml <- read_nml(paste0(ii,'/', list.files(ii, pattern = 'nml', include.dirs = T)))
  H <- abs(eg_nml$morphometry$H - max(eg_nml$morphometry$H))
  A <- eg_nml$morphometry$A

  # here's the promised input function (see R/helper.R), you can add the volume and
  # temperature conversion there
  input.values <- input(wtemp = data, H = H, A = A)
  
  input.values$wind <- wind
  
  proc.obs <- preprocess_obs(obs,input.values = input.values, H, A)
  
  theta = 1.08
  fsed_stratified = 0.01 *100
  fsed_not_stratified  =  0.0002
  nep_stratified = 0.1
  nep_not_stratified = 0
  min_stratified = 0.1
  min_not_stratified = 0
  
  parameters <- c(theta, fsed_stratified, fsed_not_stratified, nep_stratified, nep_not_stratified,
              min_stratified, min_not_stratified, max(H))

  system.time( results_oxygen <- run_model(bc = input.values, params = parameters,
                              ini = c(10, 8, 2), times = seq(from = 1, to = nrow(input.values), by = 1))
  )
  system.time(
  o2<- calc_do(input.values = input.values,fsed_stratified ,
               fsed_not_stratified,
               nep_stratified,
               nep_not_stratified,
               min_stratified ,
               min_not_stratified, wind)
  )
  
  plot(results_oxygen[,3] / input.values$vol_hypo / 1000)
  plot(o2$o2_hypo / input.values$vol_hypo / 1000)
  
  input.values$o2_epil <- o2[,"o2_epil"]
  input.values$o2_hypo <- o2[,"o2_hypo"]
  input.values$o2_total <- o2[,"o2_total"]

  input.values$year <- year(input.values$datetime)
  input.values$doy <- yday(input.values$datetime)

  
  test_data<-compare_predict_versus_observed(obs,input.values) 
  test_data$year <- year(test_data$day)
  test_data$doy <- yday(test_data$day)
  
  # fit = calc_rmse(test_data)
  fit = calc_fit(input.values , proc.obs )
  


  observed <- data.frame('time' = input.values$datetime[proc.obs[1,]],
                            'year' = input.values$year[proc.obs[1,]],
                            'doy' = input.values$doy[proc.obs[1,]],
                            'epi' =proc.obs[3,],
                            'hypo' = proc.obs[4,],
                         'epi_sim' = input.values$o2_epil[proc.obs[1,]]/input.values$vol_epil[proc.obs[1,]]/1000,
                         'hypo_sim' = input.values$o2_hypo[proc.obs[1,]]/input.values$vol_hypo[proc.obs[1,]]/1000)
  
  g1 <- ggplot(input.values) +
    geom_point(aes(doy, (o2_total/total_vol/1000), col = 'Total')) +
    geom_point(aes(doy, (o2_epil/vol_epil/1000), col = 'Epi')) +
    geom_point(aes(doy, (o2_hypo/vol_hypo/1000), col = 'Hypo')) +
    ylim(0,20)+
    facet_wrap(~year) +
    theme_bw() +
    geom_point(data = observed, aes(doy, epi, col = 'Obs_Epi'), size =2, alpha = 0.5) +
    geom_point(data = observed, aes(doy, hypo, col = 'Obs_Hypo'), size =2, alpha = 0.5)
  ggsave(file = paste0(ii,'/oxymodel.png'), g1, dpi=300, width=216,height=150,units='mm')
  