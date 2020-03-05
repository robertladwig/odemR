#' Solve ODE for simple odem model
#'
#' Extracts time (from date column) and space (aka depth) information
#'
#' @param bc matrix; Water temperatures (rows correspond to time, cols to depth)
#' @param params matrix; Water temperatures (rows correspond to time, cols to depth)
#' @param ini matrix; Water temperatures (rows correspond to time, cols to depth)
#' @param times vector having time information
#' @return results for oxgen concentrations
#' @export
#' @import deSolve 
#' @import LakeMetabolizer
run_model <- function(bc, params, ini, times){
  
  theta <- params[1]
  fsed_stratified = params[2]
  fsed_not_stratified  =  params[3]
  nep_stratified = params[4]
  nep_not_stratified = params[5]
  min_stratified = params[6]
  min_not_stratified = params[7]
  max_H = params[8]
  
  simple_oxygen <- function(t, y, params){
    
    K600 <- k.cole.base(wind(t))
    
    theta_total <- theta^(temp_total(t) - 20)
    
    NEP <- nep_not_stratified * volume_total(t) * theta_total
    
    MINER <- min_not_stratified * volume_total(t) * theta_total
    
    kO2 <- k600.2.kGAS.base(k600=K600, temperature = temp_total(t), gas='O2') # velocity value m/d?
    o2sat <- o2.at.sat.base(temp= temp_total(t), altitude = 300) * 1000 # mg O2/L 
    Fatm <- kO2 * (o2sat * volume_total(t) - y[1]) / max_H# mg/m * m/d = mg/d
    
    Fsed <- fsed_not_stratified * y[1] / max_H * theta_total # mg/m * m/d = mg/d
    
    if (is.na(therm_depth(t))){
      dOt <- NEP + Fatm + MINER - Fsed # units make sense bc every term is actually multiplied
      
      dOe <- (y[1] * volume_epi(t)) / volume_total(t) #/input.values$t.total[day]*input.values$vol_epil[day]
      
      dOh <- (y[1] * volume_hypo(t)) / volume_total(t)
    } else {
      
      theta_epil <- theta^(temp_epi(t) - 20)
      theta_hypo <-  theta^(temp_hypo(t) - 20)
      
      ##epil = epil.O2[i-1]+NEP[i]+Fatm[i]
      NEP_epil <- nep_stratified * volume_epi(t) * theta_epil
      
      kO2_epil <- k600.2.kGAS.base(k600 = K600, temperature = temp_epi(t), gas='O2')
      o2sat_epil<-o2.at.sat.base(temp = temp_epi(t), altitude = 300) * 1000
      Fatm_epil <- kO2_epil * (o2sat_epil * volume_epi(t) - y[2] ) / therm_depth(t)
      
      
      volumechange_epi = volume_epi(t) - volume_epi(t-1)  #in m^3
      volumechange_epi_proportion =  volumechange_epi /  volume_epi(t-1)
      Fepi <-  volumechange_epi_proportion * y[2]
      
      dOe <- NEP_epil + Fatm_epil + Fepi
      # print((o2_data[day,"o2_epil"]/input.values$vol_epil[day])/1000)
      
      ##hypo = hypo_o2[i-1]+Fhypo[i] - Fsed[i] 
      volumechange_hypo =  volume_hypo(t) -volume_hypo(t-1)  #in m^3
      volumechange_hypo_proportion =  volumechange_hypo / volume_hypo(t-1)
      Fhypo <- volumechange_hypo_proportion * y[3]
      
      MINER_hypo <- min_stratified * volume_hypo(t) * theta_hypo
      
      Fsed <- fsed_stratified * y[3]/(max_H - therm_depth(t) ) * theta_hypo
      
      dOh <- Fhypo - Fsed + MINER_hypo#+ NEP_hypo +Fatm_hypo
      
      dOt <- y[2] + y[3]
    }
    
    return(list(c(dOt, dOe, dOh)))
  }
  
  therm_depth <- approxfun(x = seq(from = 1, to = nrow(bc), by = 1), y = bc$td.depth, method = "linear", rule = 2)
  temp_total <- approxfun(x = seq(from = 1, to = nrow(bc), by = 1), y = bc$t.total, method = "linear", rule = 2)
  temp_epi <- approxfun(x = seq(from = 1, to = nrow(bc), by = 1), y = bc$t.epil, method = "linear", rule = 2)
  temp_hypo <- approxfun(x = seq(from = 1, to = nrow(bc), by = 1), y = bc$t.hypo, method = "linear", rule = 2)
  volume_total <- approxfun(x = seq(from = 1, to = nrow(bc), by = 1), y = bc$total_vol, method = "linear", rule = 2)
  volume_epi <- approxfun(x = seq(from = 1, to = nrow(bc), by = 1), y = bc$vol_epil, method = "linear", rule = 2)
  volume_hypo <- approxfun(x = seq(from = 1, to = nrow(bc), by = 1), y = bc$td.vol_hypodepth, method = "linear", rule = 2)
  wind <- approxfun(x = seq(from = 1, to = nrow(bc), by = 1), y = bc$wind, method = "linear", rule = 2)

  
  out <- ode(times = times, y = ini, func = simple_oxygen, parms = params, method = 'rk4')
  

}