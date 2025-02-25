modeldose <- function(modelname, stress, tspan, initial_states, param_names, param_values){
  pars.nostim = stats::setNames(as.list(param_values), param_names)
  inistate = initial_states
  out <-  deSolve::ode(
    y = c(S = stress, inistate),
    times = tspan,
    func = match.fun(modelname),
    parms = pars.nostim
  )
  return(out)
}

# This function simulates the ODE model with the selected parameters and returns the simulation.
run_model <- function(modelname, tspan, doses, ECs, initial_states, param_names, param_values){
  for (i in 1:length(doses)){
    xx <- modeldose(modelname, stress = ECs[i], tspan, initial_states, param_names, param_values)
    dose_uMadj = doses[i]
    xx = cbind(dose_uMadj, xx)
    
    if (i == 1) {
      modeldata = data.frame(xx)
    }
    else{
      modeldata = data.frame(rbind(modeldata, xx))
    }
  }
  modeldata = pivot_longer(modeldata, cols = 3:ncol(modeldata), names_to = "StateVar")
  return(modeldata)
}

compute_AIC <- function(measurements, statevars, AIC_K, modelname, doses, ECs, initial_states, param_names, param_values) {
  AIC_measurements = filter(measurements, !is.na(data4modelReal) & StateVar %in% statevars) %>%
    select(-'X')
  AIC_time = sort(unique(AIC_measurements$timeAfterExposure))

  AIC_modeldata = run_model(modelname, AIC_time, doses, ECs, initial_states, param_names, param_values) %>%
    filter(StateVar %in% statevars) %>%
    rename(timeAfterExposure = time)
    
  AIC_measurements = AIC_measurements %>%
    left_join(AIC_modeldata) %>%
    mutate(RS = (data4modelReal-value)^2)
  
  AIC_n = nrow(AIC_measurements)
  RSS = sum(AIC_measurements$RS)
  
  # print(AIC_n)
  
  if (AIC_K>AIC_n/40) {
    AIC = AIC_n*log(RSS/AIC_n) + (2*AIC_K**2 + 2*AIC_K)/(AIC_n-AIC_K-1)
    BIC = NA
    # print(paste("AICc = ", AIC))
  } else{
    AIC = AIC_n*log(RSS/AIC_n) + 2*AIC_K
    BIC = AIC_n*log(RSS/AIC_n) + log(AIC_n)*AIC_K
    # print(paste("AIC =", AIC))
    # print(paste("BIC =", BIC))
  }
  out = c(
    cost = RSS,
    AIC = AIC,
    BIC = BIC
  )
  return(out)
}
