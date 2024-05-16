#Domaća zadaća - Klobucharov model, Jan Šubarić

rm(list=ls())

#korisne funkcije
convert_UTC_to_seconds <- function(utc_time) {
  posix_time <- as.POSIXct(utc_time, tz = "UTC")
  seconds_since_epoch <- as.numeric(posix_time)
  return(seconds_since_epoch)
}

deg2rad <- function(degrees) {
  return(degrees * pi / 180)
}

calculate_time <- function(lambda, t_user) {
  t <- 43200 * lambda / pi + t_user
  while (t > 86400) {
    t <- t - 86400
  }
  
  while (t < 0) {
    t <- t + 86400
  }
  return(t)
}

calculate_amplitude <- function(phi_m, alpha) {
  amplitude <- 0
  for (n in 0:3) {
    amplitude <- amplitude + alpha[n+1] * ((phi_m / pi) ^ n)
  }
  
  if (amplitude < 0){
    return(0)
  } 
  
  return (amplitude)
}

calculate_period <- function(phi_m, beta) {
  period <- 0
  for (n in 0:3) {
    period <- period + beta[n + 1] * ((phi_m / pi) ^ n)
  }
  
  if (period > 72000) {
    return(72000)
  }
  
  return(period)
}

calculate_delay <- function(A, X, F_i) {
  if(abs(X) < pi/2) {
    I = (5e-09 + A * cos(X)) * F_i
  }else{
    I = 5e-09 * F_i
  }
  
  return(I) 
}

### UNOS PODATAKA
variable_names <- c("t_gps", "E", "A", "phi_u", "lambda_u", "alpha", "beta")
values <- list()
for (var in variable_names) {
  if (var == "t_gps"){
    input <- readline(paste("Unesite", var, "vrijednost (UTC format 'YYYY-MM-DD HH:MM:SS'): ", sep = " "))
  } else if (var %in% c("alpha", "beta")){
    input <- readline(paste("Unesite", var, "vrijednosti odvojene razmakom: ", sep = " "))
  }
  else{
    input <- readline(paste("Unesite", var, "vrijednost (deg): ", sep = " "))
  }
  
  if (var %in% c("alpha", "beta")) {
    values[[var]] <- as.numeric(strsplit(input, " ")[[1]])
  } else if (var == "t_gps") {
    values[[var]] <- convert_UTC_to_seconds(input)
  } else {
    values[[var]] <- deg2rad(as.numeric(input))
  }
}


#definirane konstante
RE <- 6378
h <- 350
c = 2.99792458e08
phi_p <- deg2rad(78.3)
lambda_p <- deg2rad(291)

### Klobucharov model
klobuchar_model <- function(phi_u, phi_p, lambda_user, E, A, alpha, beta, t_user)
{
  # Izračun kuta PSI 
  psi <- pi / 2 - E - asin((RE / (RE + h)) * cos(E))
  # Izračun geografske širine točke IPP
  phi <- asin(sin(phi_u) * cos(psi) + cos(phi_u) * sin(psi) * cos(A))
  # Izračun geografske dužine IPP
  lambda <- lambda_user + (psi * sin(A)) / cos(phi);
  # Izračun geomagnetske širine
  phi_m <- asin(sin(phi) * sin(phi_p) + cos(phi) * cos(phi_p) * cos(lambda - lambda_p))
  # Izračun lokalnog vremena u IPP
  t <- calculate_time(lambda, t_user)
  # Izračun amplitude ionosferskog kašnjenja
  A <- calculate_amplitude(phi_m, alpha)
  # Izračun perioda ionosferskog kašnjenja
  P <- calculate_period(phi_m, beta)
  # Izračun faze ionosferskog kašnjenja
  X <- 2 * pi * (t - 50400) / P
  # Izračun ionosferske funkcije preslikavanja
  F_i <- (1 - ((RE / (RE + h)) * cos(E)) ^ 2) ^ -0.5
  # Izračun vertikalnog ionosferskog kašnjenja u sekundama
  I <- calculate_delay(A, X, F_i)
  
  return(I);
}

#ispis jednog rezultata za numerički primjer
ion_delay <- klobuchar_model(values$phi_u, phi_p, values$lambda_u, values$E, values$A, values$alpha, values$beta, values$t_gps) * c
cat("Ionosfersko kašnjenje numeričkog primjera iznosi (dION):", ion_delay)

### SIMULACIJA 24-SATNOG IONOSFERSKOG KAŠNJENJA U KORACIMA OD 1 MINUTE

#vektora za pohranu kašnjenja za svaku minutu kroz 24 sata
delay_sim_minute <- numeric(24*60)

for (hour in 0:23) {
  for (minute in 0:59) {
    #izračunavanje GPS vremena za trenutnu minutu
    t_gps_minute <- values$t_gps + hour * 3600 + minute * 60
    
    #izračun kašnjenja za trenutnu minutu
    delay_sim_minute[hour * 60 + minute + 1] <- klobuchar_model(values$phi_u, phi_p, values$lambda_u, values$E, values$A, values$alpha, values$beta, t_gps_minute) * c
  }
}

#ispis rezultata simulacije
print(delay_sim_minute)
print(length(delay_sim_minute))