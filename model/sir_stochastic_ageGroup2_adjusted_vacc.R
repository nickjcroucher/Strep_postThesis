freq <- user(1) # prev model is daily but aggregated to weekly
dt <- 1/freq
initial(time) <- 0
update(time) <- (step + 1) * dt

burnin_days <- 0

# 1. PARAMETERS ################################################################
time_shift_1 <- user(0, min = 0)
# trans_time_shift_1 <- 10^(time_shift_1)
beta_0 <- user(0, min = 0)
beta_1 <- user(0, min = 0, max = 1)
# theta <- 0.19 # proportion of vaccinated children in 0-14 age group
vaccIntro <- user() # 2648 days (April 2010) from Jan 2003 time point
theta <- (time-(burnin_days+vaccIntro))/(14*365)

UK_calibration_kids <- 1.07638532472038 # FIXED (Lochen et al., 2022)
UK_calibration_adults <- 0.536936186788821 # FIXED (Lochen et al., 2022)

# stratify log_delta
log_delta1 <- user(0, min = -10, max = 1)
# rho <- user(0, min = 0, max = 5)
log_delta2 <- user(0, min = -10, max = 1)

hypo_sigma_1_day <- 15.75 # (95% CI 7.88-31.49) (Chaguza et al., 2021)

# sigma_1 range as rate 0.06349206 (95% CI 0.03175611, 0.1269036)
sigma_1 <- 1/hypo_sigma_1_day # test sigma_1 (A -> R) later
# psi <- user(0, min = 0) # Immunity differences between children & adults
sigma_2 <- 1 # Assumed acute phase, 1 day

# Ageing
# group 1 = children, group 2 = adults
age_rate[1] <- 1 / (14 * 365)
age_rate[2] <- 0

# Natural mortality
# group 1 = children, group 2 = adults
mu_0[1] <- 0
mu_0[2] <- 1 / ((80.70 - 14) * 365)
# mu_0[1] <- 0
# mu_0[2] <- 0
# mu <- 0
mu_1 <- 0 # disease-related death, no data available
pi <- 3.141593 # FIXED
omega <- user(0, min = 0, max = 1)

# Dimensions of arrays #########################################################
N_age <- 2

dim(N_ini) <- N_age
# dim(S_ini) <- N_age
dim(A_ini) <- N_age # consider adults only, kids as faction
# dim(log_A_ini) <- N_age

dim(N) <- N_age
dim(S) <- N_age
dim(A) <- N_age
dim(D) <- N_age
dim(R) <- N_age

dim(m) <- c(N_age, N_age)
dim(foi_ij) <- c(N_age, N_age)
# dim(vacc_m) <- c(N_age, N_age)
# dim(vacc) <- N_age
dim(v) <- N_age
dim(lambda) <- N_age
dim(delta) <- N_age
dim(mu_0) <- N_age
dim(age_rate) <- N_age

dim(p_Suscep) <- N_age
dim(p_Asym) <- N_age
dim(p_Dis) <- N_age
dim(p_Rec) <- N_age

dim(p_SA) <- N_age
dim(p_SR) <- N_age
dim(p_AD) <- N_age
dim(p_AR) <- N_age
dim(p_DR) <- N_age
dim(p_Dd) <- N_age
dim(p_RS) <- N_age

dim(n_Sborn) <- N_age
dim(n_Suscep) <- N_age
dim(n_SA) <- N_age
dim(n_SR) <- N_age
dim(n_Sdead) <- N_age
dim(n_Asym) <- N_age
dim(n_AD) <- N_age
dim(n_AR) <- N_age
dim(n_Adead) <- N_age
dim(n_Dis) <- N_age
dim(n_DR) <- N_age
dim(n_Dd) <- N_age
dim(n_Ddead) <- N_age
dim(n_Resist) <- N_age
dim(n_RS) <- N_age
dim(n_Rdead) <- N_age
dim(n_age_S) <- N_age
dim(n_age_A) <- N_age
dim(n_age_D) <- N_age
dim(n_age_R) <- N_age

# 2. INITIAL VALUES ############################################################
# Initial values (user-defined parameters)
N_ini[] <- user()
max_A_ini <- 0
min_A_ini <- log(1/(N_ini[1]+N_ini[2]), 10)

# directly test log_A_ini as scaled
log_A_ini <- user()
phi <- user(0, min = 0, max = 3)
A_ini[1] <- (if (N_ini[1] - as.integer(10^(log_A_ini*(max_A_ini-min_A_ini)+min_A_ini)*N_ini[1]) <= 0) 0 else
  (as.integer(10^(log_A_ini*(max_A_ini-min_A_ini)+min_A_ini)*N_ini[1])))

A_ini[2] <- (if (N_ini[2] - as.integer(10^(log_A_ini*(max_A_ini-min_A_ini)+min_A_ini)*phi*N_ini[2]) <= 0) 0 else
  (as.integer(10^(log_A_ini*(max_A_ini-min_A_ini)+min_A_ini)*phi*N_ini[2])))

# Age-structured states:
initial(S[]) <- N_ini[i] -(A_ini[i]+0+0) # D_ini = R_ini = 0
initial(A[]) <- A_ini[i]
initial(D[]) <- 0
initial(R[]) <- 0

# Initial states:
initial(N_tot) <- sum(N_ini)
initial(S_tot) <- sum(N_ini) -(sum(A_ini)+0+0)
initial(A_tot) <- sum(A_ini)
initial(D_tot) <- 0
initial(R_tot) <- 0

# make it traditional way:
initial(n_AD1_weekly) <- 0
initial(n_AD2_weekly) <- 0

# 3. UPDATES ###################################################################
# age-structured contact matrix featured in lambda:
# https://mrc-ide.github.io/odin.dust/articles/sir_models.html
N[] <- S[i] + A[i] + D[i] + R[i]

m[, ] <- user() # age-structured contact matrix

# coverage*efficacy*proportion of kids 2y.o. (from 0-14)
# vacc_m[1, 1] <- 0 #0.9*0.862*theta # child->child
# vacc_m[1, 2] <- 0 #0.9*0.862*theta  # adult->child
# vacc_m[2, 1] <- 0
# vacc_m[2, 2] <- 0

# vacc_eff <- user(0, min = 0, max = 1) # previously 0.862
# vacc_m[1, ] <- 0.9*0.862*theta # child->child & adult -> child
# vacc_m[2, ] <- 0

# vacc must be defined as coverage*efficacy*proportion of kids 2y.o.*theta (theta as gradual vaccination)
# 0.9*efficacy*0.19*theta
vacc <- user(0, min = 0, max = 1)
v[1] <- (if (time >= (burnin_days+vaccIntro)*freq)
  vacc*theta
  else
    vacc*0
  )
v[2] <- vacc*0

# additional time steps for beta_1 (2 years)
# difractions based on PCV7 era 
beta_diff <- 1 #user(0, min = 0, max = 1)

# beta_diff for pre-PCV7 only
beta <- (if (time < 0) beta_0 else 
  (if (time >= (burnin_days+(1339))*freq) # time <= (burnin_days+1461)*freq && 
  (beta_0*((1+beta_1*cos(2*pi*((time_shift_1*(365))+time)/(365))))) else
    (beta_0*((1+beta_1*beta_diff*cos(2*pi*((time_shift_1*(365))+time)/(365)))))))

# foi_ij[, ] <- (if (time >= (burnin_days+vaccIntro)*freq)
#   beta * m[i, j] * (((A[j] + D[j])/N[j]) * (1 - vacc_m[i, j]))
#   else
#     beta * m[i, j] * (((A[j] + D[j])/N[j]))
# )

foi_ij[, ] <- beta * m[i, j] * (((A[j] + D[j])/N[j]))

lambda[] <- sum(foi_ij[i, ])

delta[1] <- (10^(log_delta1))*UK_calibration_kids
delta[2] <- (10^(log_delta2))*UK_calibration_adults
# delta[2] <- (10^(log_delta1))*rho #*UK_calibration_adults

# sigma_1[1] <- hypo_sigma_1 # test no A -> R in kids
# sigma_1[2] <- psi*hypo_sigma_1

# Cumulative hazard
p_Suscep[] <- lambda[i]+mu_0[i]+age_rate[i]+v[i]
p_Asym[] <- delta[i]+sigma_1+mu_0[i]+age_rate[i]
p_Dis[] <- sigma_2+mu_1+mu_0[i]+age_rate[i]
p_Rec[] <- omega+mu_0[i]+age_rate[i]

p_SA[] <- (if (lambda[i]/(lambda[i]+mu_0[i]+age_rate[i]+v[i]) <= 0) 0 else 
  (lambda[i]/(lambda[i]+mu_0[i]+age_rate[i]+v[i])))
p_SR[] <- (if (v[i]/(mu_0[i]+age_rate[i]+v[i]) <= 0) 0 else 
  (v[i]/(mu_0[i]+age_rate[i]+v[i])))

p_AD[] <- (if (delta[i]/(delta[i]+sigma_1+mu_0[i]+age_rate[i]) <= 0) 0 else 
  (delta[i]/(delta[i]+sigma_1+mu_0[i]+age_rate[i])))
p_AR[] <- (if (sigma_1/(sigma_1+mu_0[i]+age_rate[i]) <= 0) 0 else 
  (sigma_1/(sigma_1+mu_0[i]+age_rate[i])))

p_DR[] <- (if (sigma_2/(sigma_2+mu_1+mu_0[i]+age_rate[i]) <= 0) 0 else 
  (sigma_2/(sigma_2+mu_1+mu_0[i]+age_rate[i])))
p_Dd[] <- (if (mu_1/(mu_1+mu_0[i]+age_rate[i]) <= 0) 0 else 
  (mu_1/(mu_1+mu_0[i]+age_rate[i])))

p_RS[] <- (if (omega/(omega+mu_0[i]+age_rate[i]) <= 0) 0 else 
  (omega/(omega+mu_0[i]+age_rate[i])))


# Draws for numbers changing between compartments
# Leaving S
n_Suscep[] <- rbinom(S[i], 1 - exp(-p_Suscep[i]*dt))
n_SA[] <- rbinom(n_Suscep[i], p_SA[i])
n_SR[] <- rbinom((n_Suscep[i] - n_SA[i]), p_SR[i])
n_age_S[1] <- if (i == 1) rbinom((n_Suscep[i] - n_SA[i] - n_SR[i]),
                                 age_rate[i]/(mu_0[i] + age_rate[i])) else 0
n_age_S[2] <- 0
n_Sdead[] <- n_Suscep[i] - n_SA[i] - n_SR[i] - n_age_S[i]

# Leaving A
n_Asym[] <- rbinom(A[i], 1- exp(-p_Asym[i]*dt))
n_AD[] <- rbinom(n_Asym[i], p_AD[i])
n_AR[] <- rbinom((n_Asym[i] - n_AD[i]), p_AR[i])
n_age_A[1] <- if (i == 1) rbinom((n_Asym[i] - n_AD[i] - n_AR[i]),
                                 age_rate[i]/(mu_0[i] + age_rate[i])) else 0
n_age_A[2] <- 0
n_Adead[] <- n_Asym[i] - n_AD[i] - n_AR[i] - n_age_A[i]

# Leaving D
n_Dis[] <- rbinom(D[i], 1- exp(-p_Dis[i]*dt))
n_DR[] <- rbinom(n_Dis[i], p_DR[i])
n_Dd[] <- rbinom((n_Dis[i] - n_DR[i]), p_Dd[i])
n_age_D[1] <- if (i == 1) rbinom(n_Dis[i] - n_DR[i] - n_Dd[i],
                                 age_rate[i]/(mu_0[i] + age_rate[i])) else 0
n_age_D[2] <- 0
n_Ddead[] <- n_Dis[i] - n_DR[i] - n_Dd[i] - n_age_D[i]

# Leaving R
n_Resist[] <- rbinom(R[i], 1- exp(-p_Rec[i]*dt)) # RS is considered 0 in both age groups
n_RS[] <- rbinom(n_Resist[i], p_RS[i])
n_age_R[1] <- if (i == 1) rbinom(n_Resist[i] - n_RS[i],
                                 age_rate[i]/(mu_0[i] + age_rate[i])) else 0
n_age_R[2] <- 0
n_Rdead[] <- n_Resist[i] - n_RS[i] - n_age_R[i]

# Equations for transitions between compartments by age group
n_Sborn[] <- n_Sdead[i] + n_Adead[i] + n_Dd[i] + n_Ddead[i] + n_Rdead[i]
born <- sum(n_Sborn)

update(S[1]) <- S[1] + (born + n_RS[1]) - (n_SA[1] + n_Sdead[1] + n_age_S[1] + n_SR[1])
update(S[2]) <- S[2] + (n_age_S[1] + n_RS[2]) - (n_SA[2] + n_Sdead[2] + n_age_S[2] + n_SR[2])

update(A[1]) <- A[1] + n_SA[1] - (n_AD[1] + n_AR[1] + n_Adead[1] + n_age_A[1])
update(A[2]) <- A[2] + n_SA[2] + n_age_A[1] - (n_AD[2] + n_AR[2] + n_Adead[2])

update(D[1]) <- D[1] + n_AD[1] - (n_DR[1] + n_Dd[1] + n_Ddead[1] + n_age_D[1])
update(D[2]) <- D[2] + n_AD[2] + n_age_D[1] - (n_DR[2] + n_Dd[2] + n_Ddead[2])

update(R[1]) <- R[1] + (n_AR[1] + n_DR[1] + n_SR[1]) - (n_RS[1] + n_Rdead[1] + n_age_R[1])
update(R[2]) <- R[2] + (n_AR[2] + n_DR[2] + n_age_R[1]+ n_SR[2]) - (n_RS[2] + n_Rdead[2])

# Core equations of the transitions
update(N_tot) <- sum(N)
update(S_tot) <- sum(S)
update(A_tot) <- sum(A)
update(D_tot) <- sum(D)
update(R_tot) <- sum(R)
# based on tutorial: https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model

# that "little trick" previously explained in https://github.com/mrc-ide/dust/blob/master/src/sir.cpp for cumulative incidence:
# based on tutorial: https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model
update(n_AD1_weekly) <- if (step %% (7*freq) == 0) n_AD[1] else n_AD1_weekly + n_AD[1]
update(n_AD2_weekly) <- if (step %% (7*freq) == 0) n_AD[2] else n_AD2_weekly + n_AD[2]

