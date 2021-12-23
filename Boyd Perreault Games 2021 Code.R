# Code for Boyd, R. and C. Perreault, 2022. Evolution of social learning with payoff and content bias, Games. 
#Last Modified: December 23rd 2021. 

##################################################################################################
####---General structure of the code (see article for more complete description of the code)---###
##################################################################################################
#Each generation, the following happens:
# 1. The state of the environment switches randomly
# 2. The individuals get a social cue, an environmental cue, and a payoff cue:
# 3. Individuals combine those three cues using the optimal decision rule analytically derived above to chose a behavior.
# 4. Viability selection occurs. The baseline fitness is W. Individuals with the favored behavior, given the current state of the environment, get a fitness benefit d. Reproduction is based on relative fitness, which is calculated by dividing individual fitness by the maximum fitness in the population. The fitness of individuals relative to the maximum fitness in the population is used as a vector of weights in the sample function in the R language base package in order to sample the individuals that reproduce and transmit their alleles to the next generation.
# 5. Mutations in G and g alleles in the next generation occur with probability M. The values of the mutant alleles are drawn from a normal distribution with mean equal the allele of the parent and standard deviation msd. G and g are unlinked and mutate independently.

######################
###---Parameters---###
######################

pop<- 10000 # Population size
n<- c(12)   # Number of social models observed by agents
env_cue_quality<-c(0.01, 0.001) # Mean of the distribution of environmental cue
sd_env_cue <- 1 # Standard deviation of the environmental cue distribution. 
env_change<- c(0.1, 0.01) # Probability that the state of the environment changes. Here we treat environmental state 0 as the environmental state 2 in the paper. Behavior 1 is favored in env 1 and behavior 0 is favored in env 0
var_payoff_cue <- c(25) # Variance in the payoff cue
m <- 0.05 # Mutation rate. Probability that an individual's genes mutate at a given time step
m_sd <- 0.5 # Magnitude of mutation. Mutation is modeled as a normal distribution with mean = current genotype (G or g value) and standard deviation = m_sd
f <- 0.5 # Fitness benefit associated with favored behavior

running_average_window<- 5000 # Specifies how many environmental changes over which moving average of g or G is calculated. For instance, running_average_window 100 means 100 env shifts. 
history_slope_fitting<- 2000 # The simulation is stopped by fitting a lm to median(G/g)~time. The parameter history_slope_fitting specifies how many previous generations are considered when fitting the model
tape_median_G<-0 # Used to record median values of G in the population 
tape_median_g<-0 # Used to record median values of g in the population
tape_average_fitness<-0

n_sims<-length(n)*length(env_change)*length(env_cue_quality)*length(var_payoff_cue) # How many simulations are to be run for a particular set of parameters
results<- mat.or.vec(n_sims,7)
colnames(results)<-c("n_models", "env_cue_quality", "rate_env_change", "variance_payoff","G","g", "mean_fitness")

sim<- 1 # Keeps track of number of simulations (unique set of parameters) have been run so far

######################
###---SIMULATION---###
######################

for (n_models in n){
  for (mu_env_cue in env_cue_quality) {
    for (prob_env_change in env_change) {
      for (variance_payoff_cue in var_payoff_cue) {
        
        ###########################
        ###---Setting objects---###
        ###########################
        env_state <- 1 # The current state of the environment. Starts at 1. Varies between 1 and 0.
        d <- -1 # The expected payoff of behavior 1 in Env 1 is mu + d. The payoff cue is normally distributed with mean -d in env 1 and with mean d in env 0
        g <- rep(1,pop) # Small g is the genotype that specifies how environmental information is weighted
        G <- rep(1,pop) # Big G is the genotype that specifies how payoff information is weighted
        behavior <-rep(1,pop) # The behavior of each individuals in the population. Everyone starts with behavior 1
        observation_window<-running_average_window/prob_env_change
        
        ##############################
        ###---Run the simulation---###
        ##############################
        stop = 0
        tick<- 1
        while (stop==0){  # Keep the simulation running until equilibrium values of G and g have been found
          
          ################################
          ###---ENVIRONMENTAL CHANGE---###
          ################################
          if (runif(1) <= prob_env_change) # The environment switches state with probability prob_env_change 
          {env_state <- +!env_state} & {mu_env_cue <- mu_env_cue*-1} & {d<- d*-1} # If environment changes, then change (1) the state of the environment; (2) the mean of the environmental cue distribution mu_env_cue; (3) the mean of payoff cue d
          
          ###########################################
          ###---AGENTS GET AN ENVIRONMENTAL CUE---###
          ###########################################
          env_cues <- rnorm(pop, mean = mu_env_cue, sd = sd_env_cue) # An array of environmental cue. One cue per agent
          
          ###################################
          ###---AGENTS GET A SOCIAL CUE---###
          ###################################
          social_cues <- rbinom(pop, size = n_models, prob = sum(behavior)/pop) # An array containing the number of models with behavior 1 observed by each individual
          
          ###########################################
          ###---AGENTS GET A PAYOFF CUE (x2-x1)---###
          ###########################################
          variance_denominator_a <- 1/social_cues                   # The first side of the additive formula for variance of the payoff cue (1/j), i.e. variance = v(1/j + 1/n-j)
          variance_denominator_a[variance_denominator_a==Inf] <-0   # Avoid Inf value when no models have behavior 1 (j == 0). Replaces Inf by 0
          variance_denominator_b <- 1/(n_models-social_cues)        # One side of the denominator of variance for payoff cue (1/n-j). Calculated separately because it can generate Inf value when j==n
          variance_denominator_b[variance_denominator_b==Inf] <-0   # Avoid Inf value when all models observed show behavior 1 (i.e. j == n)
          var_payoff <- variance_payoff_cue*(variance_denominator_a + variance_denominator_b) # A matrix of variance for the payoff cue. Payoff cue is normally distributed mean d/-d and variance nv/j(n-j)
          payoff_cues <- rnorm(pop, mean = d, sd = sqrt(var_payoff)) # Creates an array of payoff cue. One per individual 
          
          ######################################
          ###---AGENT DECIDE ON A BEHAVIOR---### 
          ######################################
          var_models <- (social_cues*(n_models-social_cues))/social_cues  # This is the part of the decision rule that is mulmedian(Gtiplied by G and (x2-x1). Calculated first to eliminate Inf that occur when j = 0
          var_models[is.nan(var_models)]<-0                               # Replace Inf values when j = 0 with zeros 
          behavior <- social_cues - n_models/2 > (((G * var_models * payoff_cues) - g*env_cues)/2) # Creates a array with the chosen behavior (TRUE ==1, FALSE ==0) 
          
          #####################
          ###---SELECTION---###
          #####################
          q <- behavior==env_state # Array that specifies whether individuals have favored (True) or unfavored (False) behavior
          w <-(q+f)/max(q+f) # w is the relative fitness of individuals
          G <- sample(G, size = pop, replace = TRUE, prob = w)  # Next generation of G. This replaces the array of G with new values picked randomly from the previous generation with probability w (i.e. relative fitness)
          g <- sample(g, size = pop, replace = TRUE, prob = w)  # Next generation of g. This replaces the array of g with new values picked randomly from the previous generation with probability w (i.e. relative fitness)
          
          ####################
          ###---MUTATION---###
          ####################
          mutate_G<-rbinom(pop,1,m) # An array that specifies whether an agent's G mutates. 
          mutate_g<-rbinom(pop,1,m) # An array that specifies whether an agents's g mutates
          G[mutate_G==1]<-rnorm(length(G[mutate_G==1]), mean =G[mutate_G==1], sd = m_sd) # Updates genotype. Mutation is model as normal distribution with mean = current G value and sd = m_sd
          g[mutate_g==1]<-rnorm(length(g[mutate_g==1]), mean =g[mutate_g==1], sd = m_sd) # Updates genotype. Mutation is model as normal distribution with mean = current g value and sd = m_sd
          
          ##################################
          ###---IS SIMULATION STOPPED?---###
          ##################################
      
          # Save average fitness, median G and median g on rolling tapes
          tape_average_fitness[history_slope_fitting+1]<-mean(w) # Store average fitness on a rolling recording tape. Append new value at the end of the tape
          tape_average_fitness<-tape_average_fitness[-1] 
          tape_median_G[history_slope_fitting+1]<-median(G) 
          tape_median_g[history_slope_fitting+1]<-median(g)
          tape_median_G <- tape_median_G[-1] 
          tape_median_g <- tape_median_g[-1] 
          
          if (tick>observation_window){  # If enough time steps have elapsed 
            lm_G <-lm(tape_median_G~c(seq(from = 1, to = history_slope_fitting, by = 1))) # Fit a linear model to the median values on tape. Median ~ Time
            lm_g <-lm(tape_median_g~c(seq(from = 1, to = history_slope_fitting, by = 1))) # Fit a linear model to the median values on tape. Median ~ Time
            slope_G <- summary(lm_G)$coefficients[2,1] # Extract slope for lm object
            slope_g <- summary(lm_g)$coefficients[2,1] # Extract slope for lm object
            if (slope_G < 0.001 &&  slope_g< 0.001) {stop<-1} # If median values of G and g have been stable enough (small slope), stop the simulation
          }
          tick<- tick+1}
        
#########################
###---Store Results---###
#########################
results[sim,1]<-n_models
results[sim,2]<-mu_env_cue
results[sim,3]<-prob_env_change
results[sim,4]<-variance_payoff_cue
results[sim,5]<-median(G)
results[sim,6]<-median(g)
results[sim,7]<-mean(tape_average_fitness)
sim<-sim+1
}}}} 
