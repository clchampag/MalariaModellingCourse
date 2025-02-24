########################################################################################
######  Exploring Malaria Transmission Models in R with agent-based models     #########
########################################################################################
# 
# This exercise aims to explore different aspects of malaria transmission models using R. 
# It involves simulating the spread of malaria using individual-based models, 
# visualizing the results, and analyzing the impact of various parameters on the dynamics
# of the disease. 
#
########################################################################################

# Load required packages
library(ggplot2)
library(gridExtra)
library(grid)
set.seed(as.integer(Sys.time()))

population_size = 100
initial_infected = 1
transmission_probability = 0.1
recovery_probability = 0.05
total_steps = 100

# Function to initialize population
initialize_population <- function(population_size, initial_infected) {
  population <- data.frame(ID = 1:population_size, status = rep('S', population_size))
  population$status[1:initial_infected] <- 'I'
  return(population)
}

# Function to simulate disease spread
simulate_transmission <- function(population, transmission_probability, recovery_probability, total_steps) {
  n_infected <- numeric(total_steps + 1)
  n_infected[1] <- sum(population$status == 'I')
  
  totpop <- data.frame(t(population$status))
  
  for (step in 1:total_steps) {
    for (i in 1:nrow(population)) {
      if (population$status[i] == 'I') {
        
        for (j in 1:nrow(population)) {
          if (population$status[j] == 'S' && runif(1) < transmission_probability) {
            population$status[j] <- 'I'
          }
        }
        if (runif(1) < recovery_probability) {
          population$status[i] <- 'R'
        }
      }
    }
    n_infected[step + 1] <- sum(population$status == 'I')
    totpop=rbind(totpop,as.character(t(population$status)))
  }
  return(totpop)
}

population <- initialize_population(population_size, initial_infected)
simul <- simulate_transmission(population, transmission_probability, recovery_probability, total_steps)
View(simul)
