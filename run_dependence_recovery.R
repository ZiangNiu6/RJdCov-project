# This is a script running the ICA simulation
source("dependence_recovery/Test_function.R")

########################### run the simulations ################################
# Gaussian simulation
source("dependence_recovery/Gaussian-simulation.R")

# Student simulation
source("dependence_recovery/Student-simulation.R")

########################### summarize and plot #################################
# summarize the result
source("dependence_recovery/Gaussian_result.R")
source("dependence_recovery/Student-result.R")

# plot (Figure 8 in the supplement)
source("dependence_recovery/plot.R")