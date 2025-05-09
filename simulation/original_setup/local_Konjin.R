# This is a Rscript that focus on Konjin alternative simulation
library(mvtnorm)
library(ggplot2)  
library(dplyr)
library(tibble)
library(grid)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
set.seed(1)
n <- 300
B <- 200
source("simulation/Test_function.R")
# define the matrix function
signal_mat <- function(delta){
  A_d <- matrix(0, 6, 6)
  for (i in 1:3) {
    for (j in 1:3) {
      if (i == j){
        A_d[(2*(i-1)+1):(2*i), (2*(j-1)+1):(2*j)] <- (1 - delta) * diag(2)
      }else{
        A_d[(2*(i-1)+1):(2*i), (2*(j-1)+1):(2*j)] <- delta * diag(2)
      }
    }
  }
  return(A_d)
}

# data generation
generate_data <- function(type, sigma){
  
  Y <- matrix(0, nrow = n, ncol = 6)
  # loop over mvn, mvc and mvt
  switch (type,
    mvn = {
      Y_initial <- rmvnorm(3 * n, sigma = sigma)
    },
    mvc = {
      Y_initial <- rmvnorm(3 * n, sigma = sigma)^3
    },
    mvt = {
      Y_initial <- mvtnorm::rmvt(3 * n, sigma = sigma, df = 5)
    }
  )
  
  # output the matrix
  Y <- matrix(t(Y_initial), ncol = 6, byrow = TRUE)
  return(Y)
}

# generate reference data
emprical <- gensamdistrjdcov(n, dim_list = rep(2, 3), niter=1000)
sigma <- katlabutils::generate_cov_ar1(0.5, 2)
# loop over type
type_list <- c("mvn", "mvc", "mvt")
for (type in type_list) {
  
  # construct signal list and empty p-value array
  delta <- switch (type,
                   mvn = {
                     seq(0, 1.6, by = 0.4) / sqrt(n)
                   },
                   mvc = {
                     seq(0, 0.4, by = 0.1) / sqrt(n)
                   },
                   mvt = {
                     seq(0, 1, by = 0.25) / sqrt(n)
                   }
  )
  
  # construct empty p-value array
  p_value <- array(NA, dim = c(B, 2, length(delta)),
                   dimnames = list(
                     reps = 1:B,
                     method = c("RJdCov", "JdCov"),
                     signal = delta
                   ))
  
  # loop over delta
  for (signal in delta) {
    
    # construct the A_d matrix
    A_d <- signal_mat(delta = signal)
    
    for (k in 1:B) {
      
      # generate data
      Y <- generate_data(type = type, sigma = sigma)
      
      # Konijn alternative
      X <- Y%*%A_d
      X_list <- list(X[, 1:2], X[, 3:4], X[, 5:6])
      
      # compute RJdCov test
      statis <- computestatisticjdcov(X, dim_list = rep(2, 3))
      p_value[k, "RJdCov", as.character(signal)] <- (sum(emprical >= statis) + 1) / (length(emprical) + 1)
      
      # compute JdCov test
      statis <- jdcov.test(X_list, stat.type = "U", alpha = 0.05, B = 500)
      p_value[k, "JdCov", as.character(signal)] <- statis$p.value
      
      # print the p-value
      print(p_value[k,,as.character(signal)])
    }
  }
  
  # save the matrix
  saveRDS(p_value, sprintf("simulation/results/local_%s.rds", type))
  
}

########################## plot for local power ################################

TEXTWIDTH = 6.0689
TEXTHEIGHT = 9.33476
alpha <- 0.05

# read the results
local_power_mvn <- readRDS("simulation/results/local_mvn.rds")
local_power_mvc <- readRDS("simulation/results/local_mvc.rds")
local_power_mvt <- readRDS("simulation/results/local_mvt.rds")

# transform to tibble, compute the rejection rate and rescale the signal
rescaled_power_mvn <- as_tibble(as.table(local_power_mvn)) |>
  rename(p_value = n) |>
  group_by(signal, method) |>
  summarise(
    power = mean(p_value <= alpha)
  ) |>
  ungroup() |>
  mutate(signal = sqrt(n) * as.numeric(signal),
         a = 1)
rescaled_power_mvc <- as_tibble(as.table(local_power_mvc)) |>
  rename(p_value = n) |>
  group_by(signal, method) |>
  summarise(
    power = mean(p_value <= alpha)
  ) |>
  ungroup() |>
  mutate(signal = sqrt(n) * as.numeric(signal),
         a = 3)
rescaled_power_mvt <- as_tibble(as.table(local_power_mvt)) |>
  rename(p_value = n) |>
  group_by(signal, method) |>
  summarise(
    power = mean(p_value <= alpha)
  ) |>
  ungroup() |>
  mutate(signal = sqrt(n) * as.numeric(signal),
         a = 5)

# rbind the final results
local_power <- rbind(rescaled_power_mvn, rescaled_power_mvc, rescaled_power_mvt)

# create the parameter tibble
variable_parameters = bind_rows(
  dplyr::tibble(a = c(1, 3, 5), n = n,  
                variable_setting = paste0("a = ", a),
                fixed_setting = sprintf("n = %d", n))) %>%
  dplyr::mutate(setting = row_number(), 
                variable_setting = factor(variable_setting),
                fixed_setting = factor(fixed_setting))

# change levels
levels(variable_parameters$variable_setting) <- c("Multivariate Gaussian", "Gaussian copula", "Multivariate t-distribution") 

# append the parameter grid to the tibble
local_power <- local_power |>
  dplyr::left_join(variable_parameters, by = "a")

# Gaussian figure
power_plot_g = local_power |>
  dplyr::filter(a == 1) |>
  ggplot(aes(x = signal, y = power, color = method)) + 
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = seq(0,1,0.25), minor_breaks = seq(0, 1, 0.25), limits = c(0.02, NA)) +
  facet_wrap(variable_setting~., scales = "free") +
  geom_point(size = 1) +
  geom_line() +
  scale_color_manual(values = c("#90ee90", "red"))+
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) + 
  labs(colour = "test type")

# Gaussian Copula figure
power_plot_gc = local_power |>
  dplyr::filter(a == 3) |>
  ggplot(aes(x = signal, y = power, color = method)) + 
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(minor_breaks = seq(0, 1, 0.25), limits = c(0.02, NA))+
  facet_wrap(variable_setting~., scales = "free") +
  geom_point(size = 1) +
  geom_line() +
  scale_color_manual(values = c("#90ee90", "red"))+
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) + 
  labs(colour = "test type")

# Multivariate t-distribution Figure
power_plot_t = local_power |>
  dplyr::filter(a == 5) |>
  ggplot(aes(x = signal, y = power, color = method)) + 
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(minor_breaks = seq(0, 1, 0.25), limits = c(0.02, NA))+
  facet_wrap(variable_setting~., scales = "free") +
  geom_point(size = 1) +
  geom_line() +
  scale_color_manual(values = c("#90ee90", "red"))+
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) + 
  labs(colour = "test type")

# get the legend
legend <- get_plot_component(
  local_power |>
    dplyr::filter(a == 5) |>
    ggplot(aes(x = signal, y = power, color = method)) + 
    scale_x_continuous(expand = c(0.01,0.01)) +
    facet_wrap(variable_setting~., scales = "free") +
    geom_point(size = 1) +
    geom_line() +
    scale_color_manual(values = c("#90ee90", "red"))+
    theme_bw() + 
    theme(legend.position = "bottom",
          legend.key.size = unit(0.4, 'cm')) + 
    labs(color = "test type"),'guide-box', return_all = TRUE)

# combine plots together
plot <- plot_grid(power_plot_g, power_plot_gc, power_plot_t, ncol = 3, align = "v")

# create common x and y labels
y.grob <- textGrob("Power", 
                   gp=gpar(col="black"), rot=90)

x.grob <- textGrob("Signal strength", 
                   gp=gpar(col="black"))

g <- grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob, nrow=1),
                  plot_grid(legend[[3]], nrow = 1), nrow=2,heights=c(7, 1))


# save the plot
ggsave(
  filename = "simulation/new_figures/original/local_power_comparison.pdf",
  plot = g,
  device = "pdf",
  width = TEXTWIDTH,
  height = 0.25*TEXTHEIGHT
)

######################### inidividual plot #####################################

# MVN plot
pdf(file = "simulation/new_figures/original/local1.pdf",   # The directory you want to save the file in
    width =  0.7*TEXTWIDTH, # The width of the plot in inches
    height = 0.5*TEXTHEIGHT)

plot(unique(rescaled_power_mvn$signal), 
     rescaled_power_mvn |> filter(method == "RJdCov") |> dplyr::select(power) |> dplyr::pull(), 
     type = "o", pch = 19, col = "red",
     main="MVN", ylab="", xlab = "",ylim = c(0.05,1), yaxt="n",xaxt="n")
lines(unique(rescaled_power_mvn$signal), 
      rescaled_power_mvn |> filter(method == "JdCov") |> dplyr::select(power) |> dplyr::pull(), 
      type = "o", pch = 19, col = "blue")
legend("bottomright", legend=c("Proposed", "JdCov"),
       col=c("red", "blue"), lty = 1, cex=0.7)
axis(side = 2, at = c(seq(0, 1, 0.2)), cex.axis = 0.7, mgp = c(1,0.6,0))
axis(side = 1, at = unique(rescaled_power_mvn$signal), cex.axis = 0.7, mgp = c(1,0.6,0))
mtext("Signal strength", side=1, line = 1.6, cex=1.2)
mtext("Power", side=2, line = 1.6, cex=1.2)

dev.off()



# MVC plot
pdf(file = "simulation/new_figures/original/local2.pdf",   # The directory you want to save the file in
    width =  0.7*TEXTWIDTH, # The width of the plot in inches
    height = 0.5*TEXTHEIGHT)

plot(unique(rescaled_power_mvc$signal), 
     rescaled_power_mvc |> filter(method == "RJdCov") |> dplyr::select(power) |> dplyr::pull(), 
     type = "o", pch = 19, col = "red",
     main="MVN Copula", ylab="", xlab = "",ylim = c(0.05,1), yaxt="n",xaxt="n")
lines(unique(rescaled_power_mvc$signal), 
      rescaled_power_mvc |> filter(method == "JdCov") |> dplyr::select(power) |> dplyr::pull(), 
      type = "o", pch = 19, col = "blue")
legend("bottomright", legend=c("Proposed", "JdCov"),
       col=c("red", "blue"), lty = 1, cex=0.7)
axis(side = 2, at = c(seq(0, 1, 0.2)), cex.axis = 0.7, mgp = c(1,0.6,0))
axis(side = 1, at = unique(rescaled_power_mvc$signal), cex.axis = 0.7, mgp = c(1,0.6,0))
mtext("Signal strength", side=1, line = 1.6, cex=1.2)
mtext("Power", side=2, line = 1.6, cex=1.2)

dev.off()


# MVT plot
pdf(file = "simulation/new_figures/original/local3.pdf",   # The directory you want to save the file in
    width =  0.7*TEXTWIDTH, # The width of the plot in inches
    height = 0.5*TEXTHEIGHT)

plot(unique(rescaled_power_mvt$signal), 
     rescaled_power_mvt |> filter(method == "RJdCov") |> dplyr::select(power) |> dplyr::pull(), 
     type = "o", pch = 19, col = "red",
     main="MVT", ylab="", xlab = "",ylim = c(0.05,1), yaxt="n",xaxt="n")
lines(unique(rescaled_power_mvt$signal), 
      rescaled_power_mvt |> filter(method == "JdCov") |> dplyr::select(power) |> dplyr::pull(), 
      type = "o", pch = 19, col = "blue")
legend("bottomright", legend=c("Proposed", "JdCov"),
       col=c("red", "blue"), lty = 1, cex=0.7)
axis(side = 2, at = c(seq(0, 1, 0.2)), cex.axis = 0.7, mgp = c(1,0.6,0))
axis(side = 1, at = unique(rescaled_power_mvt$signal), cex.axis = 0.7, mgp = c(1,0.6,0))
mtext("Signal strength", side=1, line = 1.6, cex=1.2)
mtext("Power", side=2, line = 1.6, cex=1.2)

dev.off()






