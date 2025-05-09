TEXTWIDTH = 6.0689
TEXTHEIGHT = 9.33476

# create figure folder
figure_dir <- "ICA/figure"
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir)
  cat("Directory created:", figure_dir, "\n")
} else {
  cat("Directory already exists:", figure_dir, "\n")
}


# load the result
library(tibble)
library(dplyr)
library(purrr)
library(tidyr)
result_ica <- readRDS("ICA/results/ICA_results.rds")

# extend the results
data_to_plot <- result_ica |> mutate(
  grid_id = 1:n(),
  result = map(result, ~ as_tibble(.))
) %>%
  unnest_wider(result, names_sep = "_") |> 
  unnest(c(result_proposed, result_MT)) |>
  pivot_longer(cols = c(result_proposed, result_MT),
               values_to = "ICA_results",
               names_to = "method") |>
  dplyr::select(grid_id, method, ICA_results)

# print the result
summarized_results <- data_to_plot |> dplyr::group_by(grid_id, method) |> 
  dplyr::summarise(mean = mean(ICA_results), sd = sd(ICA_results)) |> 
  dplyr::ungroup() |> 
  tidyr::pivot_wider(id_cols = grid_id, names_from = method, values_from = c(mean, sd)) |> 
  as.data.frame() |> 
  round(3)

summarized_results[, c("mean_result_proposed", "sd_result_proposed", 
                       "mean_result_MT", "sd_result_MT")]
  

## plot code
pdf(file = sprintf("%s/ica_box_plot1.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)
name_vec <- rep(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"),
                each = 2)

boxplot(ICA_results~method*grid_id, 
        data = as.data.frame(data_to_plot |> filter(grid_id %in% 1:6)), 
        notch=FALSE, xlim = c(0.8, 12.1),
        outline=FALSE,
        col=(c("blue","darkgreen")),
        ylab="", xlab = "", yaxt="n",xaxt="n",
        main="ICA estimation error")
axis(side = 2, at=c(seq(0, 0.6, 0.1)), cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=c(seq(1, 12, 1)), labels = name_vec[1:12],cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Distributions", side=1, line = 1.6, cex=0.7)
mtext("Estimation error", side=2, line = 1.6, cex=0.7)


dev.off()

pdf(file = sprintf("%s/ica_box_plot2.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)


boxplot(ICA_results~method*grid_id, 
        data = as.data.frame(data_to_plot |> filter(grid_id %in% 7:12)), 
        notch=FALSE,
        outline=FALSE,
        col=(c("blue","darkgreen")), xlim = c(0.8, 12.1),
        ylab="", xlab = "", yaxt="n",xaxt="n",
        main="ICA estimation error")
axis(side = 2, at=c(seq(0, 0.5, 0.1)),cex.axis=0.7, mgp=c(1,0.5,0))
axis(side = 1, at=c(seq(1, 12, 1)), labels = name_vec[13:24],cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Distributions", side=1, line = 1.6, cex=0.7)
mtext("Estimation error", side=2, line = 1.6, cex=0.7)


dev.off()
