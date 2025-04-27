
rm(list=ls())

source(paste(getwd(),"/src/environment.R",sep=""))

data_MSE <- read.table(file=paste(data_dir,"/sim2_MSE_data.txt",sep=""))
data_coeff <- read_rds(file=paste(data_dir,"/sim2_coeff_data.rds",sep=""))
param <- read_rds(file=paste(data_dir,"/sim2_param.rds",sep=""))

attach(data_coeff)
attach(param)

data_MSE <- data_MSE %>%
  gather(key = "Method", value = "MSE")


method_mapping <- c("naive_MSE" = "Naive Forecast",
                    "pcr_MSE" = "PCR Forecast",
                    "pls_MSE" = "PLS Forecast",
                    "spcr_MSE" = "SPCR Forecast")

data_MSE <- data_MSE %>%
  mutate(Method = method_mapping[Method])

method_order <- c("Naive Forecast", "PCR Forecast", "PLS Forecast", "SPCR Forecast")

data_MSE$Method <- factor(data_MSE$Method, levels = method_order)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#method_colors <- c("Naive Forecast" = "maroon4",
#"PCR Forecast" = "orangered",
#"PLS Forecast" = "lightskyblue",
#"SPCR Forecast" = "seashell")

filename_boxplot_sim2 <- paste(figures_dir,"/boxplot_sim2.pdf",sep="")

boxplot_sim2 <- ggplot(data_MSE, aes(x = Method, y = MSE, fill=Method)) +
  geom_boxplot() +
  scale_fill_manual(values = cbPalette) +
  labs(title = "Mean Squared Error of All Methods",
       x = NULL,
       y = "Mean Squared Error") +
  theme_grey() +
  theme(axis.text.x = element_blank(), legend.title = element_blank())

print(boxplot_sim2)

ggsave(filename = filename_boxplot_sim2, boxplot_sim2, width = 10, height = 6)

mean_coeff_pcr <- colMeans(coeff_pcr)
mean_coeff_pls <- colMeans(coeff_pls)
mean_coeff_spcr <- colMeans(coeff_spcr)

data_pcr <- data.frame(Index = seq_along(mean_coeff_pcr), Coeff_Values = mean_coeff_pcr)
data_pcr$dataset <- "Principal Component Regression"

# Create a column indicating the highlighting category
data_pcr$Highlight <- ifelse(data_pcr$Index <= p1, "p1",
                             ifelse(data_pcr$Index <= p1 + p2, "p2.1", "p2.2"))

data_pls <- data.frame(Index = seq_along(mean_coeff_pls), Coeff_Values = mean_coeff_pls)
data_pls$dataset <- "Partial Least Squares Regression"

# Create a column indicating the highlighting category
data_pls$Highlight <- ifelse(data_pls$Index <= p1, "p1",
                             ifelse(data_pls$Index <= p1 + p2, "p2.1", "p2.2"))

data_spcr <- data.frame(Index = seq_along(mean_coeff_spcr), Coeff_Values = mean_coeff_spcr)
data_spcr$dataset <- "Supervised Principal Component Regression"

# Create a column indicating the highlighting category
data_spcr$Highlight <- ifelse(data_spcr$Index <= p1, "p1",
                              ifelse(data_spcr$Index <= p1 + p2, "p2.1", "p2.2"))

# Combine the data into one dataframe
combined_data_coeff <- rbind(data_pcr, data_pls, data_spcr)
# Specify the order of levels in the "dataset" column
combined_data_coeff$dataset <- factor(combined_data_coeff$dataset,
                                      levels = c("Principal Component Regression",
                                                 "Partial Least Squares Regression",
                                                 "Supervised Principal Component Regression"))

combined_plot_coeff <- ggplot(combined_data_coeff,aes(x = Index, y = Coeff_Values, group = dataset,
                                                      shape = Highlight)) +
  geom_point(size=1) +
  scale_shape_manual(values = c(2,0,4)) +
  facet_wrap(~ dataset, scales = "free_y", ncol = 1) +
  labs(title = "Mean Coefficient Values for All Methods",
       x = "Index",
       y = "Coefficient Values") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",  # Set legend position to bottom
        legend.box = "horizontal")  # Remove legend margin)

print(combined_plot_coeff)

filename_combined_plot_coeff <- paste(figures_dir,"/combined_plot_coeff_sim2.pdf",sep="")

ggsave(filename = filename_combined_plot_coeff, combined_plot_coeff, width = 10, height = 6)

detach(data_coeff)
detach(param)
