# Load required libraries
library(readxl)
library(minpack.lm)
library(ggplot2)
library(stringr)
library(dplyr)

# Set working directory and read data
setwd("/Users/tangyuqiao/Desktop/EAI/downstream/Data Analysis w Tugce/time trend/")
file_path <- "grouped_data_TEST.xlsx"
data <- read_excel(file_path)

# Select target compound and groups
compound_name <- "6PPD-quinone (6PPD-Q)"
selected_groups <- c("TN_Dark_Dry", "TN_Light_Dry", "TN_Dark_Wet", "TN_Light_Wet")
compound_data <- subset(data, name == compound_name & Group %in% selected_groups)

# Data preprocessing
compound_data$ElapsedTime <- as.numeric(str_extract(compound_data$Condition, "\\d+\\.?\\d*(?=d)"))
compound_data$Concentration <- as.numeric(compound_data$Concentration)
compound_data <- compound_data %>%
  filter(!is.na(ElapsedTime), !is.na(Concentration),
         is.finite(ElapsedTime), is.finite(Concentration),
         ElapsedTime >= 0, Concentration >= 0) %>%
  arrange(ElapsedTime)

# Define models
monotonic_model <- function(time, C0, tau) {
  C0 * exp(-time / tau)
}

non_monotonic_model <- function(time, C0, Cx, tau1, tau2) {
  C0 + Cx * (-exp(-time / tau1) + exp(-time / tau2))
}

# Initialize results storage
results <- data.frame(
  Group = character(),
  Model = character(),
  C0 = numeric(),
  Cx = numeric(),
  Tau1 = numeric(),
  Tau2 = numeric(),
  Sigma1 = numeric(),
  Sigma2 = numeric(),
  stringsAsFactors = FALSE
)

# Create plot
p <- ggplot() +
  labs(x = "Elapsed time (days)", y = "Concentration (μg/g)",
       title = paste("Degradation Kinetics of", compound_name))

# Iterate through groups
for (i in seq_along(selected_groups)) {
  group <- selected_groups[i]
  group_data <- subset(compound_data, Group == group)
  
  if (nrow(group_data) < 3) {
    message(paste("Skipping", group, "- Insufficient data points"))
    next
  }
  
  # check the trend type
  is_monotonic <- all(diff(group_data$Concentration) <= 0)
  fit_success <- FALSE
  
  if (!is_monotonic && nrow(group_data) >= 3) {
    # **LOESS smoothy only for non-monotonic data**
    span_value <- ifelse(nrow(group_data) < 10, 1, 0.75)  # fit the mini datasets
    loess_fit <- loess(Concentration ~ ElapsedTime, data = group_data, span = span_value)
    time_interpolated <- seq(min(group_data$ElapsedTime), max(group_data$ElapsedTime), length.out = 100)
    smooth_concentration <- predict(loess_fit, newdata = data.frame(ElapsedTime = time_interpolated))
    
    # initial data for non-monotonic data
    group_data_interpolated <- data.frame(
      ElapsedTime = time_interpolated,
      Concentration = smooth_concentration,
      Group = group
    )
    
    # **fit the non-monotonic model**
    C0_start <- max(group_data_interpolated$Concentration, na.rm = TRUE)
    Cx_start <- C0_start - min(group_data_interpolated$Concentration, na.rm = TRUE)
    
    tryCatch({
      fit <- nlsLM(
        Concentration ~ non_monotonic_model(ElapsedTime, C0, Cx, tau1, tau2),
        data = group_data_interpolated,
        start = list(C0 = C0_start, Cx = Cx_start, tau1 = 2, tau2 = 15),
        lower = c(0, 0, 0.1, 1),  
        upper = c(Inf, Inf, 30, 50),
        control = list(maxiter = 10000, tol = 1e-6)
      )
      
      # extract parameters
      coeffs <- coef(summary(fit))
      results <- rbind(results, data.frame(
        Group = group,
        Model = "Non-monotonic",
        C0 = coeffs["C0", "Estimate"],
        Cx = coeffs["Cx", "Estimate"],
        Tau1 = coeffs["tau1", "Estimate"],
        Tau2 = coeffs["tau2", "Estimate"],
        Sigma1 = coeffs["tau1", "Std. Error"],
        Sigma2 = coeffs["tau2", "Std. Error"]
      ))
      fit_success <- TRUE
      
      # **generazi pre curve**
      time_range <- seq(0, max(group_data_interpolated$ElapsedTime), length.out = 100)
      pred_df <- data.frame(
        ElapsedTime = time_range,
        Concentration = predict(fit, newdata = data.frame(ElapsedTime = time_range)),
        Group = group
      )
      
      # **add curve to  ggplot**
      p <- p + geom_line(data = pred_df, aes(x = ElapsedTime, y = Concentration, color = factor(Group)))
      
    }, error = function(e) {
      message(paste("Non-monotonic model failed for", group, ":", e$message))
    })
  }
  
  # **fit the monotnic model**
  if (is_monotonic) {
    C0_start <- max(group_data$Concentration, na.rm = TRUE)
    tryCatch({
      fit <- nlsLM(
        Concentration ~ monotonic_model(ElapsedTime, C0, tau),
        data = group_data,
        start = list(C0 = C0_start, tau = 5),
        control = list(maxiter = 1000)
      )
      
      # extract parameters
      coeffs <- coef(summary(fit))
      results <- rbind(results, data.frame(
        Group = group,
        Model = "Monotonic",
        C0 = coeffs["C0", "Estimate"],
        Cx = NA,
        Tau1 = coeffs["tau", "Estimate"],
        Tau2 = NA,
        Sigma1 = coeffs["tau", "Std. Error"],
        Sigma2 = NA
      ))
      
      # **generize pre curve**
      time_range <- seq(0, max(group_data$ElapsedTime), length.out = 100)
      pred_df <- data.frame(
        ElapsedTime = time_range,
        Concentration = predict(fit, newdata = data.frame(ElapsedTime = time_range)),
        Group = group
      )
      
      # **add predict curve ggplot**
      p <- p + geom_line(data = pred_df, aes(x = ElapsedTime, y = Concentration, color = factor(Group)))
      
    }, error = function(e) {
      message(paste("Monotonic model failed for", group, ":", e$message))
    })
  }
  
  # **using different shape for the origin concentration**
  p <- p + geom_point(data = group_data, aes(x = ElapsedTime, y = Concentration, color = factor(Group), shape = factor(Group)), size = 3)
}

p <- p +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold")
  ) +
  guides(color = guide_legend(title = "Group"), shape = guide_legend(title = "Group"))

# out put the plot and the result
print(results)
file_name <- paste0(gsub("[^A-Za-z0-9]", "_", compound_name), ".png")
ggsave(file_name, p, width = 10, height = 6, dpi = 300)
