# load R packages
library(readxl)
library(minpack.lm)
library(ggplot2)
library(stringr)
library(dplyr)
library(openxlsx)

# seting working directory
setwd("/Users/tangyuqiao/Desktop/EAI/downstream/Data Analysis w Tugce/time trend/naturaltrend/")
file_path <- "groupdata.xlsx"
data <- read_excel(file_path)

# select groups
selected_groups <- c("TNU Dry", "TNU Wet")

# define models
monotonic_model <- function(time, C0, tau) {
  C0 * exp(-time / tau)
}

non_monotonic_model <- function(time, C0, Cx, tau1, tau2) {
  C0 + Cx * (-exp(-time / tau1) + exp(-time / tau2))
}

monotonic_increasing_model <- function(time, C0, Cf, tau) {
  (C0 - Cf) * exp(-time / tau) + Cf
}

# initialize results storage
results <- data.frame(
  Compound = character(),
  Group = character(),
  Model = character(),
  C0 = numeric(),
  Cx = numeric(),
  Cf = numeric(),
  Tau1 = numeric(),
  Tau2 = numeric(),
  Sigma1 = numeric(),
  Sigma2 = numeric(),
  stringsAsFactors = FALSE
)

# loop through all compounds
for (compound_name in unique(data$name)) {
  compound_data <- subset(data, name == compound_name & Group %in% selected_groups)
  
  # data processing(unnecessary if you already done)
  compound_data$ElapsedTime <- as.numeric(str_extract(compound_data$Condition, "\\d+(?=w)"))
  compound_data$Concentration <- as.numeric(compound_data$Concentration)
  compound_data$Std <- as.numeric(compound_data$std)
  
  compound_data <- compound_data %>%
    filter(!is.na(ElapsedTime), !is.na(Concentration), !is.na(Std),
           is.finite(ElapsedTime), is.finite(Concentration),
           ElapsedTime >= 0, Concentration >= 0) %>%
    arrange(ElapsedTime)
  
  # creat plot project
  p <- ggplot() +
    labs(x = "Elapsed time (weeks)", y = "Abundance (Peak area)",
         title = paste("Degradation Kinetics of", compound_name))
  
  for (group in selected_groups) {
    group_data <- subset(compound_data, Group == group)
    
    if (nrow(group_data) < 3 || !0 %in% group_data$ElapsedTime) {
      message(paste("Skipping", group, "- Insufficient data or no 0w measurement"))
      next
    }
    
    # set C0=t0
    C0_initial <- group_data$Concentration[group_data$ElapsedTime == 0][1]
    
    # check the trend of concentration changes
    concentration_diff <- diff(group_data$Concentration)
    is_monotonic_decreasing <- all(concentration_diff <= 0)
    is_monotonic_increasing <- all(concentration_diff >= 0)
    is_non_monotonic <- !(is_monotonic_decreasing || is_monotonic_increasing)
    
    # add error bar
    p <- p + geom_errorbar(data = group_data,
                           aes(x = ElapsedTime, ymin = Concentration - Std, ymax = Concentration + Std,
                               color = factor(Group)), width = 0.2)
    
    smooth_time <- seq(min(group_data$ElapsedTime), max(group_data$ElapsedTime), length.out = 100)
    fitted <- FALSE
    
    # fit 3 different models
    if (is_monotonic_decreasing) {
      tryCatch({
        fit <- nlsLM(Concentration ~ monotonic_model(ElapsedTime, C0, tau),
                     data = group_data,
                     start = list(C0 = C0_initial, tau = max(group_data$ElapsedTime) / 2),
                     control = list(maxiter = 1000))
        
        coeffs <- coef(fit)
        fit_summary <- summary(fit)
        sigma1 <- fit_summary$coefficients["tau", "Std. Error"]
        
        smooth_values <- monotonic_model(smooth_time, coeffs["C0"], coeffs["tau"])
        p <- p + geom_line(data = data.frame(ElapsedTime = smooth_time, Concentration = smooth_values, Group = group),
                           aes(x = ElapsedTime, y = Concentration, color = Group))
        
        results <- rbind(results, data.frame(
          Compound = compound_name,
          Group = group,
          Model = "monotonic decreasing",
          C0 = coeffs["C0"],
          Cx = NA,
          Cf = NA,
          Tau1 = coeffs["tau"],
          Tau2 = NA,
          Sigma1 = sigma1,
          Sigma2 = NA
        ))
      }, error = function(e) { 
        message(paste("Monotonic decreasing model failed for", group, ":", e$message))
      })
    }  
    
    
    if (is_monotonic_increasing) {
      tryCatch({
        fit <- nlsLM(Concentration ~ monotonic_increasing_model(ElapsedTime, C0, Cf, tau),
                     data = group_data,
                     start = list(C0 = C0_initial, Cf = max(group_data$Concentration, na.rm = TRUE),
                                  tau = max(group_data$ElapsedTime) / 2),
                     control = list(maxiter = 1000))
        
        coeffs <- coef(fit)
        fit_summary <- summary(fit)
        sigma1 <- fit_summary$coefficients["tau", "Std. Error"]
        
        smooth_values <- monotonic_increasing_model(smooth_time, coeffs["C0"], coeffs["Cf"], coeffs["tau"])
        p <- p + geom_line(data = data.frame(ElapsedTime = smooth_time, Concentration = smooth_values, Group = group),
                           aes(x = ElapsedTime, y = Concentration, color = Group))
        
        results <- rbind(results, data.frame(
          Compound = compound_name,
          Group = group,
          Model = "monotonic increasing",
          C0 = coeffs["C0"],
          Cx = NA,
          Cf = coeffs["Cf"],
          Tau1 = coeffs["tau"],
          Tau2 = NA,
          Sigma1 = sigma1,
          Sigma2 = NA
        ))
      }, error = function(e) { 
        message(paste("Monotonic increasing model failed for", group, ":", e$message))
      })
    }  
    
    
    if (!fitted && is_non_monotonic) {
      tryCatch({
        fit <- nlsLM(Concentration ~ non_monotonic_model(ElapsedTime, C0, Cx, tau1, tau2),
                     data = group_data,
                     start = list(C0 = C0_initial, Cx = max(group_data$Concentration) - C0_initial,
                                  tau1 = max(group_data$ElapsedTime) / 4, tau2 = max(group_data$ElapsedTime) / 2),
                     control = list(maxiter = 2000))
        
        coeffs <- coef(fit)
        fit_summary <- summary(fit)
        sigma1 <- fit_summary$coefficients["tau1", "Std. Error"]
        sigma2 <- fit_summary$coefficients["tau2", "Std. Error"]
        
        smooth_values <- non_monotonic_model(smooth_time, coeffs["C0"], coeffs["Cx"], coeffs["tau1"], coeffs["tau2"])
        p <- p + geom_line(data = data.frame(ElapsedTime = smooth_time, Concentration = smooth_values, Group = group),
                           aes(x = ElapsedTime, y = Concentration, color = Group))
        
        results <- rbind(results, data.frame(
          Compound = compound_name,
          Group = group,
          Model = "non-monotonic",
          C0 = coeffs["C0"],
          Cx = coeffs["Cx"],
          Cf = NA,
          Tau1 = coeffs["tau1"],
          Tau2 = coeffs["tau2"],
          Sigma1 = sigma1,
          Sigma2 = sigma2
        ))
      }, error = function(e) { 
        message(paste("Non-monotonic model failed for", group, ":", e$message))
      })
    }  
    
    p <- p + geom_point(data = group_data, aes(x = ElapsedTime, y = Concentration, color = factor(Group), shape = factor(Group)), size = 3)
  }
  
  p <- p + theme(legend.position = "right") +
    guides(color = guide_legend(title = "Group"), shape = guide_legend(title = "Group"))
  
  file_name <- paste0(gsub("[^A-Za-z0-9]", "_", compound_name), ".png")
  ggsave(file_name, p, width = 10, height = 6, dpi = 300)
}

# output the results to Excel
write.xlsx(results, "kinetic_results.xlsx", rowNames = FALSE)

print("Results saved to kinetic_results.xlsx")
