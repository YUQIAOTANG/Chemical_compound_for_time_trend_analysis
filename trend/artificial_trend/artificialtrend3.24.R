library(readxl)
library(minpack.lm)
library(ggplot2)
library(stringr)
library(dplyr)
library(openxlsx)

setwd("your_file_path")
file_path <- "yourfile.xlsx"
data <- read_excel(file_path)

#select the group(cluster in 3 different samples)
selected_groups <- c("TNU_D_Dry", "TNU_L_Dry", "TNU_D_Wet", "TNU_L_Wet")

# set the color group
group_visual <- list(
  "TNU_D_Dry" = list(color = "#1f77b4", linetype = "solid", shape = 16),
  "TNU_L_Dry" = list(color = "#ff7f0e", linetype = "longdash", shape = 17),
  "TNU_D_Wet" = list(color = "#2ca02c", linetype = "dotted", shape = 15),
  "TNU_L_Wet" = list(color = "#d62728", linetype = "dotdash", shape = 18)
)

# model function(3 different models)based on the paper
monotonic_model <- function(time, C0, tau) {
  C0 * exp(-time / tau)
}

non_monotonic_model <- function(time, C0, Cx, tau1, tau2) {
  C0 + Cx * (-exp(-time / tau1) + exp(-time / tau2))
}

monotonic_increasing_model <- function(time, C0, Cf, tau) {
  (C0 - Cf) * exp(-time / tau) + Cf
}

# dataframe of the results including time constant and std
results <- data.frame(
  Compound = character(), Group = character(), Model = character(),
  C0 = numeric(), Cx = numeric(), Cf = numeric(), 
  Tau1 = numeric(), Tau2 = numeric(), 
  Sigma1 = numeric(), Sigma2 = numeric(),
  stringsAsFactors = FALSE
)

# data extract
for (compound_name in unique(data$name)) {
  compound_data <- subset(data, name == compound_name & Group %in% selected_groups)
  
  compound_data$ElapsedTime <- as.numeric(str_extract(compound_data$Condition, "[0-9.]+(?=d)"))
  compound_data$Concentration <- as.numeric(compound_data$Concentration)
  compound_data$std <- as.numeric(compound_data$std)
  
  compound_data <- compound_data %>%
    filter(!is.na(ElapsedTime), !is.na(Concentration), !is.na(std),
           is.finite(ElapsedTime), is.finite(Concentration),
           ElapsedTime >= 0, Concentration >= 0) %>%
    arrange(ElapsedTime)
  
  # innitialize the plot 
  p <- ggplot() +
    labs(x = "Elapsed time (days)", y = "Abundance (Peak area)",
         title = paste("Degradation Kinetics of", compound_name)) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.key.width = unit(1.5, "cm"),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_color_manual(
      name = "Condition",
      values = sapply(group_visual[selected_groups], function(x) x$color),
      breaks = selected_groups
    ) +
    scale_linetype_manual(
      name = "Condition",
      values = sapply(group_visual[selected_groups], function(x) x$linetype),
      breaks = selected_groups
    ) +
    scale_shape_manual(
      name = "Condition",
      values = sapply(group_visual[selected_groups], function(x) x$shape),
      breaks = selected_groups
    )
  
  # loop for all compounds
  for (group in selected_groups) {
    group_data <- subset(compound_data, Group == group)
    
    if (nrow(group_data) < 3 || !0 %in% group_data$ElapsedTime) {
      message(paste("Skipping", group, "- Insufficient data or no 0w measurement"))
      next
    }
    
    # set the starting and max data
    C0_fixed <- group_data$Concentration[group_data$ElapsedTime == 0][1]
    Cf_fixed <- group_data$Concentration[which.max(group_data$ElapsedTime)]
    
    concentration_diff <- diff(group_data$Concentration)
    is_monotonic_decreasing <- all(concentration_diff <= 0)
    is_monotonic_increasing <- all(concentration_diff >= 0)
    is_non_monotonic <- !(is_monotonic_decreasing || is_monotonic_increasing)
    
    # set the errror bar
    p <- p + geom_errorbar(
      data = group_data,
      aes(x = ElapsedTime, ymin = Concentration - std, ymax = Concentration + std,
          color = Group),  # reflect to the group varible
      width = 0.2
    )
    
    smooth_time <- seq(min(group_data$ElapsedTime), max(group_data$ElapsedTime), length.out = 100)
    
    # model determine
    if (is_monotonic_decreasing) {
      fit <- tryCatch(nlsLM(
        Concentration ~ C0_fixed * exp(-ElapsedTime / tau),
        data = group_data,
        start = list(tau = max(group_data$ElapsedTime) / 2),
        control = list(maxiter = 1000)),
        error = function(e) NULL)
      
      if (!is.null(fit)) {
        coeffs <- coef(fit)
        smooth_values <- monotonic_model(smooth_time, C0_fixed, coeffs["tau"])
        p <- p + geom_line(
          data = data.frame(ElapsedTime = smooth_time, Concentration = smooth_values, Group = group),
          aes(x = ElapsedTime, y = Concentration, color = Group, linetype = Group)  
        )
        
        # results store
        results <- rbind(results, data.frame(
          Compound = compound_name,
          Group = group,
          Model = "monotonic decreasing (C0 fixed)",
          C0 = C0_fixed, Cx = NA, Cf = NA,
          Tau1 = coeffs["tau"], Tau2 = NA,
          Sigma1 = summary(fit)$coefficients["tau", "Std. Error"], Sigma2 = NA))
      }
    } else if (is_monotonic_increasing) {
      fit <- tryCatch(nlsLM(
        Concentration ~ (C0_fixed - Cf_fixed) * exp(-ElapsedTime / tau) + Cf_fixed,
        data = group_data,
        start = list(tau = max(group_data$ElapsedTime) / 2),
        control = list(maxiter = 1000)),
        error = function(e) NULL)
      
      if (!is.null(fit)) {
        coeffs <- coef(fit)
        smooth_values <- monotonic_increasing_model(smooth_time, C0_fixed, Cf_fixed, coeffs["tau"])
        p <- p + geom_line(
          data = data.frame(ElapsedTime = smooth_time, Concentration = smooth_values, Group = group),
          aes(x = ElapsedTime, y = Concentration, color = Group, linetype = Group)
        )
        
        results <- rbind(results, data.frame(
          Compound = compound_name, Group = group, Model = "monotonic increasing (C0,Cf fixed)",
          C0 = C0_fixed, Cx = NA, Cf = Cf_fixed,
          Tau1 = coeffs["tau"], Tau2 = NA,
          Sigma1 = summary(fit)$coefficients["tau", "Std. Error"], Sigma2 = NA))
      }
    } else if (is_non_monotonic) {
      large_tau_threshold <- 2 * (max(group_data$ElapsedTime) - min(group_data$ElapsedTime))
      
      fit <- tryCatch(nlsLM(
        Concentration ~ C0_fixed + Cx * (-exp(-ElapsedTime / tau1) + exp(-ElapsedTime / tau2)),
        data = group_data,
        start = list(Cx = Cf_fixed - C0_fixed,
                     tau1 = max(group_data$ElapsedTime)/4,
                     tau2 = max(group_data$ElapsedTime)/2),
        control = list(maxiter = 2000)),
        error = function(e) NULL)
      
      correlation <- cor(group_data$ElapsedTime, group_data$Concentration)
      
      if (!is.null(fit) && (abs(correlation) < 0.7)) {
        coeffs <- coef(fit)
        smooth_values <- non_monotonic_model(smooth_time, C0_fixed, coeffs["Cx"], coeffs["tau1"], coeffs["tau2"])
        p <- p + geom_line(
          data = data.frame(ElapsedTime = smooth_time, Concentration = smooth_values, Group = group),
          aes(x = ElapsedTime, y = Concentration, color = Group, linetype = Group)
        )
        
        results <- rbind(results, data.frame(
          Compound = compound_name, Group = group, Model = "non-monotonic (C0 fixed)",
          C0 = C0_fixed, Cx = coeffs["Cx"], Cf = NA,
          Tau1 = coeffs["tau1"], Tau2 = coeffs["tau2"],
          Sigma1 = summary(fit)$coefficients["tau1", "Std. Error"],
          Sigma2 = summary(fit)$coefficients["tau2", "Std. Error"]))
      } else {
        # fall back due to unreasonable results
        if (correlation >= 0.7) {
          mono_fit <- tryCatch(nlsLM(
            Concentration ~ (C0_fixed - Cf_fixed) * exp(-ElapsedTime / tau) + Cf_fixed,
            data = group_data,
            start = list(tau = max(group_data$ElapsedTime) / 2),
            control = list(maxiter = 1000)),
            error = function(e) NULL)
          
          if (!is.null(mono_fit)) {
            coeffs_mono <- coef(mono_fit)
            smooth_values <- monotonic_increasing_model(smooth_time, C0_fixed, Cf_fixed, coeffs_mono["tau"])
            
            p <- p + geom_line(
              data = data.frame(ElapsedTime = smooth_time, Concentration = smooth_values, Group = group),
              aes(x = ElapsedTime, y = Concentration, color = Group),  
              linetype = "dashed"  # set the special line type for fall back data
            )
            
            results <- rbind(results, data.frame(
              Compound = compound_name, Group = group,
              Model = "monotonic increasing (fallback)",
              C0 = C0_fixed, Cx = NA, Cf = Cf_fixed,
              Tau1 = coeffs_mono["tau"], Tau2 = NA,
              Sigma1 = summary(mono_fit)$coefficients["tau", "Std. Error"], Sigma2 = NA))
          }
        } else if (correlation <= -0.7) {
          mono_fit <- tryCatch(nlsLM(
            Concentration ~ C0_fixed * exp(-ElapsedTime / tau),
            data = group_data,
            start = list(tau = max(group_data$ElapsedTime) / 2),
            control = list(maxiter = 1000)),
            error = function(e) NULL)
          
          if (!is.null(mono_fit)) {
            coeffs_mono <- coef(mono_fit)
            smooth_values <- monotonic_model(smooth_time, C0_fixed, coeffs_mono["tau"])
            
            p <- p + geom_line(
              data = data.frame(ElapsedTime = smooth_time, Concentration = smooth_values, Group = group),
              aes(x = ElapsedTime, y = Concentration, color = Group), 
              linetype = "dashed" 
            )
            
            results <- rbind(results, data.frame(
              Compound = compound_name, Group = group,
              Model = "monotonic decreasing (fallback)",
              C0 = C0_fixed, Cx = NA, Cf = NA,
              Tau1 = coeffs_mono["tau"], Tau2 = NA,
              Sigma1 = summary(mono_fit)$coefficients["tau", "Std. Error"], Sigma2 = NA))
          }
        }
      }
    }
    
    # 数据点（修复映射）
    p <- p + geom_point(
      data = group_data,
      aes(x = ElapsedTime, y = Concentration, color = Group, shape = Group),
      size = 3
    )
  }
  
  # merge the plots
  p <- p + guides(
    color = guide_legend(
      title = "Condition",
      override.aes = list(
        linetype = sapply(group_visual[selected_groups], function(x) x$linetype),
        shape = sapply(group_visual[selected_groups], function(x) x$shape)
      )
    ),
    linetype = "none",  
    shape = "none"     
  )
    
    # save the results 
    file_name <- paste0(gsub("[^A-Za-z0-9]", "_", compound_name), ".png")
    ggsave(file_name, p, width = 10, height = 6, dpi = 300)
}

write.xlsx(results, "kinetic_results.xlsx", rowNames = FALSE)
print("Results saved to kinetic_results.xlsx")