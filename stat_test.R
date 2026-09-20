# Clear environment
rm(list=ls())

library(tidyverse)
library(circular)

## The R script for Žužić et al. manuscript

# Set path to the working directory containing prepared data set

getCurrentFileLocation <-  function()
{
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col = value, into = c("key", "value"), sep = "=", fill = 'right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file) == 0) {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}
setwd(getCurrentFileLocation())

# --- Helper Function to Extract Circular Peak Errors ---
# This function isolates the top 25% of errors (Q3 and above) to represent the "peak error times"
get_circular_peaks <- function(storm_num, dfpv) {
  file_name <- paste0("afullframe_", storm_num, "new.csv")
  
  if (!file.exists(file_name)) {
    stop(paste("File not found:", file_name))
  }
  
  df <- read_csv(file_name, show_col_types = FALSE)
  
  # Calculate decimal hour of the day
  # total seconds accumulates over days; modulo 86400 resets it daily
  df <- df %>%
    mutate(
      sec_of_day = `total seconds` %% 86400,
      hour_decimal = sec_of_day / 3600
    )
  
  # Filter for the highest errors (disturbed periods)
  qp <- quantile(df$`horizontal(deg)`, dfpv, na.rm = TRUE)
  df_peaks <- df %>% filter(`horizontal(deg)` >= qp)
  
  # Convert to a circular object (24-hour clock)
  circ_data <- circular(df_peaks$hour_decimal, 
                        type = "angles", 
                        units = "hours", 
                        template = "clock24", 
                        modulo = "asis")
  return(circ_data)
}

for (dfpv1 in c(0, 0.25, 0.5, 0.75)) {
  # --- 1. Load the Circular Data for All Locations ---
  cat("Loading data and converting to circular objects...\n")
  storm1 <- get_circular_peaks("", dfpv1)
  storm2 <- get_circular_peaks("Darwin_", dfpv1)
  print(dfpv1)
  
  # --- 2. Rayleigh Test for Uniformity ---
  # Tests if the errors are randomly distributed across the 24 hours (Null Hypothesis)
  # or if they have a significant directional (temporal) cluster (Alternative Hypothesis).
  cat("\n======================================================\n")
  cat("RAYLEIGH TEST FOR UNIFORMITY (p < 0.05 indicates clustering)\n")
  cat("======================================================\n")
  
  print(rayleigh.test(storm1))
  print(rayleigh.test(storm2))
  
  # --- 3. Watson-Williams Tests (Circular ANOVA) ---
  # Compares the mean angular directions (peak error times) between different locations.
  cat("\n======================================================\n")
  cat("WATSON-WILLIAMS TESTS (Comparing Subclasses)\n")
  cat("======================================================\n")
  
  
  cat("\n--- Test: Iqaluit vs. Darwin ---\n")
  cat("Hypothesis: They belong to the SAME subclass (Expect p > 0.05)\n")
  ww_1_vs_2 <- watson.williams.test(list(storm1, storm2))
  print(ww_1_vs_2)
  
  # --- 4. Mardia-Watson-Wheeler Tests (Non-Parametric) ---
  # The MWW test evaluates whether two or more circular distributions are identical.
  # It is robust against low concentration parameters and non-von Mises distributions.
  
  cat("\n======================================================\n")
  cat("MARDIA-WATSON-WHEELER TESTS (Non-Parametric Comparisons)\n")
  cat("======================================================\n")
  
  cat("\n--- Test: Iqaluit vs. Darwin ---\n")
  cat("Hypothesis: Distributions are identical (Expect p > 0.05 if same subclass)\n")
  mww_1_vs_2 <- watson.wheeler.test(list(storm1, storm2))
  print(mww_1_vs_2)
  break
}