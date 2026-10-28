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

set.seed(42) # Ensure reproducible random tie-breaking in watson.wheeler.test

for (dfpv1 in c(0, 0.75)) {
  # 0 = full dataset, 0.75 = top 25% peak errors (Q3+)
  cat("\n######################################################\n")
  cat("QUANTILE THRESHOLD:", dfpv1, "\n")
  cat("######################################################\n")
  
  storm1 <- get_circular_peaks("", dfpv1)
  storm2 <- get_circular_peaks("Darwin_", dfpv1)
  
  # Circular Summary Parameters (Mean Direction & Resultant Length R-bar)
  cat("Iqaluit Circular Mean (hrs):", as.numeric(mean.circular(storm1)) %% 24, 
      "Resultant Length (rho):", rho.circular(storm1), "\n")
  cat("Darwin Circular Mean (hrs):", as.numeric(mean.circular(storm2)) %% 24, 
      "Resultant Length (rho):", rho.circular(storm2), "\n")
  
  # Rayleigh Test for Uniformity
  print(rayleigh.test(storm1))
  print(rayleigh.test(storm2))
  
  # Watson-Williams Test
  print(watson.williams.test(list(storm1, storm2)))
  
  # Mardia-Watson-Wheeler Test
  print(watson.wheeler.test(list(storm1, storm2)))
}