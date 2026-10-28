# Clear environment
rm(list=ls())

## The R script for Žužić et al. manuscript

# Set path to the working directory containing prepared data set

library(tidyverse)
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

nks <- c("hr", "len", "avg")
rdim <- list(hr = 0, len = 0, avg = 2)

storm_df <- read_csv("afullframe_new.csv", show_col_types = FALSE)
storm_df2 <- read_csv("afullframe_Darwin_new.csv", show_col_types = FALSE)
stats_storm <- storm_df %>%
  group_by(hr) %>%
  summarise(
    len = n(),
    avg = mean(`horizontal(deg)`, na.rm = TRUE),
  ) %>%
  complete(hr = 0:23, fill = list(len=0, avg=0))
stats_storm2 <- storm_df2 %>%
  group_by(hr) %>%
  summarise(
    len = n(),
    avg = mean(`horizontal(deg)`, na.rm = TRUE),
  ) %>%
  complete(hr = 0:23, fill = list(len=0, avg=0))

# LaTeX table print loop (year)
for (vix in 3:length(nks)) {
  col_name <- nks[vix]
  cat(sprintf("$%s$ $%s$ $%s$ $%s$ $%s$\n", col_name, which.max(stats_storm[[col_name]]) - 1, 
              round(max(stats_storm[[col_name]]), rdim[[col_name]]), 
              which.min(stats_storm[[col_name]]) - 1, round(min(stats_storm[[col_name]]), rdim[[col_name]])))
}
for (vix in 3:length(nks)) {
  col_name <- nks[vix]
  cat(sprintf("$%s$ $%s$ $%s$ $%s$ $%s$\n", col_name, which.max(stats_storm2[[col_name]]) - 1, 
              round(max(stats_storm2[[col_name]]), rdim[[col_name]]), 
              which.min(stats_storm2[[col_name]]) - 1, round(min(stats_storm2[[col_name]]), rdim[[col_name]])))
}
ni1 <- c()
ni2 <- c()
for (h in 0:23) {
  h1 <- sprintf("%02d", h)
  pvals <- paste0("$", h1, "$")
  for (vix in 3:length(nks)) {
    col_name <- nks[vix]
    val <- stats_storm[[col_name]][h + 1]
    ni1 <- c(ni1, val)
    strv <- as.character(round(val, rdim[[col_name]]))
    s1 <- as.character(round(max(stats_storm[[col_name]]), rdim[[col_name]]))
    s2 <- as.character(round(min(stats_storm[[col_name]]), rdim[[col_name]]))
    
    if (strv == s1) pvals <- c(pvals, paste0("$\\mathbf{", strv, "}$"))
    else if (strv == s2) pvals <- c(pvals, paste0("$\\underline{\\mathbf{", strv, "}}$"))
    else pvals <- c(pvals, paste0("$", strv, "$"))
  }
  for (vix in 3:length(nks)) {
    col_name <- nks[vix]
    val <- stats_storm2[[col_name]][h + 1]
    ni2 <- c(ni2, val)
    strv <- as.character(round(val, rdim[[col_name]]))
    s1 <- as.character(round(max(stats_storm2[[col_name]]), rdim[[col_name]]))
    s2 <- as.character(round(min(stats_storm2[[col_name]]), rdim[[col_name]]))
    
    if (strv == s1) pvals <- c(pvals, paste0("$\\mathbf{", strv, "}$"))
    else if (strv == s2) pvals <- c(pvals, paste0("$\\underline{\\mathbf{", strv, "}}$"))
    else pvals <- c(pvals, paste0("$", strv, "$"))
  }
  cat(paste(pvals, collapse = " & "), " \\\\ \\hline\n")
}
print(ni1)
print(ni2)
print(wilcox.test(ni1, ni2))
print(t.test(ni1, ni2))


# 1. Evaluate Median (Q2) across hours to show outlier robustness
stats_med1 <- storm_df %>%
  group_by(hr) %>%
  summarise(med = median(`horizontal(deg)`, na.rm = TRUE)) %>%
  complete(hr = 0:23, fill = list(med = 0))

stats_med2 <- storm_df2 %>%
  group_by(hr) %>%
  summarise(med = median(`horizontal(deg)`, na.rm = TRUE)) %>%
  complete(hr = 0:23, fill = list(med = 0))

cat("\n--- Hourly Median Comparison (Unfiltered) ---\n")
cat("Iqaluit Median Range:", min(stats_med1$med), "to", max(stats_med1$med), 
    "Overall Median:", median(storm_df$`horizontal(deg)`, na.rm = TRUE), "\n")
cat("Darwin Median Range:", min(stats_med2$med), "to", max(stats_med2$med), 
    "Overall Median:", median(storm_df2$`horizontal(deg)`, na.rm = TRUE), "\n")
print(wilcox.test(stats_med1$med, stats_med2$med))
print(t.test(stats_med1$med, stats_med2$med))

# 2. Evaluate Filtered (< 10 m) Hourly Means to remove >1000 km outliers
filt_threshold <- 10
stats_filt1 <- storm_df %>%
  filter(`horizontal(deg)` < filt_threshold) %>%
  group_by(hr) %>%
  summarise(avg_filt = mean(`horizontal(deg)`, na.rm = TRUE)) %>%
  complete(hr = 0:23, fill = list(avg_filt = 0))

stats_filt2 <- storm_df2 %>%
  filter(`horizontal(deg)` < filt_threshold) %>%
  group_by(hr) %>%
  summarise(avg_filt = mean(`horizontal(deg)`, na.rm = TRUE)) %>%
  complete(hr = 0:23, fill = list(avg_filt = 0))

cat("\n--- Filtered (< 10 m) Hourly Mean Comparison ---\n")
cat("Iqaluit Filtered Mean:", mean(stats_filt1$avg_filt), 
    "Min:", min(stats_filt1$avg_filt), "(Hr", which.min(stats_filt1$avg_filt)-1, ")",
    "Max:", max(stats_filt1$avg_filt), "(Hr", which.max(stats_filt1$avg_filt)-1, ")\n")
cat("Darwin Filtered Mean:", mean(stats_filt2$avg_filt), 
    "Min:", min(stats_filt2$avg_filt), "(Hr", which.min(stats_filt2$avg_filt)-1, ")",
    "Max:", max(stats_filt2$avg_filt), "(Hr", which.max(stats_filt2$avg_filt)-1, ")\n")
print(wilcox.test(stats_filt1$avg_filt, stats_filt2$avg_filt))
print(t.test(stats_filt1$avg_filt, stats_filt2$avg_filt))

# 3. Quantify proportion of extreme outliers in Iqaluit
n_total_iqal <- sum(!is.na(storm_df$`horizontal(deg)`))
n_gt10_iqal  <- sum(storm_df$`horizontal(deg)` >= 10, na.rm = TRUE)
n_gt1000_iqal <- sum(storm_df$`horizontal(deg)` >= 1000, na.rm = TRUE)
cat(sprintf("\nIqaluit Outliers >= 10m: %d (%.4f%%), >= 1000m: %d (%.4f%%)\n %d",
            n_gt10_iqal, 100 * n_gt10_iqal / n_total_iqal,
            n_gt1000_iqal, 100 * n_gt1000_iqal / n_total_iqal, n_total_iqal))

# 4. Inspect timestamps (DOY, Date, UTC hour) of extreme outliers and cross-check with GFZ Kp index
outliers_iqal <- storm_df %>%
  filter(`horizontal(deg)` >= 1000) %>%
  mutate(
    DOY = doy + 1,
    kp_block = (hr %/% 3) * 3
  )

outlier_summary <- outliers_iqal %>%
  group_by(y, m, d, DOY, hr, kp_block) %>%
  summarise(
    n_epochs = n(),
    max_err_m = max(`horizontal(deg)`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(DOY, hr)

cat("\n--- Iqaluit Extreme Outlier (>= 1000 m) Dates & UTC Hours ---\n")
print(outlier_summary, n = 100)

# Cross-reference with GFZ Potsdam Kp index (https://kp.gfz.de/en/data)
kp_file <- "kp_2014.txt"
if (!file.exists(kp_file)) {
  tryCatch({
    download.file("https://kp.gfz.de/app/files/Kp_ap_Ap_SN_F107_since_1932.txt",
                  kp_file, quiet = TRUE)
  }, error = function(e) {
    message("Could not auto-download GFZ Kp file; place Kp_ap_Ap_SN_F107_since_1932.txt as kp_2014.txt")
  })
}

if (file.exists(kp_file)) {
  kp_lines <- readLines(kp_file)
  kp_lines <- kp_lines[!grepl("^#", kp_lines)]
  kp_raw <- read.table(text = kp_lines, header = FALSE)
  kp_2014 <- kp_raw %>%
    filter(V1 == 2014) %>%
    select(y = V1, m = V2, d = V3,
           Kp0 = V8, Kp3 = V9, Kp6 = V10, Kp9 = V11,
           Kp12 = V12, Kp15 = V13, Kp18 = V14, Kp21 = V15, Ap = V24) %>%
    pivot_longer(cols = starts_with("Kp"), names_to = "block_str", values_to = "Kp") %>%
    mutate(kp_block = as.numeric(sub("Kp", "", block_str))) %>%
    group_by(y, m, d) %>%
    mutate(Kp_daily_max = max(Kp, na.rm = TRUE)) %>%
    ungroup()
  
  outliers_kp <- outliers_iqal %>%
    left_join(kp_2014, by = c("y", "m", "d", "kp_block"))
  
  n_out <- nrow(outliers_kp)
  n_storm <- sum(outliers_kp$Kp >= 5, na.rm = TRUE)
  n_active <- sum(outliers_kp$Kp_daily_max >= 4, na.rm = TRUE)
  n_quiet <- sum(outliers_kp$Kp_daily_max < 4, na.rm = TRUE)
  
  cat("\n--- GFZ Kp Verification of Iqaluit Extreme Outliers (>= 1000 m) ---\n")
  cat(sprintf("Total >= 1000 m epochs: %d across %d date-hour intervals\n", n_out, nrow(outlier_summary)))
  cat(sprintf("Storm conditions (Kp >= 5): %d epochs (%.2f%%)\n", n_storm, 100 * n_storm / n_out))
  cat(sprintf("Active/Substorm conditions (Kp_daily_max >= 4): %d epochs (%.2f%%)\n", n_active, 100 * n_active / n_out))
  cat(sprintf("Quiet conditions (Kp_daily_max < 4, processing artifacts): %d epochs (%.2f%%)\n", n_quiet, 100 * n_quiet / n_out))
  
  cat(sprintf("Active/Storm days (Kp_daily_max >= 4): %d epochs (%.2f%%), of which %d epochs (%.2f%%) have 3h Kp >= 5\n", 
              n_active, 100 * n_active / n_out, n_storm, 100 * n_storm / n_out))
  
  # Show exact storm intervals (Kp >= 5)
  cat("\n--- Extreme Outlier Intervals with Kp >= 5 ---\n")
  print(outliers_kp %>%
          filter(Kp >= 5) %>%
          group_by(y, m, d, DOY, hr, Kp, Kp_daily_max) %>%
          summarise(n_epochs = n(), max_err_m = max(`horizontal(deg)`), .groups = "drop"))
  
  # Show exact storm intervals (3h Kp >= 5 or daily Kp_max >= 5)
  cat("\n--- Extreme Outlier Intervals with 3h Kp >= 5 or Kp_daily_max >= 5 ---\n")
  print(outliers_kp %>%
          filter(Kp >= 5 | Kp_daily_max >= 5) %>%
          group_by(y, m, d, DOY, hr, Kp, Kp_daily_max) %>%
          summarise(n_epochs = n(), max_err_m = max(`horizontal(deg)`), .groups = "drop"))
}