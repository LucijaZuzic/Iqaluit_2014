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