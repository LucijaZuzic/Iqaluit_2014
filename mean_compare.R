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