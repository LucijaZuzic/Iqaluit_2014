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

library(dplyr)
library(tidyr)
library(lubridate)
library(geosphere)
library(ggplot2)
library(readr)
library(stringr)
library(extrafont) # Added extrafont

# (Optional: run loadfonts() if your PDF output isn't rendering the font correctly)
loadfonts(device = "pdf", quiet = TRUE)

# 1. Initialization and reference coordinates

ref_coords <- list(
  lon = 131.132744,
  lat = -12.843697,
  height = 125.1
)

# Shared theme
bigger_theme <- theme_classic(base_size = 20, base_family = "DejaVu Sans") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 20, color = "black"),       # Pure black text
    axis.ticks = element_line(color = "black", linewidth = 1),  # Pure black tick marks
    axis.line = element_line(color = "black", linewidth = 1),   # Pure black axis lines
    legend.text = element_text(size = 20, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

full_data <- list()
case_data <- list()
date_data <- list()

# 2. File parsing and distance computation (vectorized)
storm_dir <- "Darwin_2014_POS"
stormnum <- 1
case_data[[as.character(stormnum)]] <- list()
date_data[[as.character(stormnum)]] <- list()

for (shortfile in 1:365) {
  file_path <- paste0(storm_dir,"/darw",shortfile,"0.pos")
  filename <- basename(file_path)
  
  # Read the file header to extract column names
  raw_lines <- readLines(file_path, n = 20)
  header_line <- raw_lines[8] # skip = 7
  
  # 1. Clean the header and remove the "%" and "GPST" keywords
  header_line <- gsub("%|GPST", "", header_line) 
  # 2. Trim excess spaces
  header_line <- trimws(gsub("\\s+", " ", header_line)) 
  
  # 3. Split into remaining column names and drop any accidental empty strings
  data_cols <- unlist(strsplit(header_line, " "))
  data_cols <- data_cols[data_cols != ""] 
  
  # 4. Prepend date and time to perfectly match the 15 data columns
  final_col_names <- c("date", "time", data_cols)
  
  # Read the data block using the perfectly aligned names
  df_raw <- read_table(file_path, skip = 15, col_names = final_col_names, show_col_types = FALSE)
  
  if (nrow(df_raw) == 0) next
  
  mvals <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)

  # Vectorized timestamp parsing
  df <- df_raw %>%
    mutate(
      datetime = ymd_hms(paste(gsub("/", "-", date), time)),
      d = day(datetime),
      m = month(datetime),
      y = year(datetime),
      hr = hour(datetime),
      min = minute(datetime),
      sec = second(datetime),
      doy = mvals[m] + d - 1,
      `total seconds day` = hr * 3600 + min * 60 + sec,
      `total seconds` = doy * 24 * 3600 + `total seconds day`
    )
  
  # Vectorized geodesic distance calculations
  df <- df %>%
    mutate(
      # Point 1 is a simple vector, Point 2 uses cbind() to create an N x 2 matrix
      dy = distGeo(c(ref_coords$lon, ref_coords$lat), 
                    cbind(ref_coords$lon, `latitude(deg)`)),
      
      dx = distGeo(c(ref_coords$lon, ref_coords$lat), 
                    cbind(`longitude(deg)`, ref_coords$lat)),
      
      `absolute height(m)` = abs(ref_coords$height - `height(m)`),
      `absolute latitude(deg)` = abs(dy),
      `absolute longitude(deg)` = abs(dx),
      
      `horizontal(deg)` = distGeo(c(ref_coords$lon, ref_coords$lat), 
                                  cbind(`longitude(deg)`, `latitude(deg)`)),
      
      ang = atan2(dy, dx) / pi * 180,
      `horizontal angle(deg)` = ang,
      `absolute horizontal angle(deg)` = abs(ang),
      `minimum horizontal angle(deg)` = ifelse(abs(ang) > 90, 180 - abs(ang), abs(ang))
    ) %>%
    select(-date, -time, -datetime, -dy, -dx, -ang)
  
  # Store in memory hierarchies
  full_data[[length(full_data) + 1]] <- df
  case_data[[as.character(stormnum)]][[length(case_data[[as.character(stormnum)]]) + 1]] <- df
  date_data[[as.character(stormnum)]][[shortfile]] <- df
  
  print(paste(stormnum, filename, ncol(df), nrow(df)))
}

# Combine lists into dataframes
full_df <- bind_rows(full_data)
write_csv(full_df, "afullframe_Darwin_new.csv")

if ("GPST" %in% names(full_df)) {
  print(length(full_df$GPST))
}

makeplot <- TRUE
printpart <- TRUE

rdim <- list(hr = 0, len = 0, Q1 = 3, Q2 = 3, avg = 2, Q3 = 2, max = 3)
nks <- c("hr", "len", "Q1", "Q2", "avg", "Q3", "max")

# 3. Storm day generation (statistics and plots)
for (stormnum in names(case_data)) {
  storm_df <- bind_rows(case_data[[stormnum]])
  
  # Plotting variables
  maxofs <- 365
  
  # 4. Storm aggregation generation (statistics and plots)
  
  stats_storm <- storm_df %>%
    group_by(hr) %>%
    summarise(
      len = n(),
      Q1 = quantile(`horizontal(deg)`, 0.25, na.rm = TRUE),
      Q2 = median(`horizontal(deg)`, na.rm = TRUE),
      avg = mean(`horizontal(deg)`, na.rm = TRUE),
      Q3 = quantile(`horizontal(deg)`, 0.75, na.rm = TRUE),
      max = max(`horizontal(deg)`, na.rm = TRUE)
    ) %>%
    complete(hr = 0:23, fill = list(len=0, Q1=0, Q2=0, avg=0, Q3=0, max=0))
  
  # LaTeX table print loop (year)
  if (printpart) {
    for (vix in 3:length(nks)) {
      col_name <- nks[vix]
      cat(sprintf("$%s$ $%s$ $%s$ $%s$ $%s$\n", col_name, which.max(stats_storm[[col_name]]) - 1, 
                  round(max(stats_storm[[col_name]]), rdim[[col_name]]), 
                  which.min(stats_storm[[col_name]]) - 1, round(min(stats_storm[[col_name]]), rdim[[col_name]])))
    }
    for (h in 0:23) {
      h1 <- sprintf("%02d", h)
      pvals <- paste0("$", h1, "$")
      for (vix in 3:length(nks)) {
        col_name <- nks[vix]
        val <- stats_storm[[col_name]][h + 1]
        strv <- as.character(round(val, rdim[[col_name]]))
        s1 <- as.character(round(max(stats_storm[[col_name]]), rdim[[col_name]]))
        s2 <- as.character(round(min(stats_storm[[col_name]]), rdim[[col_name]]))
        
        if (strv == s1) pvals <- c(pvals, paste0("$\\mathbf{", strv, "}$"))
        else if (strv == s2) pvals <- c(pvals, paste0("$\\underline{\\mathbf{", strv, "}}$"))
        else pvals <- c(pvals, paste0("$", strv, "$"))
      }
      cat(paste(pvals, collapse = " & "), " \\\\ \\hline\n")
    }
  }
  
  # Make plot (year)
  if (makeplot) {
    secticks <- c(0, 31 * 24 * 3600, 59 * 24 * 3600, 90 * 24 * 3600, 120 * 24 * 3600, 151 * 24 * 3600, 181 * 24 * 3600, 212 * 24 * 3600, 243 * 24 * 3600, 273 * 24 * 3600, 304 * 24 * 3600, 334 * 24 * 3600)
    seclabs <- c("1", "32", "60", "91", "121", "152", "182", "213", "244", "274", "305", "335")
    
    filtval <- 10
    tffilt <- storm_df[storm_df$`horizontal(deg)` < filtval,]
    cat(length(storm_df$`horizontal(deg)`), length(tffilt$`horizontal(deg)`), "\n")
    
    plot_title <- paste0("Time series of horizontal positioning errors [m]\nfor Darwin, NT, Australia (2014)")
    
    p <- ggplot(storm_df, aes(x = `total seconds`, y = `horizontal(deg)`)) +
      geom_line(color = "blue", linewidth = 1) +
      scale_x_continuous(breaks = secticks, labels = seclabs, limits = c(0, max(secticks) + 31 * 24 * 3600), expand = c(0, 0)) +
      labs(title = plot_title, x = "DOY", y = "Horizontal\npositioning\nerrors [m]") +
      bigger_theme
    
    ggsave(paste0("afullframe_Darwin_time_R.pdf"), plot = p, width = 10, height = 4.5, device = cairo_pdf)
    ggsave(paste0("afullframe_Darwin_time_R.png"), plot = p, width = 10, height = 4.5)
    
    stats_stormf <- tffilt %>%
      group_by(hr) %>%
      summarise(
        len = n(),
        Q1 = quantile(`horizontal(deg)`, 0.25, na.rm = TRUE),
        Q2 = median(`horizontal(deg)`, na.rm = TRUE),
        avg = mean(`horizontal(deg)`, na.rm = TRUE),
        Q3 = quantile(`horizontal(deg)`, 0.75, na.rm = TRUE),
        max = max(`horizontal(deg)`, na.rm = TRUE)
      ) %>%
      complete(hr = 0:23, fill = list(len=0, Q1=0, Q2=0, avg=0, Q3=0, max=0))
    
    # LaTeX table print loop (year)
    if (printpart) {
      for (vix in 3:length(nks)) {
        col_name <- nks[vix]
        cat(sprintf("$%s$ $%s$ $%s$ $%s$ $%s$\n", col_name, which.max(stats_stormf[[col_name]]) - 1, 
                    round(max(stats_stormf[[col_name]]), rdim[[col_name]]), 
                    which.min(stats_stormf[[col_name]]) - 1, round(min(stats_stormf[[col_name]]), rdim[[col_name]])))
      }
      for (h in 0:23) {
        h1 <- sprintf("%02d", h)
        pvals <- paste0("$", h1, "$")
        for (vix in 3:length(nks)) {
          col_name <- nks[vix]
          val <- stats_stormf[[col_name]][h + 1]
          strv <- as.character(round(val, rdim[[col_name]]))
          s1 <- as.character(round(max(stats_stormf[[col_name]]), rdim[[col_name]]))
          s2 <- as.character(round(min(stats_stormf[[col_name]]), rdim[[col_name]]))
          
          if (strv == s1) pvals <- c(pvals, paste0("$\\mathbf{", strv, "}$"))
          else if (strv == s2) pvals <- c(pvals, paste0("$\\underline{\\mathbf{", strv, "}}$"))
          else pvals <- c(pvals, paste0("$", strv, "$"))
        }
        cat(paste(pvals, collapse = " & "), " \\\\ \\hline\n")
      }
    }
    
    write_csv(tffilt, "afullframe_Darwin_filter_new.csv")
    
    dff <- read_csv("afullframe_Darwin_filter_new.csv", show_col_types = FALSE)
    
    secticks <- c(0, 31 * 24 * 3600, 59 * 24 * 3600, 90 * 24 * 3600, 120 * 24 * 3600, 151 * 24 * 3600, 181 * 24 * 3600, 212 * 24 * 3600, 243 * 24 * 3600, 273 * 24 * 3600, 304 * 24 * 3600, 334 * 24 * 3600)
    seclabs <- c("1", "32", "60", "91", "121", "152", "182", "213", "244", "274", "305", "335")
    
    plot_title <- paste0("Time series of filtered horizontal positioning errors [m]\n(< ", filtval, ") for Darwin, NT, Australia (2014)")
    
    p <- ggplot(dff, aes(x = `total seconds`, y = `horizontal(deg)`)) +
      geom_line(color = "blue", linewidth = 1) +
      scale_x_continuous(breaks = secticks, labels = seclabs, limits = c(0, max(secticks) + 31 * 24 * 3600), expand = c(0, 0)) +
      labs(title = plot_title, x = "DOY", y = "Horizontal\npositioning\nerrors [m]") +
      bigger_theme
    
    ggsave(paste0("afullframe_Darwin_filter_R.pdf"), plot = p, width = 10, height = 4.5, device = cairo_pdf)
    ggsave(paste0("afullframe_Darwin_filter_R.png"), plot = p, width = 10, height = 4.5)
  }
}