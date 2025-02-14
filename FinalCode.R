#load libraries 
library(tidyr)
library(sf)
library(raster)
library(tmap)
library(tmaptools)
library(spatialEco)
library(RStoolbox)
library(purrr)
library(dplyr)
library(terra)
library(ggplot2)
library(rnaturalearth)
library(lubridate)
library(scales)
library(stringr)
library(viridis)
library(cowplot)

#POT FISHING EFFORT 
#Grab OB data
ob5 <- OBOrig_Proc %>%
  select(YEAR, sector, gear, HAUL_ID, AVG_LAT, AVG_LONG, TOTAL_HOOKS, DMONTH, VESSEL) %>%
  filter(gear == 'Pot') %>%
  drop_na() %>%
  mutate(month = DMONTH) %>%
  group_by(YEAR, month, sector, AVG_LAT, AVG_LONG) %>%
  summarise(Sets = n_distinct(HAUL_ID), Unique_Vessels = n_distinct(VESSEL)) %>%
  ungroup()

#Gab EM data 
em5 <- EMOrig_Proc %>%
  select(YEAR, HAUL_ID, AVG_LAT, AVG_LONG, sector, gear, GearPerSet, SET_DATE, VESSEL) %>%
  filter(gear == 'Pot') %>%
  drop_na() %>%
  mutate(date_col = as.Date(SET_DATE),
         month = month(date_col)) %>%
  group_by(YEAR, month, sector, AVG_LAT, AVG_LONG) %>%
  #uses unique haul ID to generate numhber of sets and number of vessels for every unique combo of year month sector and coordiantes 
  summarise(Sets = n_distinct(HAUL_ID), Unique_Vessels = n_distinct(VESSEL)) %>%
  ungroup()

#combine all pot data into a single dataset 
obs_number_pot_sets <- rbind(ob5, em5)

# Define the resolution (2km in degrees approximately 0.018 degrees, considering 1 degree ~ 111km)
resolution <- 0.018

# Helper function to convert month numbers to month names
month_number_to_name <- function(month_number) {
  month.name[month_number]
}

# Convert month numbers to month names in the original dataframe, make sure variables are numeric
obs_number_pot_sets <- obs_number_pot_sets %>%
  mutate(
    YEAR = as.numeric(YEAR),
    month = as.numeric(month),
    AVG_LAT = as.numeric(AVG_LAT),
    AVG_LONG = as.numeric(AVG_LONG),
    Sets = as.numeric(Sets),
    month_name = month_number_to_name(month),
    Unique_Vessels = as.numeric(Unique_Vessels)
  ) %>%
  filter(YEAR > 2013) %>% #only interested in data from 2014-2023
  #filter out outlier value to be able  match CRS later on 
  filter(AVG_LONG >-670)

# Modify the data frame processing function to track unique vessel counts for confidentiality later on and to regroup data into 2km resolution boxes
create_data_frames <- function(obs_number_pot_sets, resolution) {
  df_list_sets <- list()
  
  unique_combinations <- obs_number_pot_sets %>%
    select(YEAR, month_name, sector) %>%
    distinct()
  
  for (i in 1:nrow(unique_combinations)) {
    year <- unique_combinations$YEAR[i]
    month_name <- unique_combinations$month_name[i]
    sector <- unique_combinations$sector[i]
    
    list_name <- paste(year, month_name, sector, sep = "_")
    
    temp_df <- obs_number_pot_sets %>%
      filter(YEAR == year, month_name == month_name, sector == sector)
    
    if (nrow(temp_df) == 0) {
      print(paste("No data for:", list_name))
      next
    }
    
    # Aggregate by resolution and keep track of unique vessels
    temp_df <- temp_df %>%
      mutate(
        lat_bin = floor(AVG_LAT / resolution) * resolution,
        long_bin = floor(AVG_LONG / resolution) * resolution
      ) %>%
      group_by(sector, month_name, YEAR, lat_bin, long_bin) %>%
      summarise(Sets = sum(Sets, na.rm = TRUE), 
                Unique_Vessels = sum(Unique_Vessels, na.rm = TRUE), .groups = 'drop')
    
    df_list_sets[[list_name]] <- temp_df
  }
  
  return(df_list_sets)
}

# Calling the function
df_list_sets <- create_data_frames(obs_number_pot_sets, resolution)

# Function to filter list of dataframes by exact sector name to seperate sectors 
filter_list_by_exact_name <- function(df_list, sector_name) {
  # Match dataframes by name pattern
  pattern <- paste0("_", sector_name, "$")
  matching_names <- names(df_list)[grepl(pattern, names(df_list))]
  filtered_list <- df_list[matching_names]
  
  # Ensure that each dataframe only contains the specified sector
  filtered_list <- lapply(filtered_list, function(df) {
    df[df$sector == sector_name, , drop = FALSE]  # Filter rows by sector
  })
  
  return(filtered_list)
}

# Filtering by sector
sets_cs <- filter_list_by_exact_name(df_list_sets, "Catch Shares")
sets_csem <- filter_list_by_exact_name(df_list_sets, "Catch Shares EM")
sets_le <- filter_list_by_exact_name(df_list_sets, "Limited Entry Sablefish")
sets_oa <- filter_list_by_exact_name(df_list_sets, "OA Fixed Gear")

# Combine all filtered dataframes
all_sets <- c(sets_cs, sets_csem, sets_le, sets_oa)

# Function to check if the column exists in a dataframe which will be incorporated into split_by_month
column_exists <- function(df, column_name) {
  return(column_name %in% colnames(df))
}

# Function to filter each dataframe by unique month_names and rename them
split_by_month <- function(df_list) {
  new_list <- list()
  
  for (name in names(df_list)) {
    df <- df_list[[name]]
    
    # Check if the 'month_name' column exists
    if (!column_exists(df, "month_name") | !column_exists(df, "YEAR") | !column_exists(df, "sector")) {
      warning(paste("Dataframe", name, "does not contain the required columns: 'month_name', 'YEAR', 'sector'"))
      next
    }
    
    unique_months <- unique(df$month_name)
    
    for (month in unique_months) {
      df_filtered <- df %>% filter(month_name == !!month)
      # Extract the first value of YEAR and sector assuming they are the same for the entire dataframe
      year <- df_filtered$YEAR[1]
      sector <- df_filtered$sector[1]
      new_name <- paste(year, month, sector, sep = "_")
      new_list[[new_name]] <- df_filtered
    }
  }
  
  return(new_list)
}

# Split the combined sets by month
all_sets_split <- split_by_month(all_sets)

# Function to convert month names to numbers for ordering, will be used in next function
month_name_to_number <- function(month_name) {
  match(tolower(month_name), tolower(month.name))
}

# Function to reorganize a list of dataframes by Year, Month, and Sector
reorganize_by_date_and_sector <- function(df_list) {
  # Get the names of the dataframes
  df_names <- names(df_list)
  
  # Split the names into Year, Month, and Sector components
  df_info <- do.call(rbind, strsplit(df_names, "_"))
  colnames(df_info) <- c("Year", "Month", "Sector")
  
  # Convert to a data frame for easier manipulation
  df_info <- as.data.frame(df_info, stringsAsFactors = FALSE)
  
  # Convert Year to numeric
  df_info$Year <- as.numeric(df_info$Year)
  
  # Convert Month to numeric for proper ordering
  df_info$MonthNumber <- month_name_to_number(df_info$Month)
  
  # Order the data frames by Year, Month, and Sector
  order_indices <- order(df_info$Year, df_info$Sector, df_info$MonthNumber)
  
  # Reorder the list based on the sorted indices
  df_list_ordered <- df_list[order_indices]
  
  # Rename the dataframes in the list to reflect their order
  new_names <- df_names[order_indices]
  names(df_list_ordered) <- new_names
  
  return(df_list_ordered)
}

# Reorganize by date and sector (apply as needed)
sets_monthly_df <- reorganize_by_date_and_sector(all_sets_split)

#grab observer coverage file - this was a CSV file I made myself that organized observer coverate rate by year and sector
#this only contained observed coverages for the LE and OA sectors of pot and hook-and-line fisheries 
obs.cov<-read.csv('C:/Users/michaela.melanson/Downloads/GF Observer Coverage - Rfile (3).csv')

#grab pot and rename typo in spreadsheet
obs.cov.pot<-obs.cov %>%
  filter(Gear == 'Pot') %>%
  mutate(Sector = ifelse(Sector == "Limted Entry", "Limited Entry Sablefish", Sector)) %>%
  mutate(Sector = ifelse(Sector == "Catch Shares - EM", "Catch Shares EM", Sector))

# Function to scale up the sets based on observer coverage 
scale_up_sets <- function(df, obs.cov.pot) {
  if (nrow(df) == 0) {
    return(df)
  }
  
  # Extract year and sector from the dataframe
  year <- df$YEAR[1]
  sector <- df$sector[1]
  
  # Find the matching observation rate
  rate <- obs.cov.pot %>%
    filter(Year == year, Sector == sector) %>%
    pull(Rate)
  
  if (length(rate) == 0) {
    warning(paste("No observation rate found for Year:", year, "Sector:", sector))
    rate <- 1
  }
  
  # Scale up the sets
  df <- df %>%
    mutate(Scaled_Sets = Sets / rate)
  
  return(df)
}

#scale up based on observer coverage 
scaled_sets_monthly_df <- lapply(sets_monthly_df, scale_up_sets, obs.cov.pot = obs.cov.pot)

#function to combine list of dataframes into a single dataframe
combine_dataframes <- function(df_list) {
  combined_df <- do.call(rbind, lapply(names(df_list), function(name) {
    df <- df_list[[name]]
    parts <- strsplit(name, "_")[[1]]
    df$Year <- parts[1]
    df$Month <- parts[2]
    df$Sector <- parts[3]
    return(df)
  }))
  return(combined_df)
}

# Combine all dataframes in the list into a single dataframe
scaled_sets_pots <- combine_dataframes(scaled_sets_monthly_df)

#summarise for plot to combine Catch Shares and Catch Shares EM into a single sector and summarise by total and individual sector by year
scaled_sets_pots_summary <- scaled_sets_pots %>%
  mutate(sector = if_else(sector %in% c("Catch Shares", "Catch Shares EM"), "Catch Shares", sector)) %>%
  group_by(Year, sector) %>%
  summarise(Total_Sets = sum(Scaled_Sets), .groups = "drop") %>%
  bind_rows(
    scaled_sets_pots %>%
      group_by(Year) %>%
      summarise(sector = "All Sectors", Total_Sets = sum(Scaled_Sets), .groups = "drop")
  )

# Create the line and dot plot displays estimated pot sets by sector from 2014-2023 
ggplot(scaled_sets_pots_summary, aes(x = Year, y = Total_Sets, color = sector, group=sector)) +
  geom_line() +           # Add lines
  geom_point() +          # Add dots
  labs(title = "Estimated Pot Sets",
       x = "Year", y = "Pot Sets",
       color = 'Sector') +
  scale_color_manual(values = c("Catch Shares" = "red", "Limited Entry Sablefish" = "green", "OA Fixed Gear" = "blue", "All Sectors" = "black"),
                     labels = c("CS" = "Catch Shares", "LE" = "Limited Entry", 
                                "OA" = "Open Access", "Total" = "Total")) +
  theme_minimal() 

#transfrom pot dataframes into rasters for co-occurrence/overlap work and to set CRS for coastline data for mapping
#df to spatraster function to preserve unique vessels and number of sets 
df_to_rast <- function(df, resolution) {
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  # Create the spatial extent based on the lat/long bins
  rast_sets <- rast(xmin = min(df$long_bin) - resolution, xmax = max(df$long_bin) + resolution,
                    ymin = min(df$lat_bin) - resolution, ymax = max(df$lat_bin) + resolution,
                    resolution = resolution)
  rast_vessels <- rast_sets  # Same structure for `Unique_Vessels`
  
  # Initialize both rasters with NA values
  values(rast_sets) <- NA
  values(rast_vessels) <- NA
  
  # Populate the rasters with `Scaled_Sets` and `Unique_Vessels`
  for (j in 1:nrow(df)) {
    cell <- cellFromXY(rast_sets, cbind(df$long_bin[j], df$lat_bin[j]))
    rast_sets[cell] <- df$Scaled_Sets[j]
    rast_vessels[cell] <- df$Unique_Vessels[j]
  }
  
  # Combine both layers into a single multi-layer raster
  combined_rast <- c(rast_sets, rast_vessels)
  names(combined_rast) <- c("Scaled_Sets", "Unique_Vessels")
  
  return(combined_rast)
}

# Convert each dataframe in the list to SpatRaster
scaled_raster_list <- lapply(scaled_sets_monthly_df, df_to_rast, resolution = resolution)
scaled_raster_list <- reorganize_by_date_and_sector(scaled_raster_list)

# Set the names of the rasters in the list to be the same as the dataframes
names(scaled_raster_list) <- names(scaled_sets_monthly_df)

#confidential plot for pot sets on the US West Coast map by year and month 

#combine all scaled up pot effort into a single dataframe 
combined_potsets_year<- bind_rows(scaled_sets_monthly_df, .id='source')

# Adjust the bin size (decimal precision) as needed for your confidentiality rule
bin_coordinates <- function(coord, bin_size) {
  return(round(coord / bin_size) * bin_size)
}

#assigning new bin size which is ~10km
bin_size <- 0.1 

# Apply the binning function to latitude (y) and longitude (x) columns to generate 10km boxes instead of 2km
combined_potsets_year$binned_lon <- bin_coordinates(combined_potsets_year$long_bin, bin_size)
combined_potsets_year$binned_lat <- bin_coordinates(combined_potsets_year$lat_bin, bin_size)

# Group by binned coordinates and summing the amount of unique sets and vessel numbers for each grid point
aggregated_df_values <- combined_potsets_year %>%
  group_by(YEAR, sector, binned_lon, binned_lat) %>%
  summarise(
    Total_Unique_Vessels = sum(Unique_Vessels, na.rm = TRUE),  # Sum up unique vessels
    Sets_Value = sum(Scaled_Sets, na.rm = TRUE)  # Sum the overlapping values (e.g., product)
  ) %>%
  ungroup()

#this is the filter to abide by confidentiality requirements and filters out any grid cells that have 3 or less unqiue vessels 
aggregated_df_values_filtered <- aggregated_df_values %>% filter(Total_Unique_Vessels >= 3)

#for plotting purposes, grabs the unique number of years and how many unique years there are 
unique_years <- unique(aggregated_df_values_filtered$YEAR)
n_years <- length(aggregated_df_values_filtered)

# Set the total width of the plot and how much to shift each year
x_shift <- 4  # Adjust this to control spacing between slivers

# Define constant longitudes to display for every sliver
displayed_longitudes <- c(-125:-116)

# Find global min and max latitudes across all years to standardize the latitude range
global_min_lat <- min(aggregated_df_values_filtered$binned_lat)
global_max_lat <- max(aggregated_df_values_filtered$binned_lat)

# Create a new data frame with shifted coordinates for each year
shifted_data <- do.call(rbind, lapply(seq_along(unique_years), function(i) {
  year <- unique_years[i]
  
  # Filter the data for the specific year
  year_data <- aggregated_df_values_filtered %>%
    filter(YEAR == year) %>%
    mutate(
      # Shift longitudes to create "slivers"
      long_bin_shifted = binned_lon + (i - 1) * x_shift
    )
  
  return(year_data)
}))

# Create a summary dataset with one entry per year for the year label
year_labels <- shifted_data %>%
  group_by(YEAR) %>%
  summarise(
    long_bin_shifted = mean(long_bin_shifted),  # Center the label in the sliver
    lat_bin = global_max_lat + 2,  # Position the label even higher than before
    year_label = as.character(YEAR)
  )

#grab US West Coast coastline information for the plot 
states<- map_data('state', region = c('california','oregon','washington'))
coastline <- ne_download(scale = 10, type = "coastline", category = "physical", returnclass = "sf")
coastline <- st_transform(coastline, crs = st_crs(scaled_raster_list$`2019_August_Catch Shares`))
x_limits<- c(-126,-115)
y_limits<-c(32,49)
states_clipped<- states %>%
  filter(long>=x_limits[1] & long<=x_limits[2])
bbox <- st_bbox(c(xmin = x_limits[1], xmax = x_limits[2], ymin = y_limits[1], ymax = y_limits[2]), crs = st_crs(coastline))
coastline_clipped <- st_crop(coastline, bbox)

# Shift the coastline coordinates as well
shifted_coastline <- do.call(rbind, lapply(seq_along(unique_years), function(i) {
  coastline_shifted <- coastline_clipped %>%
    mutate(
      geometry = st_geometry(coastline_clipped) + c((i - 1) * x_shift, 0)  # Shift the coastline longitude
    )
  
  return(coastline_shifted)
}))

# Plot the coastlines and data by transparency 
ggplot() +
  # Plot the coastlines (shifted)
  geom_sf(data = shifted_coastline, color = "black", fill = NA) +
  
  # Plot the fishing activity (shifted)
  geom_point(data = shifted_data, aes(x = long_bin_shifted, y = binned_lat, alpha = Sets_Value), color = 'red') +
  scale_alpha_continuous(name = "Number of Pot Sets", range = c(0.1, 1)) + # Adjust the transparency range if needed
  
  # Labels for each year (one per sliver, using summary data)
  geom_text(data = year_labels, aes(x = long_bin_shifted, y = lat_bin, label = year_label),
            hjust = 0.5, size = 4, color = "black") +  # Only one label per sliver
  
  # Repeat the same longitude labels for each sliver
  scale_x_continuous(breaks = rep(displayed_longitudes, times = n_years), labels = rep(displayed_longitudes, times = n_years)) +
  
  # Axis and legend settings
  scale_size_continuous(name = "Pot Sets", range = c(1, 4)) +  # Relabel size scale as "Pot Sets"
  labs(x = "Longitude", y = "Latitude", title = "Pot Effort by Year") +
  
  # Standardize the y-axis across all slivers by setting the global min/max latitudes for both the data and labels
  coord_sf(xlim = c(min(shifted_data$long_bin_shifted), max(shifted_data$long_bin_shifted) + x_shift),  # Extend x-axis by x_shift to accommodate last sliver
           ylim = c(global_min_lat, global_max_lat + 2.5),  # Move the year labels even higher by expanding ylim
           expand = FALSE) +  # Remove gray space by avoiding coordinate expansion
  
  guides(alpha = guide_legend(override.aes = list(size = 6))) +  # Make legend points larger
  
  # Remove x-axis numbers for slivers or replace with repeated longitudes
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove original x-axis text
    axis.ticks.x = element_blank(),  # Remove original x-axis ticks
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # Reduce plot margins
    legend.position = "bottom"
  )

#confidentiality plot by month for pot effort

#order months chronologically not alphabetically for the plot 
month_levels <- c('January','February','March','April','May','June','July','August','September','October','November','December')

# Clean the month_name column (trim spaces) and convert to factor with levels in chronological order
combined_potsets_year <- combined_potsets_year %>%
  mutate(month_name = str_trim(month_name),  # Remove leading/trailing spaces
         month_name = factor(month_name, levels = month_levels, ordered = TRUE))  # Reorder correctly

# Reorder the data based on the month_name factor levels
combined_potsets_year <- combined_potsets_year %>%
  arrange(month_name)  # Arrange the data by month_name to ensure proper order

# Manually reorder the unique months based on the factor levels
unique_months <- levels(factor(combined_potsets_year$month_name, levels = month_levels))

# Set the number of slivers (equal to the number of unique months)
n_months <- length(unique_months)

# Group by binned coordinates thus can only plot by value
aggregated_df_values <- combined_potsets_year %>%
  group_by(month_name, sector, binned_lon, binned_lat) %>%
  summarise(
    Total_Unique_Vessels = sum(Unique_Vessels, na.rm = TRUE),  # Sum up unique vessels
    Sets_Value = sum(Scaled_Sets, na.rm = TRUE)  # Sum the overlapping values (e.g., product)
  ) %>%
  ungroup()

#confidentiality filter 
aggregated_df_values_filtered <- aggregated_df_values %>% filter(Total_Unique_Vessels >= 3)

# Set the total width of the plot and how much to shift each year
x_shift <- 4  # Adjust this to control spacing between slivers

# Create a new data frame with shifted coordinates for each month
shifted_data <- do.call(rbind, lapply(seq_along(unique_months), function(i) {
  month <- unique_months[i]
  
  # Filter the data for the specific month
  month_data <- aggregated_df_values_filtered %>%
    filter(month_name == month) %>%
    mutate(
      # Shift longitudes to create "slivers"
      long_bin_shifted = binned_lon + (i - 1) * x_shift
    )
  
  return(month_data)
}))

month_labels <- shifted_data %>%
  group_by(month_name) %>%
  summarise(
    long_bin_shifted = mean(long_bin_shifted),  # Center the label in the sliver
    lat_bin = global_max_lat + 3,  # Position the labels higher for better visibility
    month_label = as.character(month_name)
  )

# Shift the coastline coordinates as well
shifted_coastline <- do.call(rbind, lapply(seq_along(unique_months), function(i) {
  coastline_shifted <- coastline_clipped %>%
    mutate(
      geometry = st_geometry(coastline_clipped) + c((i - 1) * x_shift, 0)  # Shift the coastline longitude
    )
  
  return(coastline_shifted)
}))

# Plot the coastlines and data by transparency 
ggplot() +
  # Plot the coastlines (shifted)
  geom_sf(data = shifted_coastline, color = "black", fill = NA) +
  
  # Plot the fishing activity (shifted)
  geom_point(data = shifted_data, aes(x = long_bin_shifted, y = binned_lat, alpha = Sets_Value), color = 'red') +
  scale_alpha_continuous(name = "Number of Pot Sets", range = c(0.1, 1)) + # Adjust the transparency range if needed
  
  # Labels for each year (one per sliver, using summary data)
  geom_text(data = month_labels, aes(x = long_bin_shifted, y = lat_bin, label = month_label),
            hjust = 0.4, size = 4, color = "black", angle = 45) +  # Only one label per sliver
  
  # Repeat the same longitude labels for each sliver
  #scale_x_continuous(breaks = rep(displayed_longitudes, times = n_months), labels = rep(displayed_longitudes, times = n_months)) +
  
  # Axis and legend settings
  scale_size_continuous(name = "Pot Sets", range = c(1, 4)) +  # Relabel size scale as "Pot Sets"
  labs(x = "Longitude", y = "Latitude", title = "Pot Effort by Month") +
  
  # Standardize the y-axis across all slivers by setting the global min/max latitudes for both the data and labels
  coord_sf(xlim = c(min(shifted_data$long_bin_shifted)-x_shift, max(shifted_data$long_bin_shifted) + x_shift),  # Extend x-axis by x_shift to accommodate last sliver
           ylim = c(global_min_lat, global_max_lat + 6),  # Move the year labels even higher by expanding ylim
           expand = FALSE) +  # Remove gray space by avoiding coordinate expansion
  
  guides(alpha = guide_legend(override.aes = list(size = 6))) +  # Make legend points larger
  
  # Remove x-axis numbers for slivers or replace with repeated longitudes
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove original x-axis text
    axis.ticks.x = element_blank(),  # Remove original x-axis ticks
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # Reduce plot margins
    legend.position = "bottom"
  )

#HOOK-AND-LINE DATA
#Here I will just gather and clean the HKL Data, the above confidentiality plots can also be easily applied to HKL and therefore I will not repeat again in the code for the sake of saving space 
#grab observer data 
ob6 <- OBOrig_Proc %>%
  select(YEAR, sector, gear, HAUL_ID, AVG_LAT, AVG_LONG, TOTAL_HOOKS, DMONTH, VESSEL) %>%
  filter(gear == 'Hook & Line') %>%
  drop_na() %>%
  mutate(month = DMONTH) %>%
  group_by(YEAR, month, sector, AVG_LAT, AVG_LONG) %>%
  summarise(Sets = n_distinct(HAUL_ID), Unique_Vessels = n_distinct(VESSEL)) %>%
  ungroup()

#grab EM data 
em6 <- EMOrig_Proc %>%
  select(YEAR, HAUL_ID, AVG_LAT, AVG_LONG, sector, gear, GearPerSet, SET_DATE, VESSEL) %>%
  filter(gear == 'Hook & Line') %>%
  drop_na() %>%
  mutate(date_col = as.Date(SET_DATE),
         month = month(date_col)) %>%
  group_by(YEAR, month, sector, AVG_LAT, AVG_LONG) %>%
  summarise(Sets = n_distinct(HAUL_ID), Unique_Vessels = n_distinct(VESSEL)) %>%
  ungroup()

#combine EM and OB Data to combine all HKL data 
obs_number_hkl_sets<- rbind(ob6,em6)

# Convert month numbers to month names in the original dataframe, make variables numeric, and count unique vessels 
obs_number_hkl_sets <- obs_number_hkl_sets %>%
  mutate(
    YEAR = as.numeric(YEAR),
    month = as.numeric(month),
    AVG_LAT = as.numeric(AVG_LAT),
    AVG_LONG = as.numeric(AVG_LONG),
    Sets = as.numeric(Sets),
    month_name = month_number_to_name(month),
    Unique_Vessels = as.numeric(Unique_Vessels)
  ) %>%
  filter(YEAR > 2013) #onyl grabbing data from 2014-2023 

# combining all hkl effort into the list of dataframes 
df_list_sets_hkl <- create_data_frames(obs_number_hkl_sets, resolution)

# Filtering by sector
sets_cs <- filter_list_by_exact_name(df_list_sets_hkl, "Catch Shares")
sets_le <- filter_list_by_exact_name(df_list_sets_hkl, "Limited Entry Sablefish")
sets_dtl<- filter_list_by_exact_name(df_list_sets_hkl, 'LE Fixed Gear DTL')
sets_oa <- filter_list_by_exact_name(df_list_sets_hkl, "OA Fixed Gear")

# Combine all filtered dataframes
all_sets_hkl <- c(sets_cs, sets_dtl, sets_le, sets_oa)

# Function to filter each dataframe by unique month_names and rename them
all_sets_split_hkl <- split_by_month(all_sets_hkl)

# reorganization
sets_monthly_df_hkl <- reorganize_by_date_and_sector(all_sets_split_hkl)

#pull HKL observer coverage rates ou tot scale up fishery data 
obs.cov.hkl<-obs.cov %>%
  filter(Gear == 'Hook and Line') %>%
  mutate(Sector = ifelse(Sector == "Limted Entry", "Limited Entry Sablefish", Sector)) %>%
  mutate(Sector = ifelse(Sector == "Catch Shares - EM", "Catch Shares EM", Sector))

scaled_sets_monthly_df_hkl <- lapply(sets_monthly_df_hkl, scale_up_sets, obs.cov.pot = obs.cov.hkl)

# Assuming scaled_sets_monthly_df is your list of dataframes
combined_hklsets_year<- bind_rows(scaled_sets_monthly_df_hkl, .id='source')

#make spatrasters for hkl 
# Convert each dataframe in the list to SpatRaster
scaled_raster_list_hkl <- lapply(scaled_sets_monthly_df_hkl, df_to_rast, resolution = resolution)
scaled_raster_list_hkl <- reorganize_by_date_and_sector(scaled_raster_list_hkl)

# Set the names of the rasters in the list to be the same as the dataframes
names(scaled_raster_list_hkl) <- names(scaled_sets_monthly_df_hkl)

#MIDWATER TRAWL GEAR
#Here I am doing the same thing that I did for HKL, I will gather and clean the data but leave the plots out 
#The only difference for trawl versus pot and hkl is that some of the variables may be slightly different in the plotting commands
#grab observer data 
ob7<- OBOrig_Proc %>%
  select(YEAR,sector,gear,HAUL_ID,HAUL_DURATION,AVG_LAT,AVG_LONG, DMONTH) %>%
  filter(gear=='Midwater Trawl') %>%
  mutate(month = DMONTH,
         HAUL_ID = as.character(HAUL_ID),
         Sets = HAUL_DURATION) %>%
  distinct(HAUL_ID, .keep_all = TRUE) %>%#filters to only keep one row per haul ID
  mutate(Unique_Vessels = 1) #since there is only one row per haul ID right now we can just assign one unqiue vessel for each entry 

#grab EM data 
em7<- EMOrig_Proc %>%
  select(YEAR,HAUL_ID,AVG_LAT,AVG_LONG,sector,gear, SET_DATE, HAUL_DURATION) %>%
  filter(gear=='Midwater Trawl') %>%
  mutate(date_col = as.Date(SET_DATE),
         month= month(date_col), 
         HAUL_ID = as.character(HAUL_ID),
         Sets = HAUL_DURATION) %>%
  distinct(HAUL_ID, .keep_all = TRUE)%>%
  mutate(Unique_Vessels = 1)

#grab at-sea data 
as7 <- ASOrig_Proc %>%
  select(YEAR, HAUL_ID, AVG_LAT, AVG_LONG, sector, gear, DURATION_IN_MIN, DEPLOYMENT_DATE, HAUL_ID) %>%
  filter(gear == 'Midwater Trawl') %>%
  mutate(date_col = as.Date(DEPLOYMENT_DATE),
         month = month(date_col),
         HAUL_ID = as.character(HAUL_ID),
         Sets = DURATION_IN_MIN/60) %>%
  distinct(HAUL_ID, .keep_all = TRUE) %>% # Convert HAUL_ID to character
  mutate(Unique_Vessels = 1)


# Combine all unique columns
all_columns <- union(names(ob7), union(names(em7), names(as7)))


# Add missing columns with NA values
for (col in setdiff(all_columns, names(ob7))) {
  ob7[[col]] <- NA
}

for (col in setdiff(all_columns, names(em7))) {
  em7[[col]] <- NA
}

for (col in setdiff(all_columns, names(as7))) {
  as7[[col]] <- NA
}

# Combining the three datasets
obs_number_trawl_dur <- bind_rows(ob7, em7, as7)

# Convert month numbers to month names in the original dataframe
obs_number_trawl_dur <- obs_number_trawl_dur %>%
  mutate(
    YEAR = as.numeric(YEAR),
    month = as.numeric(month),
    AVG_LAT = as.numeric(AVG_LAT),
    AVG_LONG = as.numeric(AVG_LONG),
    Sets = as.numeric(Sets),
    month_name = month_number_to_name(month),
    Unique_Vessels = as.numeric(Unique_Vessels)
  ) %>%
  filter(YEAR >2013) 

obs_number_trawl_dur <- obs_number_trawl_dur %>%
  filter(Sets > 0)  # Keep only positive values
#
# Modify the data frame processing function to track unique vessel counts
# Function to create a list of data frames grouped by unique combinations of year, month, and sector
create_data_frames <- function(obs_number_pot_sets, resolution) {
  # Initialize an empty list to store the data frames
  df_list_sets <- list()
  
  # Ensure month_name and sector are character columns in the original data
  obs_number_pot_sets <- obs_number_pot_sets %>%
    mutate(
      month_name = as.character(month_name),
      sector = as.character(sector)
    )
  
  # Create unique combinations of YEAR, month_name, and sector
  unique_combinations <- obs_number_pot_sets %>%
    select(YEAR, month_name, sector) %>%
    distinct()
  
  # Loop through each unique combination
  for (i in 1:nrow(unique_combinations)) {
    # Extract the current combination
    year <- unique_combinations$YEAR[i]
    month_name <- unique_combinations$month_name[i]
    sector <- unique_combinations$sector[i]
    
    # Create a name for the list entry
    list_name <- paste(year, month_name, sector, sep = "_")
    
    # Debug: Perform the filtering manually for row-by-row checking
    temp_df <- obs_number_pot_sets[obs_number_pot_sets$YEAR == year &
                                     obs_number_pot_sets$month_name == month_name &
                                     obs_number_pot_sets$sector == sector, ]
    
    # If the filtered dataframe is still incorrect, add extra checks here
    if (nrow(temp_df) == 0) {
      print(paste("No data for:", list_name))
      next
    }
    
    # Aggregate data by the specified resolution and group by the relevant fields
    temp_df <- temp_df %>%
      mutate(
        lat_bin = floor(AVG_LAT / resolution) * resolution,
        long_bin = floor(AVG_LONG / resolution) * resolution) # %>%
    # group_by(sector, month_name, YEAR, lat_bin, long_bin) %>%
    # summarize(Sets = sum(Sets), 
    #  Unique_Vessels = sum(Unique_Vessels, na.rm=TRUE),
    #   .groups = 'drop')
    
    # Store the filtered and aggregated data frame in the list
    df_list_sets[[list_name]] <- temp_df
  }
  
  # Return the list of filtered and aggregated data frames
  return(df_list_sets)
}

# Example of calling the function
df_list_trawl_dur <- create_data_frames(obs_number_trawl_dur, resolution)

# Example usage with the all_sets_split list
sets_monthly_df_trawl <- reorganize_by_date_and_sector(df_list_trawl_dur)

# Define the sectors for at-sea and shoreside
at_sea_sectors <- c("CATCHER-PROCESSOR", "MOTHERSHIP")
shoreside_sectors <- c("Midwater Hake", "Midwater Hake EM", "Midwater Rockfish", "Midwater Rockfish EM")

# Separate the data frames into two lists based on the sector
at_sea <- Filter(function(df_name) {
  any(sapply(at_sea_sectors, function(sector) grepl(sector, df_name)))
}, sets_monthly_df_trawl)

shoreside <- Filter(function(df_name) {
  any(sapply(shoreside_sectors, function(sector) grepl(sector, df_name)))
}, sets_monthly_df_trawl)

# Assuming scaled_sets_monthly_df is your list of dataframes
#combine EM sectors with their human observed counterparts 
sscombined_twlsets_year <- bind_rows(shoreside, .id = 'source')
sscombined_twlsets_year<-sscombined_twlsets_year %>%
  filter(Sets>0) %>%
  mutate(sector = case_when(
    sector == "Midwater Hake EM" ~ "Midwater Hake",
    sector == "Midwater Rockfish EM" ~ "Midwater Rockfish",
    TRUE ~ sector
  ))

#at-sea
ascombined_twlsets_year <- bind_rows(at_sea, .id = 'source')
ascombined_twlsets_year<-ascombined_twlsets_year %>%
  filter(Sets>0)

#make spatrasters for trawl effort 
df_to_rast <- function(df, resolution) {
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  # Create the spatial extent based on the lat/long bins
  rast_sets <- rast(xmin = min(df$long_bin) - resolution, xmax = max(df$long_bin) + resolution,
                    ymin = min(df$lat_bin) - resolution, ymax = max(df$lat_bin) + resolution,
                    resolution = resolution)
  rast_vessels <- rast_sets  # Same structure for `Unique_Vessels`
  
  # Initialize both rasters with NA values
  values(rast_sets) <- NA
  values(rast_vessels) <- NA
  
  # Populate the rasters with `Scaled_Sets` and `Unique_Vessels`
  for (j in 1:nrow(df)) {
    cell <- cellFromXY(rast_sets, cbind(df$long_bin[j], df$lat_bin[j]))
    rast_sets[cell] <- df$Sets[j]
    rast_vessels[cell] <- df$Unique_Vessels[j]
  }
  
  # Combine both layers into a single multi-layer raster
  combined_rast <- c(rast_sets, rast_vessels)
  names(combined_rast) <- c("Sets", "Unique_Vessels")
  
  return(combined_rast)
}

#OVERLAP SS
shoreside_effort<- combine_dataframes(shoreside)
atsea_effort<- combine_dataframes(at_sea)
# Convert each dataframe in the list to SpatRaster
scaled_raster_list_trawl <- lapply(shoreside, df_to_rast, resolution = resolution)
scaled_raster_list_trawl <- reorganize_by_date_and_sector(scaled_raster_list_trawl)
# Set the names of the rasters in the list to be the same as the dataframes
names(scaled_raster_list_trawl) <- names(shoreside)

#AS sector:
scaled_raster_list_trawlas <- lapply(at_sea, df_to_rast, resolution = resolution)
scaled_raster_list_trawlas <- reorganize_by_date_and_sector(scaled_raster_list_trawlas)
# Set the names of the rasters in the list to be the same as the dataframes
names(scaled_raster_list_trawlas) <- names(at_sea)

#SPECIES DISTRIBUTION MODELS 
#HUMPBACKS 
# Function to create raster stacks and calculate mean by month for each year
create_monthly_rasters <- function(base_path) {
  # List all subdirectories for each year
  year_dirs <- list.dirs(path = base_path, full.names = TRUE, recursive = FALSE)
  
  # Create an empty list to store the rasters
  monthly_rasters <- list()
  
  # Loop through each year directory
  for (year_dir in year_dirs) {
    year <- basename(year_dir)
    
    # List all subdirectories for each month within the year
    month_dirs <- list.dirs(path = year_dir, full.names = TRUE, recursive = FALSE)
    
    for (month_dir in month_dirs) {
      month <- basename(month_dir)
      
      # List all .grd files in the month directory
      file_paths <- list.files(path = month_dir, pattern = "\\.grd$", full.names = TRUE)
      
      if (length(file_paths) > 0) {
        # Read the raster files
        rasters <- lapply(file_paths, rast)
        
        # Stack the rasters
        raster_stack <- rast(rasters)
        
        # Calculate the mean of the raster stack
        monthly_mean <- mean(raster_stack, na.rm = TRUE)
        
        # Store the result in the list with a key formatted as "Year_Month"
        key <- paste0(year, "_", month)
        monthly_rasters[[key]] <- monthly_mean
      }
    }
  }
  
  return(monthly_rasters)
}

# Define the base path
#this is after I downloaded all of the daily .grd and .gri files that karin provided to me and placed into a folder
base_path <- 'C:/Users/michaela.melanson/Desktop/R Directory/groundfish biop/spatial analysis/humpback/karin'

# Create the monthly rasters
monthly_rasters <- create_monthly_rasters(base_path)

# Rename SpatRaster objects with month names
month_names <- c("01" = "January", "02" = "February", "03" = "March", "04" = "April", 
                 "05" = "May", "06" = "June", "07" = "July", "08" = "August", 
                 "09" = "September", "10" = "October", "11" = "November", "12" = "December")

rename_spat_rasters <- function(raster_list, month_names) {
  new_raster_list <- list()
  
  for (name in names(raster_list)) {
    name_parts <- strsplit(name, "_")[[1]]
    year <- name_parts[1]
    month <- name_parts[2]
    
    # Get the corresponding month name
    month_name <- month_names[month]
    
    # Create new name
    new_name <- paste(year, month_name, sep = "_")
    
    # Add the raster to the new list with the new name
    new_raster_list[[new_name]] <- raster_list[[name]]
  }
  
  return(new_raster_list)
}

# Rename the SpatRaster list
humpback_monthly <- rename_spat_rasters(monthly_rasters, month_names)
humpback.monthly<-humpback_monthly[37:156] #only selecting months from 2014-2023 

#plotting humpback SDMS
# Define a function to convert a list of SpatRasters to a list of data frames
spatrasters_to_dfs <- function(spat_raster_list) {
  df_list <- lapply(spat_raster_list, function(spat_raster) {
    as.data.frame(spat_raster, xy = TRUE)
  })
  return(df_list)
}

# transform spatrasters into a dataframe 
humpbackmonthdf<-  spatrasters_to_dfs(humpback.monthly)
##3-panel month plots for all years 

# Define the years and months we are interested in
years <- 2014:2023
months_of_interest <- c("April", "July", "October")

# Combine the list of dataframes into a single dataframe
combined_df <- do.call(rbind, lapply(names(humpbackmonthdf), function(name) {
  df <- humpbackmonthdf[[name]]
  
  # Extract year and month from the name (assumes "year_month" format)
  parts <- unlist(strsplit(name, "_"))
  year <- as.numeric(parts[1])
  month <- parts[2]
  
  # Add the year and month columns to the dataframe
  df$year <- year
  df$month <- month
  
  return(df)
}))

# Convert the month to a factor with levels in the correct order
combined_df <- combined_df %>%
  mutate(month = factor(month, levels = months_of_interest, ordered = TRUE))

# Filter the dataframe for the months we are interested in (April, July, September)
filtered_df <- combined_df %>% filter(month %in% months_of_interest)

# Check if the directory exists, and if not, create it
output_dir <- "C:/Users/michaela.melanson/Desktop/Groundfish Plots/rough draft edits/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through each year and generate/save a plot for the months of interest
for (year in years) {
  
  # Filter data for the specific year and the months of interest
  filtered_df_year <- filtered_df %>% filter(year == !!year)
  
  # Create the plot
  p_year <- ggplot() +
    # Add the whale density data
    geom_raster(data = filtered_df_year, aes(x = x, y = y, fill = mean)) +
    
    # Add the coastline layer
    geom_sf(data = coastline_clipped, color = "black", fill = NA) +
    
    # Use viridis color scale for whale density
    scale_fill_viridis(name = "Whale Density") +
    
    # Facet grid for the three months (April, July, September)
    facet_grid(~ month) +
    
    # Apply the correct coordinate system (assumes WGS84: EPSG 4326)
    coord_sf(crs = st_crs(4326)) +
    
    # Theme settings
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10),  # Adjust facet label size
      axis.title = element_blank(),          # Remove axis titles
      axis.text = element_blank(),           # Remove axis text
      axis.ticks = element_blank(),          # Remove axis ticks
      legend.position = "bottom"             # Place legend at bottom
    ) +
    
    # Add title
    ggtitle(paste("Humpback Whale Distribution for April, July, and October - Year", year))
  
  # Define the file path
  file_path <- paste0(output_dir, "humpback_whale_distribution_", year, ".png")
  
  # Print the file path for debugging
  print(paste("Saving plot to:", file_path))
  
  # Save the plot to the specified file
  ggsave(filename = file_path, plot = p_year, width = 10, height = 6)
}

#AVG MONTH PLOTS BY YEARS FOR HUMPBACKS 
# Define the months in the correct order for plotting
months_in_order <- c("January", "February", "March", "April", "May", "June", 
                     "July", "August", "September", "October", "November", "December")

# Combine the list of dataframes into a single dataframe
combined_df <- do.call(rbind, lapply(names(humpbackmonthdf), function(name) {
  df <- humpbackmonthdf[[name]]
  
  # Extract year and month from the name (assumes "year_month" format)
  parts <- unlist(strsplit(name, "_"))
  year <- as.numeric(parts[1])
  month <- parts[2]
  
  # Add the year and month columns to the dataframe
  df$year <- year
  df$month <- month
  
  return(df)
}))

# Convert month to a factor with levels in the correct order
combined_df <- combined_df %>%
  mutate(month = factor(month, levels = months_in_order, ordered = TRUE))

# Calculate the average whale density for each month across all years
avg_df <- combined_df %>%
  group_by(month, x, y) %>%  # Group by month and coordinates
  summarise(mean_density = mean(mean, na.rm = TRUE))  # Calculate the average whale density

# Plot the average whale density for each month in a 3x4 grid
p <- ggplot() +
  # Add the average whale density data
  geom_raster(data = avg_df, aes(x = x, y = y, fill = mean_density)) +
  
  # Add the coastline layer
  geom_sf(data = coastline_clipped, color = "black", fill = NA) +
  
  # Use viridis color scale for whale density
  scale_fill_viridis(name = "Avg Whale Density") +
  
  # Create a 3x4 grid for each month
  facet_wrap(~ month, ncol = 4) +
  
  # Apply the correct coordinate system (assumes WGS84: EPSG 4326)
  coord_sf(crs = st_crs(4326)) +
  
  # Theme settings
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),  # Adjust facet label size
    axis.title = element_blank(),          # Remove axis titles
    axis.text = element_blank(),           # Remove axis text
    axis.ticks = element_blank(),          # Remove axis ticks
    legend.position = "bottom"             # Place legend at bottom
  ) +
  
  # Add title
  ggtitle("Average Humpback Whale Density Across Months (2014-2023)")

# Save the plot to a file
output_dir <- "C:/Users/michaela.melanson/Desktop/Groundfish Plots/rough draft edits/"
ggsave(filename = paste0(output_dir, "average_humpback_whale_density_monthly_letter.png"),
       plot = p, 
       width = 8.5,      # US Letter width in inches
       height = 11,      # US Letter height in inches
       units = "in",     
       dpi = 300)

#LEATHERBACKS 
create_monthly_rasters <- function(base_path, pattern) {
  # List all .grd files in the specified path
  file_paths <- list.files(path = base_path, pattern = pattern, full.names = TRUE)
  
  # Extract year and month from file names
  file_info <- data.frame(
    file_path = file_paths,
    year = as.numeric(sub(".*_(\\d{4})-(\\d{2})-\\d{2}_mean\\.grd$", "\\1", file_paths)),
    month = as.numeric(sub(".*_(\\d{4})-(\\d{2})-\\d{2}_mean\\.grd$", "\\2", file_paths)),
    stringsAsFactors = FALSE
  )
  
  # Create an empty list to store the rasters
  monthly_rasters <- list()
  
  # Loop through each year and month to create raster stacks and calculate mean
  for (year in unique(file_info$year)) {
    for (month in unique(file_info$month)) {
      # Filter files for the specific year and month
      files_to_stack <- file_info %>%
        filter(year == !!year, month == !!month) %>%
        pull(file_path)
      
      if (length(files_to_stack) > 0) {
        # Read the raster files
        rasters <- lapply(files_to_stack, rast)
        
        # Stack the rasters
        raster_stack <- rast(rasters)
        
        # Calculate the mean of the raster stack
        monthly_mean <- mean(raster_stack, na.rm = TRUE)
        
        # Store the result in the list with a key formatted as "Year_Month"
        key <- paste0(year, "_", sprintf("%02d", month))
        monthly_rasters[[key]] <- monthly_mean
      }
    }
  }
  
  return(monthly_rasters)
}
# Define the base path and file pattern
#again this is after I have downloaded the daily .grd and .gri files into a folder on my computer 
base_path <- 'C:/Users/michaela.melanson/Desktop/R Directory/groundfish biop/spatial analysis/leatherback/allyears'
pattern <- 'lbst_\\d{4}-\\d{2}-\\d{2}_mean\\.grd$' # Adjust pattern to match your file naming convention

# Create the monthly rasters
monthly_rasters <- create_monthly_rasters(base_path, pattern)

#filter monthly rasters 2013 and up 
# Define the year and month thresholds
year_threshold <- 2013
month_threshold <- 0
# Function to filter the list based on year and month
filter_spat_rasters <- function(raster_list, year_threshold, month_threshold) {
  filtered_list <- list()
  
  for (name in names(raster_list)) {
    # Split the name to get year and month
    name_parts <- strsplit(name, "_")[[1]]
    year <- as.numeric(name_parts[1])
    month <- as.numeric(name_parts[2])
    
    # Filter condition: year > 2013 or (year == 2013 and month > month_threshold)
    if (year > year_threshold || (year > year_threshold && month > month_threshold)) {
      filtered_list[[name]] <- raster_list[[name]]
    }
  }
  
  return(filtered_list)
}

lbst_monthly <- filter_spat_rasters(monthly_rasters, year_threshold, month_threshold)

month_names <- c("January", "February", "March", "April", "May", "June", 
                 "July", "August", "September", "October", "November", "December")

# Function to rename SpatRaster objects
rename_spat_rasters <- function(raster_list, month_names) {
  new_raster_list <- list()
  
  for (name in names(raster_list)) {
    # Split the name to get year and month
    name_parts <- strsplit(name, "_")[[1]]
    year <- name_parts[1]
    month <- as.numeric(name_parts[2])
    
    # Get the corresponding month name
    month_name <- month_names[month]
    
    # Create new name
    new_name <- paste(year, month_name, sep = "_")
    
    # Add the raster to the new list with the new name
    new_raster_list[[new_name]] <- raster_list[[name]]
  }
  
  return(new_raster_list)
}

# Rename the SpatRaster list
lbst.monthly<- rename_spat_rasters(lbst_monthly, month_names)

#transform leatherback spatrasters into a dataframe 
lbstmonthdf<-  spatrasters_to_dfs(lbst.monthly)

# Define the years and months we are interested in
years <- 2014:2023
months_of_interest <- c("May", "August", "November")

# Combine the list of dataframes into a single dataframe for lbstmonthdf (Leatherback Sea Turtle)
combined_df_lbst <- do.call(rbind, lapply(names(lbstmonthdf), function(name) {
  df <- lbstmonthdf[[name]]  # Use lbstmonthdf instead of humpbackmonthdf
  
  # Extract year and month from the name (assumes "year_month" format)
  parts <- unlist(strsplit(name, "_"))
  year <- as.numeric(parts[1])
  month <- parts[2]
  
  # Add the year and month columns to the dataframe
  df$year <- year
  df$month <- month
  
  return(df)
}))

# Convert the month to a factor with levels in the correct order
combined_df_lbst <- combined_df_lbst %>%
  mutate(month = factor(month, levels = months_of_interest, ordered = TRUE))

# Filter the dataframe for the months we are interested in (April, July, September)
filtered_df_lbst <- combined_df_lbst %>% filter(month %in% months_of_interest)

# Directory to save the plots
output_dir <- "C:/Users/michaela.melanson/Desktop/Groundfish Plots/rough draft edits/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through each year and generate/save a plot for the months of interest
for (year in years) {
  
  # Filter data for the specific year and the months of interest
  filtered_df_year <- filtered_df_lbst %>% filter(year == !!year)
  
  # Create the plot
  p_year <- ggplot() +
    # Add the habitat suitability data
    geom_raster(data = filtered_df_year, aes(x = x, y = y, fill = mean)) +  
    
    # Add the coastline layer
    geom_sf(data = coastline_clipped, color = "black", fill = NA) +
    
    # Use viridis color scale for habitat suitability
    scale_fill_viridis(name = "Habitat Suitability") +
    
    # Facet grid for the three months (April, July, September)
    facet_grid(~ month) +
    
    # Apply the correct coordinate system (assumes WGS84: EPSG 4326)
    coord_sf(crs = st_crs(4326)) +
    
    # Theme settings
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10),  # Adjust facet label size
      axis.title = element_blank(),          # Remove axis titles
      axis.text = element_blank(),           # Remove axis text
      axis.ticks = element_blank(),          # Remove axis ticks
      legend.position = "bottom"             # Place legend at bottom
    ) +
    
    # Add title
    ggtitle(paste("Leatherback Sea Turtle Habitat Suitability for May, August, and November - Year", year))
  
  # Save the plot to a file
  ggsave(filename = paste0(output_dir, "lbst_habitat_suitability_", year, ".png"),
         plot = p_year, width = 10, height = 6)
}


#MONTHLY LBST ACROSS YEARS
# Define the months in the correct order for plotting
months_in_order <- c("January", "February", "March", "April", "May", "June", 
                     "July", "August", "September", "October", "November", "December")

# Combine the list of dataframes into a single dataframe
combined_df_lbst <- do.call(rbind, lapply(names(lbstmonthdf), function(name) {
  df <- lbstmonthdf[[name]]  # Use lbstmonthdf instead
  
  # Extract year and month from the name (assumes "year_month" format)
  parts <- unlist(strsplit(name, "_"))
  year <- as.numeric(parts[1])
  month <- parts[2]
  
  # Add the year and month columns to the dataframe
  df$year <- year
  df$month <- month
  
  return(df)
}))

# Convert month to a factor with levels in the correct order
combined_df_lbst <- combined_df_lbst %>%
  mutate(month = factor(month, levels = months_in_order, ordered = TRUE))

# Calculate the average habitat suitability for each month across all years
avg_df_lbst <- combined_df_lbst %>%
  group_by(month, x, y) %>%  # Group by month and coordinates
  summarise(mean_suitability = mean(mean, na.rm = TRUE)) 

# Plot the average habitat suitability for each month in a 3x4 grid
p <- ggplot() +
  # Add the average habitat suitability data
  geom_raster(data = avg_df_lbst, aes(x = x, y = y, fill = mean_suitability)) +
  
  # Add the coastline layer
  geom_sf(data = coastline_clipped, color = "black", fill = NA) +
  
  # Use viridis color scale for habitat suitability
  scale_fill_viridis(name = "Avg Habitat Suitability") +
  
  # Create a 3x4 grid for each month
  facet_wrap(~ month, ncol = 4) +
  
  # Apply the correct coordinate system (assumes WGS84: EPSG 4326)
  coord_sf(crs = st_crs(4326)) +
  
  # Theme settings
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),  # Adjust facet label size
    axis.title = element_blank(),          # Remove axis titles
    axis.text = element_blank(),           # Remove axis text
    axis.ticks = element_blank(),          # Remove axis ticks
    legend.position = "bottom"             # Place legend at bottom
  ) +
  
  # Add title
  ggtitle("Average Leatherback Sea Turtle Habitat Suitability Across Months (2014-2023)")

# Save the plot to a file
output_dir <- "C:/Users/michaela.melanson/Desktop/Groundfish Plots/rough draft edits/"
ggsave(filename = paste0(output_dir, "average_lbst_habitat_suitability_monthly_letter.png"),
       plot = p, 
       width = 8.5,      # US Letter width in inches
       height = 11,      # US Letter height in inches
       units = "in",     
       dpi = 300)


#OVERLAP ANALYSIS 
#HUMPBACK OVERLAP FOR ALL GEAR TYPES 

#POT OVERLAP WITH HUMPBACKS 
process_and_plot <- function(year, month, sector, spat_raster, lbst_raster, resolution) {
  if (!is.null(spat_raster)) {
    
    # Check and log CRS mismatch, reproject if needed for both layers in `spat_raster`
    if (!identical(crs(spat_raster), crs(lbst_raster))) {
      print(paste("CRS mismatch: Reprojecting spat_raster for", year, month, sector))
      spat_raster <- project(spat_raster, crs(lbst_raster))
    }
    
    # Check for extent overlap
    ext_overlap <- intersect(ext(spat_raster), ext(lbst_raster))
    if (is.null(ext_overlap)) {
      warning(paste("No spatial extent overlap for", year, month, sector))
      return(NULL)
    } else {
      print(paste("Extent overlap found for", year, month, sector))
    }
    
    # Separate layers in the `spat_raster`
    scaled_sets_rast <- spat_raster[[1]]  # Fishing effort layer
    unique_vessels_rast <- spat_raster[[2]]  # Unique vessels layer
    
    # Reproject `unique_vessels_rast` independently if needed
    if (!identical(crs(unique_vessels_rast), crs(lbst_raster))) {
      print(paste("CRS mismatch: Reprojecting unique_vessels_rast for", year, month, sector))
      unique_vessels_rast <- project(unique_vessels_rast, crs(lbst_raster))
    }
    
    # Resample `scaled_sets_rast` and `unique_vessels_rast` to match the habitat suitability raster
    scaled_sets_rast <- resample(scaled_sets_rast, lbst_raster, method = "bilinear")
    unique_vessels_rast <- resample(unique_vessels_rast, lbst_raster, method = "near")
    
    # Mask both rasters to keep only overlapping areas with `lbst_raster`
    overlap_raster <- mask(scaled_sets_rast, lbst_raster)
    unique_vessels_rast <- mask(unique_vessels_rast, lbst_raster)
    
    valid_overlap_cells <- sum(!is.na(values(overlap_raster)))
    valid_unique_vessel_cells <- sum(!is.na(values(unique_vessels_rast)))
    
    # print(paste("Overlap raster has", valid_overlap_cells, "valid cells"))
    # print(paste("Unique vessels raster has", valid_unique_vessel_cells, "valid cells"))
    
    if (valid_overlap_cells == 0 || valid_unique_vessel_cells == 0) {
      warning(paste("No valid data for", year, month, sector))
      return(NULL)
    }
    
    # Convert the habitat suitability `lbst_raster` to a dataframe as well
    lbst_df <- as.data.frame(lbst_raster, xy = TRUE, na.rm = TRUE)
    
    # Convert both rasters (Scaled Sets and Unique Vessels) to dataframes
    product_df <- as.data.frame(overlap_raster, xy = TRUE, na.rm = TRUE)
    unique_vessels_df <- as.data.frame(unique_vessels_rast, xy = TRUE, na.rm = TRUE)
    
    # Merge the three dataframes based on spatial coordinates (x, y)
    combined_df <- merge(product_df, unique_vessels_df, by = c("x", "y"), all.x = TRUE)
    combined_df <- merge(combined_df, lbst_df, by = c("x", "y"), all.x = TRUE, suffixes = c("_ScaledSets", "_LBST"))
    
    return(combined_df)
  } else {
    warning(paste("spat_raster is NULL for", year, month, sector))
    return(NULL)
  }
}
# Create a list to store the dataframes
product_set_month <- list()

# Loop through the year, month, and sector combinations
for (name in names(scaled_raster_list)) {
  parts <- strsplit(name, "_")[[1]]
  year <- parts[1]
  month <- parts[2]
  sector <- parts[3]
  
  spat_raster_name <- paste(year, month, sector, sep = "_")
  spat_raster <- scaled_raster_list[[spat_raster_name]]
  lbst_raster_name <- paste(year, month, sep = "_")
  lbst_raster <- humpback.monthly[[lbst_raster_name]]
  
  product_df <- process_and_plot(year, month, sector, spat_raster, lbst_raster, resolution)
  
  # Check if the result is not NULL before storing it
  if (!is.null(product_df)) {
    product_set_month[[spat_raster_name]] <- product_df
  } else {
    warning(paste("No data for", spat_raster_name))
  }
}

# Filter overlap dataframe for NA values 
has_rows <- function(df) {
  return(nrow(df) > 0)
}

# Filter the list of dataframes to keep only those with rows
filtered_product_monthly <- Filter(has_rows, product_set_month)

# Combine all dataframes in the list into a single dataframe
combined_filtered_df <- combine_dataframes(filtered_product_monthly)

conf.hbw.pot<-combined_filtered_df %>%
  mutate(product = (mean*Scaled_Sets))

# Apply the binning function to latitude (y) and longitude (x) columns to 10 km
conf.hbw.pot$binned_lon <- bin_coordinates(conf.hbw.pot$x, bin_size)
conf.hbw.pot$binned_lat <- bin_coordinates(conf.hbw.pot$y, bin_size)

# Group by binned coordinates thus can only plot by value
aggregated_df_values <- conf.hbw.pot %>%
  group_by(binned_lon, binned_lat) %>%
  summarise(
    Total_Unique_Vessels = sum(Unique_Vessels, na.rm = TRUE),  # Sum up unique vessels
    Overlap_Value = sum(product, na.rm = TRUE)  # Sum the overlapping values (e.g., product)
  ) %>%
  ungroup()

#filtering for confidentiality so no point displayed shows less than 3 vessels 
aggregated_df_values_filtered <- aggregated_df_values %>% filter(Total_Unique_Vessels >= 3)

#set coordinate limits for plots
y_limits=c(32,49)
x_limits=c(-126,-117)

#plot to display overlap on the map for leatherbacks and pot fishing, taking into consideration confidentiality
ggplot(data = aggregated_df_values_filtered) +
  geom_tile(aes(x = binned_lon, y = binned_lat, fill = Overlap_Value)) +
  geom_sf(data = coastline_clipped, color = "black", fill = NA) +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_minimal() +
  coord_sf(xlim = x_limits, ylim = y_limits, crs = st_crs(coastline_clipped)) +
  labs(title = "Pot Sets Overlap with Humpback Whale Density",
       x = "Longitude", y = "Latitude", fill = "Overlap Value") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(floor(min(x_limits)), ceiling(max(x_limits)), by = 1)) +
  scale_y_continuous(breaks = seq(floor(min(y_limits)), ceiling(max(y_limits)), by = 1))

#summary tables for overlap by month and sector and year and sector (for sumamry tables and line plots)
hbw.pot.year.sum <- conf.hbw.pot %>%
  mutate(Sector = if_else(Sector %in% c("Catch Shares", "Catch Shares EM"), "Catch Shares", Sector)) %>%
  group_by(Year, Sector) %>%
  summarise(Total_Overlap = sum(product))

hbw.pot.month.sum <- conf.hbw.pot %>%
  mutate(Sector = if_else(Sector %in% c("Catch Shares", "Catch Shares EM"), "Catch Shares", Sector)) %>%
  group_by(Month, Sector) %>%
  summarise(Total_Overlap = sum(product))

#HKL OVERLAP WITH HUMPBACKS 
# Create a list to store the dataframes
product_set_month <- list()

# Loop through the year, month, and sector combinations
for (name in names(scaled_raster_list_hkl)) {
  parts <- strsplit(name, "_")[[1]]
  year <- parts[1]
  month <- parts[2]
  sector <- parts[3]
  
  spat_raster_name <- paste(year, month, sector, sep = "_")
  spat_raster <- scaled_raster_list_hkl[[spat_raster_name]]
  lbst_raster_name <- paste(year, month, sep = "_")
  lbst_raster <- humpback.monthly[[lbst_raster_name]]
  
  product_df <- process_and_plot(year, month, sector, spat_raster, lbst_raster, resolution)
  
  # Check if the result is not NULL before storing it
  if (!is.null(product_df)) {
    product_set_month[[spat_raster_name]] <- product_df
  } else {
    warning(paste("No data for", spat_raster_name))
  }
}

# Filter the list of dataframes to keep only those with rows
filtered_product_monthly <- Filter(has_rows, product_set_month)

# Combine all dataframes in the list into a single dataframe
combined_filtered_df <- combine_dataframes(filtered_product_monthly)

conf.hbw.hkl<-combined_filtered_df %>%
  mutate(product = (mean*Scaled_Sets))

# Apply the binning function to latitude (y) and longitude (x) columns to 10 km
conf.hbw.hkl$binned_lon <- bin_coordinates(conf.hbw.hkl$x, bin_size)
conf.hbw.hkl$binned_lat <- bin_coordinates(conf.hbw.hkl$y, bin_size)

# Group by binned coordinates thus can only plot by value
aggregated_df_values <- conf.hbw.hkl %>%
  group_by(binned_lon, binned_lat) %>%
  summarise(
    Total_Unique_Vessels = sum(Unique_Vessels, na.rm = TRUE),  # Sum up unique vessels
    Overlap_Value = sum(product, na.rm = TRUE)  # Sum the overlapping values (e.g., product)
  ) %>%
  ungroup()

#filtering for confidentiality so no point displayed shows less than 3 vessels 
aggregated_df_values_filtered <- aggregated_df_values %>% filter(Total_Unique_Vessels >= 3)

#set coordinate limits for plots
y_limits=c(32,49)
x_limits=c(-126,-117)

#plot to display overlap on the map for leatherbacks and pot fishing, taking into consideration confidentiality
ggplot(data = aggregated_df_values_filtered) +
  geom_tile(aes(x = binned_lon, y = binned_lat, fill = Overlap_Value)) +
  geom_sf(data = coastline_clipped, color = "black", fill = NA) +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_minimal() +
  coord_sf(xlim = x_limits, ylim = y_limits, crs = st_crs(coastline_clipped)) +
  labs(title = "Hook-and-Line Sets Overlap with Humpback Whale Density",
       x = "Longitude", y = "Latitude", fill = "Overlap Value") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(floor(min(x_limits)), ceiling(max(x_limits)), by = 1)) +
  scale_y_continuous(breaks = seq(floor(min(y_limits)), ceiling(max(y_limits)), by = 1))

#summary tables for overlap by month and sector and year and sector (for sumamry tables and line plots)
hbw.hkl.year.sum <- conf.hbw.hkl %>%
  group_by(Year, Sector) %>%
  summarise(Total_Overlap = sum(product))

hbw.hkl.month.sum <- conf.hbw.hkl %>%
  group_by(Month, Sector) %>%
  summarise(Total_Overlap = sum(product))

#SHORESIDE TRAWL OVERLAP WITH HUMPBACKS
#new function, nothing fundamentally changes just makes is accessible to work with the trawl hours data instead of sets 
process_and_plot <- function(year, month, sector, spat_raster, lbst_raster, resolution) {
  if (!is.null(spat_raster)) {
    
    # Check and log CRS mismatch, reproject if needed for both layers in `spat_raster`
    if (!identical(crs(spat_raster), crs(lbst_raster))) {
      print(paste("CRS mismatch: Reprojecting spat_raster for", year, month, sector))
      spat_raster <- project(spat_raster, crs(lbst_raster))
    }
    
    # Check for extent overlap
    ext_overlap <- intersect(ext(spat_raster), ext(lbst_raster))
    if (is.null(ext_overlap)) {
      warning(paste("No spatial extent overlap for", year, month, sector))
      return(NULL)
    } else {
      print(paste("Extent overlap found for", year, month, sector))
    }
    
    # Separate layers in the `spat_raster`
    scaled_sets_rast <- spat_raster[[1]]  # Fishing effort layer
    unique_vessels_rast <- spat_raster[[2]]  # Unique vessels layer
    
    # Reproject `unique_vessels_rast` independently if needed
    if (!identical(crs(unique_vessels_rast), crs(lbst_raster))) {
      print(paste("CRS mismatch: Reprojecting unique_vessels_rast for", year, month, sector))
      unique_vessels_rast <- project(unique_vessels_rast, crs(lbst_raster))
    }
    
    # Resample `scaled_sets_rast` and `unique_vessels_rast` to match the habitat suitability raster
    scaled_sets_rast <- resample(scaled_sets_rast, lbst_raster, method = "bilinear")
    unique_vessels_rast <- resample(unique_vessels_rast, lbst_raster, method = "near")
    
    # Mask both rasters to keep only overlapping areas with `lbst_raster`
    overlap_raster <- mask(scaled_sets_rast, lbst_raster)
    unique_vessels_rast <- mask(unique_vessels_rast, lbst_raster)
    
    valid_overlap_cells <- sum(!is.na(values(overlap_raster)))
    valid_unique_vessel_cells <- sum(!is.na(values(unique_vessels_rast)))
    
    if (valid_overlap_cells == 0 || valid_unique_vessel_cells == 0) {
      warning(paste("No valid data for", year, month, sector))
      return(NULL)
    }
    
    # Convert the habitat suitability `lbst_raster` to a dataframe as well
    lbst_df <- as.data.frame(lbst_raster, xy = TRUE, na.rm = TRUE)
    
    # Convert both rasters (Scaled Sets and Unique Vessels) to dataframes
    product_df <- as.data.frame(overlap_raster, xy = TRUE, na.rm = TRUE)
    unique_vessels_df <- as.data.frame(unique_vessels_rast, xy = TRUE, na.rm = TRUE)
    
    # Merge the three dataframes based on spatial coordinates (x, y)
    combined_df <- merge(product_df, unique_vessels_df, by = c("x", "y"), all.x = TRUE)
    combined_df <- merge(combined_df, lbst_df, by = c("x", "y"), all.x = TRUE, suffixes = c("_Sets", "_LBST"))
    
    # Log the structure of the combined dataframe
    # print(paste("Combined DF Rows:", nrow(combined_df)))
    
    return(combined_df)
  } else {
    warning(paste("spat_raster is NULL for", year, month, sector))
    return(NULL)
  }
}
# Create a list to store the dataframes
product_set_month <- list()

# Loop through the year, month, and sector combinations
for (name in names(scaled_raster_list_trawl)) {
  parts <- strsplit(name, "_")[[1]]
  year <- parts[1]
  month <- parts[2]
  sector <- parts[3]
  
  spat_raster_name <- paste(year, month, sector, sep = "_")
  spat_raster <- scaled_raster_list_trawl[[spat_raster_name]]
  lbst_raster_name <- paste(year, month, sep = "_")
  lbst_raster <- humpback.monthly[[lbst_raster_name]]
  
  product_df <- process_and_plot(year, month, sector, spat_raster, lbst_raster, resolution)
  
  # Check if the result is not NULL before storing it
  if (!is.null(product_df)) {
    product_set_month[[spat_raster_name]] <- product_df
  } else {
    warning(paste("No data for", spat_raster_name))
  }
}

filtered_product_monthlyss <- Filter(has_rows, product_set_month)

# Combine all dataframes in the list into a single dataframe
combined_filtered_ss <- combine_dataframes(filtered_product_monthlyss)

conf.hbw.trawlss<-combined_filtered_ss %>%
  mutate(product = (mean*Sets))

# Apply the binning function to latitude (y) and longitude (x) columns
conf.hbw.trawlss$binned_lon <- bin_coordinates(conf.hbw.trawlss$x, bin_size)
conf.hbw.trawlss$binned_lat <- bin_coordinates(conf.hbw.trawlss$y, bin_size)

#getting new values for grid cells 
aggregated_df_values <- conf.hbw.trawlss %>%
  group_by(binned_lon, binned_lat) %>%
  summarise(
    Total_Unique_Vessels = sum(Unique_Vessels, na.rm = TRUE),  # Sum up unique vessels
    Overlap_Value = sum(product, na.rm = TRUE)  # Sum the overlapping values (e.g., product)
  ) %>%
  ungroup()

#filtering for confidentiality 
aggregated_df_values_filtered <- aggregated_df_values %>% filter(Total_Unique_Vessels >= 3)

#setting x and y limits for the plot 
y_limits=c(40,49)
x_limits=c(-126,-122)

#plotting on a map the overlap between leatherbacks and ss trawl effort 
ggplot(data = aggregated_df_values_filtered) +
  geom_tile(aes(x = binned_lon, y = binned_lat, fill = Overlap_Value)) +
  geom_sf(data = coastline_clipped, color = "black", fill = NA) +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_minimal() +
  coord_sf(xlim = x_limits, ylim = y_limits, crs = st_crs(coastline_clipped)) +
  labs(title = " SS Midwater Trawl Overlap with Humpback Whale Density",
       x = "Longitude", y = "Latitude", fill = "Overlap Value") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(floor(min(x_limits)), ceiling(max(x_limits)), by = 1)) +
  scale_y_continuous(breaks = seq(floor(min(y_limits)), ceiling(max(y_limits)), by = 1))


hbw.ss.year.sum<-conf.hbw.trawlss %>%
  group_by(Year, Sector) %>%
  summarise(prod = sum(product))%>%
  pivot_wider(names_from = Sector, values_from = prod)
write.csv(hbw.ss.year.sum, 'C:/Users/michaela.melanson/Downloads/year.sectorss2.csv') 

hbw.ss.month.sum<-conf.hbw.trawlss%>%
  group_by(Month,Sector) %>%
  summarise(prod = sum(product)) %>%
  pivot_wider(names_from = Sector, values_from = prod)
write.csv(hbw.ss.month.sum, 'C:/Users/michaela.melanson/Downloads/month.sectorss2.csv') 

# For combined_filtered_ss dataset
shoreside_proportions_hbw <- conf.hbw.trawlss %>%
  mutate(Latitude_Group = ifelse(y < 46, "Below 46", "Above 46")) %>%  # Categorize latitudes
  group_by(Sector, Latitude_Group) %>%  # Group by Sector and Latitude group
  summarize(Total_Overlap = sum(product, na.rm = TRUE)) %>%  # Sum the Sets column within each group
  mutate(Proportion = Total_Overlap / sum(Total_Overlap)) %>%  # Calculate proportion within each sector
  ungroup() %>%
  mutate(
    Overall_Overlap = sum(Total_Overlap),  # Calculate total overlap across all groups
    Overall_Proportion = Total_Overlap / Overall_Overlap  # Calculate the proportion across all Latitude_Groups
  )
#AT-SEA TRAWL WITH HUMPBACKS 
# Create a list to store the dataframes
product_set_month <- list()

# Loop through the year, month, and sector combinations
for (name in names(scaled_raster_list_trawlas)) {
  parts <- strsplit(name, "_")[[1]]
  year <- parts[1]
  month <- parts[2]
  sector <- parts[3]
  
  spat_raster_name <- paste(year, month, sector, sep = "_")
  spat_raster <- scaled_raster_list_trawlas[[spat_raster_name]]
  lbst_raster_name <- paste(year, month, sep = "_")
  lbst_raster <- humpback.monthly[[lbst_raster_name]]
  
  product_df <- process_and_plot(year, month, sector, spat_raster, lbst_raster, resolution)
  
  # Check if the result is not NULL before storing it
  if (!is.null(product_df)) {
    product_set_month[[spat_raster_name]] <- product_df
  } else {
    warning(paste("No data for", spat_raster_name))
  }
}

filtered_product_monthlyas <- Filter(has_rows, product_set_month)

# Combine all dataframes in the list into a single dataframe
combined_filtered_as <- combine_dataframes(filtered_product_monthlyas)

conf.hbw.trawlas<-combined_filtered_as %>%
  mutate(product = (mean*Sets))

# Apply the binning function to latitude (y) and longitude (x) columns
conf.hbw.trawlas$binned_lon <- bin_coordinates(conf.hbw.trawlas$x, bin_size)
conf.hbw.trawlas$binned_lat <- bin_coordinates(conf.hbw.trawlas$y, bin_size)

#getting new values for grid cells 
aggregated_df_values <- conf.hbw.trawlas %>%
  group_by(binned_lon, binned_lat) %>%
  summarise(
    Total_Unique_Vessels = sum(Unique_Vessels, na.rm = TRUE),  # Sum up unique vessels
    Overlap_Value = sum(product, na.rm = TRUE)  # Sum the overlapping values (e.g., product)
  ) %>%
  ungroup()

#filtering for confidentiality 
aggregated_df_values_filtered <- aggregated_df_values %>% filter(Total_Unique_Vessels >= 3)

#setting x and y limits for the plot 
y_limits=c(40,49)
x_limits=c(-126,-122)

#plotting on a map the overlap between leatherbacks and ss trawl effort 
ggplot(data = aggregated_df_values_filtered) +
  geom_tile(aes(x = binned_lon, y = binned_lat, fill = Overlap_Value)) +
  geom_sf(data = coastline_clipped, color = "black", fill = NA) +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_minimal() +
  coord_sf(xlim = x_limits, ylim = y_limits, crs = st_crs(coastline_clipped)) +
  labs(title = " AS Midwater Trawl Overlap with Humpback Whale Density",
       x = "Longitude", y = "Latitude", fill = "Overlap Value") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(floor(min(x_limits)), ceiling(max(x_limits)), by = 1)) +
  scale_y_continuous(breaks = seq(floor(min(y_limits)), ceiling(max(y_limits)), by = 1))

hbw.as.year.sum<-conf.hbw.trawlas %>%
  group_by(Year, Sector) %>%
  summarise(prod = sum(product))%>%
  pivot_wider(names_from = Sector, values_from = prod)
write.csv(hbw.as.year.sum, 'C:/Users/michaela.melanson/Downloads/hbw.as.year.sum.csv') 

hbw.as.month.sum<-conf.hbw.trawlas%>%
  group_by(Month,Sector) %>%
  summarise(prod = sum(product)) %>%
  pivot_wider(names_from = Sector, values_from = prod)
write.csv(hbw.as.year.sum, 'C:/Users/michaela.melanson/Downloads/hbw.as.year.sum.csv') 

# For combined_filtered_ss dataset
atsea_proportions_hbw <- conf.hbw.trawlas %>%
  mutate(Latitude_Group = ifelse(y < 46, "Below 46", "Above 46")) %>%  # Categorize latitudes
  group_by(Sector, Latitude_Group) %>%  # Group by Sector and Latitude group
  summarize(Total_Overlap = sum(product, na.rm = TRUE)) %>%  # Sum the Sets column within each group
  mutate(Proportion = Total_Overlap / sum(Total_Overlap)) %>%  # Calculate proportion within each sector
  ungroup() %>%
  mutate(
    Overall_Overlap = sum(Total_Overlap),  # Calculate total overlap across all groups
    Overall_Proportion = Total_Overlap / Overall_Overlap  # Calculate the proportion across all Latitude_Groups
  )
#AT-SEA TRAWL OVERLAP WITH HUMPBACKS 

#LEATHERBACK OVERLAP FOR ALL GEAR TYPES 
#POT OVERLAP WITH LEATHERBACKS
process_and_plot <- function(year, month, sector, spat_raster, lbst_raster, resolution) {
  if (!is.null(spat_raster)) {
    
    # Check and log CRS mismatch, reproject if needed for both layers in `spat_raster`
    if (!identical(crs(spat_raster), crs(lbst_raster))) {
      print(paste("CRS mismatch: Reprojecting spat_raster for", year, month, sector))
      spat_raster <- project(spat_raster, crs(lbst_raster))
    }
    
    # Check for extent overlap
    ext_overlap <- intersect(ext(spat_raster), ext(lbst_raster))
    if (is.null(ext_overlap)) {
      warning(paste("No spatial extent overlap for", year, month, sector))
      return(NULL)
    } else {
      print(paste("Extent overlap found for", year, month, sector))
    }
    
    # Separate layers in the `spat_raster`
    scaled_sets_rast <- spat_raster[[1]]  # Fishing effort layer
    unique_vessels_rast <- spat_raster[[2]]  # Unique vessels layer
    
    # Reproject `unique_vessels_rast` independently if needed
    if (!identical(crs(unique_vessels_rast), crs(lbst_raster))) {
      print(paste("CRS mismatch: Reprojecting unique_vessels_rast for", year, month, sector))
      unique_vessels_rast <- project(unique_vessels_rast, crs(lbst_raster))
    }
    
    # Resample `scaled_sets_rast` and `unique_vessels_rast` to match the habitat suitability raster
    scaled_sets_rast <- resample(scaled_sets_rast, lbst_raster, method = "bilinear")
    unique_vessels_rast <- resample(unique_vessels_rast, lbst_raster, method = "near")
    
    # Mask both rasters to keep only overlapping areas with `lbst_raster`
    overlap_raster <- mask(scaled_sets_rast, lbst_raster)
    unique_vessels_rast <- mask(unique_vessels_rast, lbst_raster)
    
    valid_overlap_cells <- sum(!is.na(values(overlap_raster)))
    valid_unique_vessel_cells <- sum(!is.na(values(unique_vessels_rast)))
    
    # print(paste("Overlap raster has", valid_overlap_cells, "valid cells"))
    # print(paste("Unique vessels raster has", valid_unique_vessel_cells, "valid cells"))
    
    if (valid_overlap_cells == 0 || valid_unique_vessel_cells == 0) {
      warning(paste("No valid data for", year, month, sector))
      return(NULL)
    }
    
    # Convert the habitat suitability `lbst_raster` to a dataframe as well
    lbst_df <- as.data.frame(lbst_raster, xy = TRUE, na.rm = TRUE)
    
    # Convert both rasters (Scaled Sets and Unique Vessels) to dataframes
    product_df <- as.data.frame(overlap_raster, xy = TRUE, na.rm = TRUE)
    unique_vessels_df <- as.data.frame(unique_vessels_rast, xy = TRUE, na.rm = TRUE)
    
    # Merge the three dataframes based on spatial coordinates (x, y)
    combined_df <- merge(product_df, unique_vessels_df, by = c("x", "y"), all.x = TRUE)
    combined_df <- merge(combined_df, lbst_df, by = c("x", "y"), all.x = TRUE, suffixes = c("_ScaledSets", "_LBST"))
    
    return(combined_df)
  } else {
    warning(paste("spat_raster is NULL for", year, month, sector))
    return(NULL)
  }
}
# Create a list to store the dataframes
product_set_month <- list()

# Loop through the year, month, and sector combinations
for (name in names(scaled_raster_list)) {
  parts <- strsplit(name, "_")[[1]]
  year <- parts[1]
  month <- parts[2]
  sector <- parts[3]
  
  spat_raster_name <- paste(year, month, sector, sep = "_")
  spat_raster <- scaled_raster_list[[spat_raster_name]]
  lbst_raster_name <- paste(year, month, sep = "_")
  lbst_raster <- lbst.monthly[[lbst_raster_name]]
  
  product_df <- process_and_plot(year, month, sector, spat_raster, lbst_raster, resolution)
  
  # Check if the result is not NULL before storing it
  if (!is.null(product_df)) {
    product_set_month[[spat_raster_name]] <- product_df
  } else {
    warning(paste("No data for", spat_raster_name))
  }
}

# Filter the list of dataframes to keep only those with rows
filtered_product_monthly <- Filter(has_rows, product_set_month)

# Combine all dataframes in the list into a single dataframe
combined_filtered_df <- combine_dataframes(filtered_product_monthly)

conf.lbst.pot<-combined_filtered_df %>%
  mutate(product = (mean*Scaled_Sets))

# Apply the binning function to latitude (y) and longitude (x) columns to 10 km
conf.lbst.pot$binned_lon <- bin_coordinates(conf.lbst.pot$x, bin_size)
conf.lbst.pot$binned_lat <- bin_coordinates(conf.lbst.pot$y, bin_size)

# Group by binned coordinates thus can only plot by value
aggregated_df_values <- conf.lbst.pot %>%
  group_by(binned_lon, binned_lat) %>%
  summarise(
    Total_Unique_Vessels = sum(Unique_Vessels, na.rm = TRUE),  # Sum up unique vessels
    Overlap_Value = sum(product, na.rm = TRUE)  # Sum the overlapping values (e.g., product)
  ) %>%
  ungroup()

#filtering for confidentiality so no point displayed shows less than 3 vessels 
aggregated_df_values_filtered <- aggregated_df_values %>% filter(Total_Unique_Vessels >= 3)

#set coordinate limits for plots
y_limits=c(32,49)
x_limits=c(-126,-117)

#plot to display overlap on the map for leatherbacks and pot fishing, taking into consideration confidentiality
ggplot(data = aggregated_df_values_filtered) +
  geom_tile(aes(x = binned_lon, y = binned_lat, fill = Overlap_Value)) +
  geom_sf(data = coastline_clipped, color = "black", fill = NA) +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_minimal() +
  coord_sf(xlim = x_limits, ylim = y_limits, crs = st_crs(coastline_clipped)) +
  labs(title = "Pot Sets Overlap with Leatherback Habitat Suitability",
       x = "Longitude", y = "Latitude", fill = "Overlap Value") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(floor(min(x_limits)), ceiling(max(x_limits)), by = 1)) +
  scale_y_continuous(breaks = seq(floor(min(y_limits)), ceiling(max(y_limits)), by = 1))

#summary tables for overlap by month and sector and year and sector (for sumamry tables and line plots)
lbst.pot.year.sum <- conf.lbst.pot %>%
  mutate(Sector = if_else(Sector %in% c("Catch Shares", "Catch Shares EM"), "Catch Shares", Sector)) %>%
  group_by(Year, Sector) %>%
  summarise(Total_Overlap = sum(product))

lbst.pot.month.sum <- conf.lbst.pot %>%
  mutate(Sector = if_else(Sector %in% c("Catch Shares", "Catch Shares EM"), "Catch Shares", Sector)) %>%
  group_by(Month, Sector) %>%
  summarise(Total_Overlap = sum(product))

#HKL OVERLAP WITH LEATHERBACKS
# Create a list to store the dataframes
product_set_month <- list()

# Loop through the year, month, and sector combinations
for (name in names(scaled_raster_list_hkl)) {
  parts <- strsplit(name, "_")[[1]]
  year <- parts[1]
  month <- parts[2]
  sector <- parts[3]
  
  spat_raster_name <- paste(year, month, sector, sep = "_")
  spat_raster <- scaled_raster_list_hkl[[spat_raster_name]]
  lbst_raster_name <- paste(year, month, sep = "_")
  lbst_raster <- lbst.monthly[[lbst_raster_name]]
  
  product_df <- process_and_plot(year, month, sector, spat_raster, lbst_raster, resolution)
  
  # Check if the result is not NULL before storing it
  if (!is.null(product_df)) {
    product_set_month[[spat_raster_name]] <- product_df
  } else {
    warning(paste("No data for", spat_raster_name))
  }
}

# Filter the list of dataframes to keep only those with rows
filtered_product_monthly <- Filter(has_rows, product_set_month)

# Combine all dataframes in the list into a single dataframe
combined_filtered_df <- combine_dataframes(filtered_product_monthly)

conf.lbst.hkl<-combined_filtered_df %>%
  mutate(product = (mean*Scaled_Sets))

# Apply the binning function to latitude (y) and longitude (x) columns to 10 km
conf.lbst.hkl$binned_lon <- bin_coordinates(conf.lbst.hkl$x, bin_size)
conf.lbst.hkl$binned_lat <- bin_coordinates(conf.lbst.hkl$y, bin_size)

# Group by binned coordinates thus can only plot by value
aggregated_df_values <- conf.lbst.hkl %>%
  group_by(binned_lon, binned_lat) %>%
  summarise(
    Total_Unique_Vessels = sum(Unique_Vessels, na.rm = TRUE),  # Sum up unique vessels
    Overlap_Value = sum(product, na.rm = TRUE)  # Sum the overlapping values (e.g., product)
  ) %>%
  ungroup()

#filtering for confidentiality so no point displayed shows less than 3 vessels 
aggregated_df_values_filtered <- aggregated_df_values %>% filter(Total_Unique_Vessels >= 3)

#set coordinate limits for plots
y_limits=c(32,49)
x_limits=c(-126,-117)

#plot to display overlap on the map for leatherbacks and pot fishing, taking into consideration confidentiality
ggplot(data = aggregated_df_values_filtered) +
  geom_tile(aes(x = binned_lon, y = binned_lat, fill = Overlap_Value)) +
  geom_sf(data = coastline_clipped, color = "black", fill = NA) +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_minimal() +
  coord_sf(xlim = x_limits, ylim = y_limits, crs = st_crs(coastline_clipped)) +
  labs(title = "Hook-and-Line Sets Overlap with Leatherback Habitat Suitability",
       x = "Longitude", y = "Latitude", fill = "Overlap Value") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(floor(min(x_limits)), ceiling(max(x_limits)), by = 1)) +
  scale_y_continuous(breaks = seq(floor(min(y_limits)), ceiling(max(y_limits)), by = 1))

#summary tables for overlap by month and sector and year and sector (for sumamry tables and line plots)
lbst.hkl.year.sum <- conf.lbst.hkl %>%
  group_by(Year, Sector) %>%
  summarise(Total_Overlap = sum(product))

lbst.hkl.month.sum <- conf.lbst.hkl %>%
  group_by(Month, Sector) %>%
  summarise(Total_Overlap = sum(product))

#SHORESIDE TRAWL WITH LEATHERBACKS 
#new function, nothing fundamentally changes just makes is accessible to work with the trawl hours data instead of sets 
process_and_plot <- function(year, month, sector, spat_raster, lbst_raster, resolution) {
  if (!is.null(spat_raster)) {
    
    # Check and log CRS mismatch, reproject if needed for both layers in `spat_raster`
    if (!identical(crs(spat_raster), crs(lbst_raster))) {
      print(paste("CRS mismatch: Reprojecting spat_raster for", year, month, sector))
      spat_raster <- project(spat_raster, crs(lbst_raster))
    }
    
    # Check for extent overlap
    ext_overlap <- intersect(ext(spat_raster), ext(lbst_raster))
    if (is.null(ext_overlap)) {
      warning(paste("No spatial extent overlap for", year, month, sector))
      return(NULL)
    } else {
      print(paste("Extent overlap found for", year, month, sector))
    }
    
    # Separate layers in the `spat_raster`
    scaled_sets_rast <- spat_raster[[1]]  # Fishing effort layer
    unique_vessels_rast <- spat_raster[[2]]  # Unique vessels layer
    
    # Reproject `unique_vessels_rast` independently if needed
    if (!identical(crs(unique_vessels_rast), crs(lbst_raster))) {
      print(paste("CRS mismatch: Reprojecting unique_vessels_rast for", year, month, sector))
      unique_vessels_rast <- project(unique_vessels_rast, crs(lbst_raster))
    }
    
    # Resample `scaled_sets_rast` and `unique_vessels_rast` to match the habitat suitability raster
    scaled_sets_rast <- resample(scaled_sets_rast, lbst_raster, method = "bilinear")
    unique_vessels_rast <- resample(unique_vessels_rast, lbst_raster, method = "near")
    
    # Mask both rasters to keep only overlapping areas with `lbst_raster`
    overlap_raster <- mask(scaled_sets_rast, lbst_raster)
    unique_vessels_rast <- mask(unique_vessels_rast, lbst_raster)
    
    valid_overlap_cells <- sum(!is.na(values(overlap_raster)))
    valid_unique_vessel_cells <- sum(!is.na(values(unique_vessels_rast)))
    
    #print(paste("Overlap raster has", valid_overlap_cells, "valid cells"))
    #print(paste("Unique vessels raster has", valid_unique_vessel_cells, "valid cells"))
    
    if (valid_overlap_cells == 0 || valid_unique_vessel_cells == 0) {
      warning(paste("No valid data for", year, month, sector))
      return(NULL)
    }
    
    # Convert the habitat suitability `lbst_raster` to a dataframe as well
    lbst_df <- as.data.frame(lbst_raster, xy = TRUE, na.rm = TRUE)
    
    # Convert both rasters (Scaled Sets and Unique Vessels) to dataframes
    product_df <- as.data.frame(overlap_raster, xy = TRUE, na.rm = TRUE)
    unique_vessels_df <- as.data.frame(unique_vessels_rast, xy = TRUE, na.rm = TRUE)
    
    # Log the row counts for all dataframes
    # print(paste("Product DF Rows:", nrow(product_df), "Unique Vessels DF Rows:", nrow(unique_vessels_df), "LBST DF Rows:", nrow(lbst_df)))
    
    # Merge the three dataframes based on spatial coordinates (x, y)
    combined_df <- merge(product_df, unique_vessels_df, by = c("x", "y"), all.x = TRUE)
    combined_df <- merge(combined_df, lbst_df, by = c("x", "y"), all.x = TRUE, suffixes = c("_Sets", "_LBST"))
    
    # Log the structure of the combined dataframe
    # print(paste("Combined DF Rows:", nrow(combined_df)))
    
    return(combined_df)
  } else {
    warning(paste("spat_raster is NULL for", year, month, sector))
    return(NULL)
  }
}
# Create a list to store the dataframes
product_set_month <- list()

# Loop through the year, month, and sector combinations
for (name in names(scaled_raster_list_trawl)) {
  parts <- strsplit(name, "_")[[1]]
  year <- parts[1]
  month <- parts[2]
  sector <- parts[3]
  
  spat_raster_name <- paste(year, month, sector, sep = "_")
  spat_raster <- scaled_raster_list_trawl[[spat_raster_name]]
  lbst_raster_name <- paste(year, month, sep = "_")
  lbst_raster <- lbst.monthly[[lbst_raster_name]]
  
  product_df <- process_and_plot(year, month, sector, spat_raster, lbst_raster, resolution)
  
  # Check if the result is not NULL before storing it
  if (!is.null(product_df)) {
    product_set_month[[spat_raster_name]] <- product_df
  } else {
    warning(paste("No data for", spat_raster_name))
  }
}

filtered_product_monthlyss <- Filter(has_rows, product_set_month)

# Combine all dataframes in the list into a single dataframe
combined_filtered_ss <- combine_dataframes(filtered_product_monthlyss)

conf.lbst.trawlss<-combined_filtered_ss %>%
  mutate(product = (mean*Sets))

# Apply the binning function to latitude (y) and longitude (x) columns
conf.lbst.trawlss$binned_lon <- bin_coordinates(conf.lbst.trawlss$x, bin_size)
conf.lbst.trawlss$binned_lat <- bin_coordinates(conf.lbst.trawlss$y, bin_size)

#getting new values for grid cells 
aggregated_df_values <- conf.lbst.trawlss %>%
  group_by(binned_lon, binned_lat) %>%
  summarise(
    Total_Unique_Vessels = sum(Unique_Vessels, na.rm = TRUE),  # Sum up unique vessels
    Overlap_Value = sum(product, na.rm = TRUE)  # Sum the overlapping values (e.g., product)
  ) %>%
  ungroup()

#filtering for confidentiality 
aggregated_df_values_filtered <- aggregated_df_values %>% filter(Total_Unique_Vessels >= 3)

#setting x and y limits for the plot 
y_limits=c(40,49)
x_limits=c(-126,-122)

#plotting on a map the overlap between leatherbacks and ss trawl effort 
ggplot(data = aggregated_df_values_filtered) +
  geom_tile(aes(x = binned_lon, y = binned_lat, fill = Overlap_Value)) +
  geom_sf(data = coastline_clipped, color = "black", fill = NA) +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_minimal() +
  coord_sf(xlim = x_limits, ylim = y_limits, crs = st_crs(coastline_clipped)) +
  labs(title = " SS Midwater Trawl Overlap with Leatherback Habitat Sutiability",
       x = "Longitude", y = "Latitude", fill = "Overlap Value") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(floor(min(x_limits)), ceiling(max(x_limits)), by = 1)) +
  scale_y_continuous(breaks = seq(floor(min(y_limits)), ceiling(max(y_limits)), by = 1))


lbst.ss.year.sum<-conf.lbst.trawlss %>%
  group_by(Year, Sector) %>%
  summarise(prod = sum(product))%>%
  pivot_wider(names_from = Sector, values_from = prod)
write.csv(lbst.ss.year.sum, 'C:/Users/michaela.melanson/Downloads/lbst.ss.year.sum.csv') 

lbst.ss.month.sum<-conf.hmw.trawlss%>%
  group_by(Month,Sector) %>%
  summarise(prod = sum(product)) %>%
  pivot_wider(names_from = Sector, values_from = prod)
write.csv(lbst.ss.month.sum, 'C:/Users/michaela.melanson/Downloads/lbst.ss.month.sum.csv') 

#AT-SEA TRAWL WITH LEATHERBACKS 
# Create a list to store the dataframes
product_set_month <- list()

# Loop through the year, month, and sector combinations
for (name in names(scaled_raster_list_trawlas)) {
  parts <- strsplit(name, "_")[[1]]
  year <- parts[1]
  month <- parts[2]
  sector <- parts[3]
  
  spat_raster_name <- paste(year, month, sector, sep = "_")
  spat_raster <- scaled_raster_list_trawlas[[spat_raster_name]]
  lbst_raster_name <- paste(year, month, sep = "_")
  lbst_raster <- lbst.monthly[[lbst_raster_name]]
  
  product_df <- process_and_plot(year, month, sector, spat_raster, lbst_raster, resolution)
  
  # Check if the result is not NULL before storing it
  if (!is.null(product_df)) {
    product_set_month[[spat_raster_name]] <- product_df
  } else {
    warning(paste("No data for", spat_raster_name))
  }
}

filtered_product_monthlyas <- Filter(has_rows, product_set_month)

# Combine all dataframes in the list into a single dataframe
combined_filtered_as <- combine_dataframes(filtered_product_monthlyas)

conf.lbst.trawlas<-combined_filtered_as %>%
  mutate(product = (mean*Sets))

# Apply the binning function to latitude (y) and longitude (x) columns
conf.lbst.trawlas$binned_lon <- bin_coordinates(conf.lbst.trawlas$x, bin_size)
conf.lbst.trawlas$binned_lat <- bin_coordinates(conf.lbst.trawlas$y, bin_size)

#getting new values for grid cells 
aggregated_df_values <- conf.lbst.trawlas %>%
  group_by(binned_lon, binned_lat) %>%
  summarise(
    Total_Unique_Vessels = sum(Unique_Vessels, na.rm = TRUE),  # Sum up unique vessels
    Overlap_Value = sum(product, na.rm = TRUE)  # Sum the overlapping values (e.g., product)
  ) %>%
  ungroup()

#filtering for confidentiality 
aggregated_df_values_filtered <- aggregated_df_values %>% filter(Total_Unique_Vessels >= 3)

#setting x and y limits for the plot 
y_limits=c(40,49)
x_limits=c(-126,-122)

#plotting on a map the overlap between leatherbacks and ss trawl effort 
ggplot(data = aggregated_df_values_filtered) +
  geom_tile(aes(x = binned_lon, y = binned_lat, fill = Overlap_Value)) +
  geom_sf(data = coastline_clipped, color = "black", fill = NA) +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_minimal() +
  coord_sf(xlim = x_limits, ylim = y_limits, crs = st_crs(coastline_clipped)) +
  labs(title = " AS Midwater Trawl Overlap with Leatherback Habitat Sutiability",
       x = "Longitude", y = "Latitude", fill = "Overlap Value") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(floor(min(x_limits)), ceiling(max(x_limits)), by = 1)) +
  scale_y_continuous(breaks = seq(floor(min(y_limits)), ceiling(max(y_limits)), by = 1))


lbst.as.year.sum<-conf.lbst.trawlas %>%
  group_by(Year, Sector) %>%
  summarise(prod = sum(product))%>%
  pivot_wider(names_from = Sector, values_from = prod)
write.csv(lbst.as.year.sum, 'C:/Users/michaela.melanson/Downloads/lbst.as.year.sum.csv') 

lbst.as.month.sum<-conf.hmw.trawlas%>%
  group_by(Month,Sector) %>%
  summarise(prod = sum(product)) %>%
  pivot_wider(names_from = Sector, values_from = prod)
write.csv(lbst.as.month.sum, 'C:/Users/michaela.melanson/Downloads/lbst.as.month.sum.csv') 

#AMENDMENT 32 ANALYSIS 
#HUMPBACKS 
# Compute the average habitat suitability for each month
compute_monthly_average <- function(month, data_list) {
  # Filter the list for the specified month
  month_dfs <- data_list[grep(month, names(data_list))]
  
  # Bind all dataframes for the same month together
  combined_df <- bind_rows(month_dfs)
  
  # Group by coordinates and calculate the mean habitat suitability
  averaged_df <- combined_df %>%
    group_by(x, y) %>%
    summarise(mean_suitability = mean(mean, na.rm = TRUE))
  
  return(averaged_df)
}

# List of months in the dataset
months <- c("January", "February", "March", "April", "May", "June", 
            "July", "August", "September", "October", "November", "December")

# Compute the average humpback whale density for each month across years (2014-2023)
monthly_avg_habitat <- map(months, ~ compute_monthly_average(.x, humpbackmonthdf))

# Name the list elements
names(monthly_avg_habitat) <- months

#this is a dataframe that I created that has fathoms associated with coordinate values along the west coast 
fathoms<-read.csv('fathom.csv')

# Convert the fathoms dataframe to an sf object
fathoms_sf <- st_as_sf(fathoms, coords = c("Longitude", "Latitude"), crs = 4326)

# Filter using the geometry column to check if coordinates fall within the specified ranges
#this is the northern california/oregon area 
fathom_area_1 <- fathoms_sf %>%
  filter(
    st_coordinates(.)[,2] >= 40.17 & st_coordinates(.)[,2] <= 46.27 &  # Latitude range
      Fathom >= 75 & Fathom <= 100
  )

#this is the central california area 
fathom_area_2 <- fathoms_sf %>%
  filter(
    st_coordinates(.)[,2] >= 34.45 & st_coordinates(.)[,2] <= 40.17 &  # Latitude range
      Fathom >= 75 & Fathom <= 125
  )

# Combine the two areas into one filtered dataset
fathom_subset <- bind_rows(fathom_area_1, fathom_area_2)

# Filter using st_coordinates to extract Latitude and Longitude for filtering
#oregon/northern california area 
fathom_area_high <- fathom_subset %>%
  filter(
    st_coordinates(.)[, 2] >= 40.17 &  # Latitude
      st_coordinates(.)[, 2] <= 46.5 &
      Fathom >= 75 & Fathom <= 100
  )

#central california area 
fathom_area_low <- fathom_subset %>%
  filter(
    st_coordinates(.)[, 2] >= 34.45 &
      st_coordinates(.)[, 2] < 40.17 &
      Fathom >= 75 & Fathom <= 125
  )

plot_habitat_with_fathoms <- function(habitat_df, fathom_sf, month, area) {
  # Ensure CRS is set for both habitat_df and fathom_sf
  if (is.na(st_crs(fathom_sf))) {
    st_crs(fathom_sf) <- 4326
  }
  
  if (!inherits(habitat_df, "sf")) {
    habitat_sf <- st_as_sf(habitat_df, coords = c("x", "y"), crs = 4326)
  } else {
    habitat_sf <- habitat_df
  }
  
  # Filter fathom_sf based on the area and desired depths
  if (area == "high") {
    fathom_sf <- fathom_sf %>% filter(Fathom %in% c(75, 100))
  } else {
    fathom_sf <- fathom_sf %>% filter(Fathom %in% c(75, 125))
  }
  
  # Extract coordinates from the sf object to use in geom_tile
  habitat_sf <- habitat_sf %>%
    mutate(x = st_coordinates(.)[, 1], y = st_coordinates(.)[, 2])
  
  ggplot() +
    geom_tile(data = habitat_sf, aes(x = x, y = y, fill = mean_suitability)) +
    scale_fill_viridis_c() +
    geom_sf(data = fathom_sf, aes(color = factor(Fathom)), fill = NA, size = 0.2, linetype = "dashed") +  # Color and label fathom lines
    scale_color_manual(values = c("orange", "red", "pink")) +  # Custom colors for the fathom lines
    theme_minimal() +
    coord_sf(
      xlim = if (area == "high") {
        c(-125.5, -123)  # X limits for high area
      } else {
        c(-125, -120)  # X limits for low area
      },
      ylim = if (area == "high") {
        c(40.17, 46.5)  # High area latitude limits
      } else {
        c(34.5, 40.17)  # Low area latitude limits
      },
      crs = st_crs(habitat_sf)
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(title = paste("Average Humpback Whale Density  -", month),
         x = "Longitude", y = "Latitude", fill = "Whale Density", color = "Fathom Depth (m)") +  # Label the fathom lines
    scale_x_continuous(breaks = if (area == "high") {
      seq(-125.5, -123, by = 0.5)  # Breaks for high area
    } else {
      seq(-125, -120, by = 1)  # Breaks for low area
    }) +
    scale_y_continuous(
      breaks = if (area == "high") {
        seq(40.17, 46.5, by = 1)
      } else {
        seq(34.5, 40.17, by = 1)
      }
    )
}

# Create plots for each month and each area
#northern california/oregon area
plots_high_hbw <- map2(monthly_avg_habitat, names(monthly_avg_habitat), ~ plot_habitat_with_fathoms(.x, fathom_area_high, .y, "high"))
#central california area 
plots_low_hbw <- map2(monthly_avg_habitat, names(monthly_avg_habitat), ~ plot_habitat_with_fathoms(.x, fathom_area_low, .y, "low"))

#then to pull the respective months you can just select from the plots_high or plots_low whatever month you want
#usinga number identified for month in a double bracket 

#LEATHERBACKS 
# Compute the average humpback whale density for each month across years (2014-2023)
monthly_avg_habitat <- map(months, ~ compute_monthly_average(.x, lbstmonthdf))

# Name the list elements
names(monthly_avg_habitat) <- months

#this is a dataframe that I created that has fathoms associated with coordinate values along the west coast 
fathoms<-read.csv('fathom.csv')

# Convert the fathoms dataframe to an sf object
fathoms_sf <- st_as_sf(fathoms, coords = c("Longitude", "Latitude"), crs = 4326)

# Filter using the geometry column to check if coordinates fall within the specified ranges
#this is the northern california/oregon area 
fathom_area_1 <- fathoms_sf %>%
  filter(
    st_coordinates(.)[,2] >= 40.17 & st_coordinates(.)[,2] <= 46.27 &  # Latitude range
      Fathom >= 75 & Fathom <= 100
  )

#this is the central california area 
fathom_area_2 <- fathoms_sf %>%
  filter(
    st_coordinates(.)[,2] >= 34.45 & st_coordinates(.)[,2] <= 40.17 &  # Latitude range
      Fathom >= 75 & Fathom <= 125
  )

# Combine the two areas into one filtered dataset
fathom_subset <- bind_rows(fathom_area_1, fathom_area_2)

# Filter using st_coordinates to extract Latitude and Longitude for filtering
#oregon/northern california area 
fathom_area_high <- fathom_subset %>%
  filter(
    st_coordinates(.)[, 2] >= 40.17 &  # Latitude
      st_coordinates(.)[, 2] <= 46.5 &
      Fathom >= 75 & Fathom <= 100
  )

#central california area 
fathom_area_low <- fathom_subset %>%
  filter(
    st_coordinates(.)[, 2] >= 34.45 &
      st_coordinates(.)[, 2] < 40.17 &
      Fathom >= 75 & Fathom <= 125
  )

plot_habitat_with_fathoms <- function(habitat_df, fathom_sf, month, area) {
  # Ensure CRS is set for both habitat_df and fathom_sf
  if (is.na(st_crs(fathom_sf))) {
    st_crs(fathom_sf) <- 4326
  }
  
  if (!inherits(habitat_df, "sf")) {
    habitat_sf <- st_as_sf(habitat_df, coords = c("x", "y"), crs = 4326)
  } else {
    habitat_sf <- habitat_df
  }
  
  # Filter fathom_sf based on the area and desired depths
  if (area == "high") {
    fathom_sf <- fathom_sf %>% filter(Fathom %in% c(75, 100))
  } else {
    fathom_sf <- fathom_sf %>% filter(Fathom %in% c(75, 125))
  }
  
  # Extract coordinates from the sf object to use in geom_tile
  habitat_sf <- habitat_sf %>%
    mutate(x = st_coordinates(.)[, 1], y = st_coordinates(.)[, 2])
  
  ggplot() +
    geom_tile(data = habitat_sf, aes(x = x, y = y, fill = mean_suitability)) +
    scale_fill_viridis_c() +
    geom_sf(data = fathom_sf, aes(color = factor(Fathom)), fill = NA, size = 0.2, linetype = "dashed") +  # Color and label fathom lines
    scale_color_manual(values = c("orange", "red", "pink")) +  # Custom colors for the fathom lines
    theme_minimal() +
    coord_sf(
      xlim = if (area == "high") {
        c(-125.5, -123)  # X limits for high area
      } else {
        c(-125, -120)  # X limits for low area
      },
      ylim = if (area == "high") {
        c(40.17, 46.5)  # High area latitude limits
      } else {
        c(34.5, 40.17)  # Low area latitude limits
      },
      crs = st_crs(habitat_sf)
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(title = paste("Average Leatherback Habitat Suitability  -", month),
         x = "Longitude", y = "Latitude", fill = "Suitability", color = "Fathom Depth (m)") +  # Label the fathom lines
    scale_x_continuous(breaks = if (area == "high") {
      seq(-125.5, -123, by = 0.5)  # Breaks for high area
    } else {
      seq(-125, -120, by = 1)  # Breaks for low area
    }) +
    scale_y_continuous(
      breaks = if (area == "high") {
        seq(40.17, 46.5, by = 1)
      } else {
        seq(34.5, 40.17, by = 1)
      }
    )
}

# Create plots for each month and each area
#northern california/oregon area
plots_high_lbst <- map2(monthly_avg_habitat, names(monthly_avg_habitat), ~ plot_habitat_with_fathoms(.x, fathom_area_high, .y, "high"))
#central california area 
plots_low_lbst <- map2(monthly_avg_habitat, names(monthly_avg_habitat), ~ plot_habitat_with_fathoms(.x, fathom_area_low, .y, "low"))

#then to pull the respective months you can just select from the plots_high or plots_low whatever month you want
#usinga number identified for month in a double bracket 