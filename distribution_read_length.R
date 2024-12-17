##############################################################################
################### READ LENGTH ANALYSIS JOSEFINE WEIß #######################
####################### IF YOU USE, PLEASE CITE ##############################
##############################################################################

# First of all you have to read in following documents for KL12

# biodiv/user/joweiss/KL12/INFO_samples_kl12.txt
# biodiv/user/joweiss/KL12/output/out.kraken2/Age_KL12.txt

#KL12
####### LOOP FOR READ COUNTS ################
library(data.table)

folder_path_info <- "/Volumes/projects/biodiv/user/joweiss/KL12/out.read_extract_Embryophyta"

# Get a list of all files in the folder
file_list_info_1 <- list.files(path = folder_path_info, full.names = TRUE)

# Filter the file list to include only files with "_merged_tax3193_extract_length" in their names
file_list_info_2 <- file_list_info_1[grep("_tax3193_extract_merged_count.txt", file_list_info_1)]

# Create an empty list to store the data tables for each file
table_list <- list()

# Loop through each file in the folder
table_list <- lapply(file_list_info_2, function(file_path) {
  # Read the data from the file (assuming it's a single value per file)
  value <- as.numeric(readLines(file_path))
  
  # Extract the file name without the path and extension
  file_name_info <- basename(file_path)
  
  # Create a data table with the file name and value
  data.table(file_name = file_name_info, value = value)
})

# Combine all the data tables into a single data table
result_table <- rbindlist(table_list, fill = TRUE)

# Print the resulting table
print(result_table)

# Write the table to a CSV file
output_file <- "/Volumes/projects/biodiv/user/joweiss/KL12/out.read_extract_Embryophyta/result_table_KL12.csv"
write.csv(result_table, file = output_file, row.names = FALSE)

#######LOOP FOR READ LENGTH ################
library(ggplot2)
library(tools)  
library(stringr)# Required for file_path_sans_ext function
library(gridExtra)



theme_set(theme(plot.title = element_text(size = 6)))



# Specify the folder path where the files are located
folder_path <- "/Volumes/projects/biodiv/user/joweiss/KL12/output/out.kraken2/out.length_extract_Embryophyta"

output_path <- "/Volumes/projects/biodiv/user/joweiss/KL12/output/out.kraken2/PLOTS"
# Get a list of all files in the specified folder
file_list <- list.files(path = folder_path, full.names = TRUE)

# Filter the file list to include only files with "_merged_tax3193_extract_length" in their names
file_list_data <- file_list[grep("_merged_tax3193_extract_length", file_list)]

# Get a list of all files with "_merged_tax3193_extract_length" in the specified folder
#file_list <- list.files(path ="/Volumes/projects/biodiv/user/joweiss/KL12/output/out.kraken2/out.length_extract_Embryophyta", pattern = "_merged_tax3193_extract_length", full.names = TRUE)

# set NA to 0 in Age_KL12
Age_KL12$age[is.na(Age_KL12$age)] <- 0


# Set up an empty list to store the plots
plot_list <- list()

# MEDIAN Calculations 
medians <- list()

# Create an empty table to store the data
result_data_for_plot_KL12 <- data.frame(age_value = numeric(),
                          V2 = numeric(),
                          file_name_partial = character(),
                          stringsAsFactors = FALSE)

# Loop through each file in the list
for (file in file_list_data) {
  # Check if the file is empty or has no content
  if (is.na(file.size(file)) || file.size(file) == 0) {
    # Skip this file and move on to the next iteration
    next
  }
  # Read the data from the file
  data <- read.table(file, header=FALSE)
  
  
  #extract the file name
  file_name_data <- file_path_sans_ext(basename(file))
 
   # Extract the first 29 characters of the file name
  file_name_partial <- substr(file_name_data, 1, 30)

  # Find the matching row in the results table based on the partial file name
  matching_rows <- result_table[str_detect(result_table$file_name, file_name_partial), ]
  
  # Check if there are matching rows
  if (nrow(matching_rows) > 0) {
    # Extract the value from the first matching row
    info_value <- matching_rows$value[1]
  } else {
    # Set info_value to a default value or handle the case when there are no matches
    info_value <- "N/A"
  }
  
  # Find the matching row in the result table for Age_KL12 file
  age_matching_row <- Age_KL12[str_detect(Age_KL12$samples, file_name_partial), ]
  

  # Extract the Age value from the matching row
  age_value <- age_matching_row$age
  
  # Check if there are matching rows
#  if (nrow(age_matching_row) > 0) {
    # Extract the value from the first matching row
 #   age_value <- age_matching_row$age[1]
#  } else {
    # Set info_value to a default value or handle the case when there are no matches
 #   age_value <- "N/A"
  #}
  
  # Add the data to the result_data table
  result_data_for_plot_KL12 <- rbind(result_data_for_plot_KL12, data.frame(age_value = age_value,
                                               V2 = data$V2,
                                               file_name_partial = file_name_partial,
                                               stringsAsFactors = FALSE))  
  # Create the ggplot objects
  plot1 <- ggplot(data, aes(x = V2)) +
    geom_histogram(binwidth = 0.5, colour = "black", fill = "white") +
    labs(title= paste("Age:", age_value),x= info_value)
  
  plot2 <- ggplot(data, aes(x = V2)) +
    geom_density() +
    labs(title= paste("Age:", age_value), x= info_value)
  
  # Generate a dynamic name for the plot
  plot_name <- paste0("X", gsub("-", "_", file_name_partial), "_plot3")
  
  # Create the plot with the dynamic name
  assign(plot_name, ggplot(data, aes(x = V2)) +
           geom_histogram(aes(y = after_stat(density)),
                          binwidth = 0.5,
                          colour = "black", fill = "white") +
           geom_density(alpha = 0.2, fill = "#FF6666") +
           labs(title = paste("Age:", age_value, gsub("-", "_", file_name_partial)), x = info_value))
  
  print(plot_name)
  
}


# Write the table to a CSV file
output_file <- "/Volumes/projects/biodiv/user/joweiss/KL12/result_data_for_plot_KL12.csv"
write.csv(result_data_for_plot_KL12, file = output_file, row.names = FALSE)



library(dplyr)

# Calculate the average read length over time
average_read_length_KL12 <- result_data_for_plot_KL12 %>%
  group_by(age_value) %>%
  summarize(average_length = mean(V2))

# Print the average read length over time
print(average_read_length_KL12)

# Remove the first three rows from the data table
average_read_length_KL12 <- average_read_length_KL12[-c(1:3), ]

# Print the updated data table
print(average_read_length_KL12)


# Create a geom_segment plot of read length over time (age_value)
plot_segments <- ggplot(average_read_length_KL12, aes(x = age_value, y = average_length, xend = age_value, yend = 0)) +
  geom_col_segs(color = "black") +
  scale_x_continuous(breaks = seq(0, 20000, by = 2000)) +
  labs(title = "Read Length over Time (KL12)", x = "Age (yrs)", y = "Mean Read Length per sample") +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(), axis.ticks = element_line())
print(plot_segments)



# Plot boxplots over continuous time
ggplot(result_data_for_plot_KL12, aes(y = V2, group= file_name_partial)) +
  geom_boxplot() +
  labs(x = "Continuous Time", y = "Measure") +
  theme_minimal()


# Convert age_value to numeric
result_data_for_plot_KL12$age_value <- as.numeric(as.character(result_data_for_plot_KL12$age_value))

# Exclude rows with age_value of 1, 2, 3, or 4
filtered_data <- Result_data_for_plot_KL12 %>%
  filter(!(age_value %in% c(1, 2, 3, 4)))


filtered_data <- Result_data_for_plot_Australia %>%
  filter(age_value != 0)

filtered_data <- Result_data_for_Plot_TobagoForReal %>%
  filter(age_value != 0)

# Plot boxplots over time with filtered data
ggplot(filtered_data, aes(x = factor(age_value), y = V2)) +
  geom_boxplot() +
  labs(x = "Age (yrs)", y = "Read Length") +
  scale_x_discrete() +
  theme_minimal()




plot_grid(X200313_A00902_A_L002_APMG_9_11_plot3, 
          X200313_A00902_A_L002_APMG_9_10_plot3,
          X200313_A00902_A_L002_APMG_9_9__plot3,
          X200313_A00902_A_L002_APMG_9_8__plot3, 
          X200313_A00902_A_L002_APMG_9_7__plot3, 
          X200313_A00902_A_L002_APMG_9_6__plot3,
          X191011_SND405_A_L007_APMG_6_6__plot3, 
          X200313_A00902_A_L002_APMG_9_5__plot3, 
          X200313_A00902_A_L002_APMG_9_3__plot3, 
          X200313_A00902_A_L002_APMG_9_2__plot3, 
          X200313_A00902_A_L002_APMG_9_1__plot3, 
          X200313_A00902_A_L001_APMG_8_12_plot3,
          X191011_SND405_A_L007_APMG_6_5__plot3, 
          X200313_A00902_A_L001_APMG_8_8__plot3,
          X200313_A00902_A_L001_APMG_8_11_plot3, 
          X200313_A00902_A_L001_APMG_8_10_plot3,
          X200313_A00902_A_L001_APMG_8_9__plot3,
          X191011_SND405_A_L007_APMG_6_4__plot3, 
          X200313_A00902_A_L001_APMG_8_6__plot3, 
          X191011_SND405_A_L007_APMG_6_3__plot3, 
          X200313_A00902_A_L001_APMG_8_5__plot3,
          X200313_A00902_A_L001_APMG_8_4__plot3,
          X191011_SND405_A_L007_APMG_6_2__plot3,
          X200313_A00902_A_L001_APMG_8_3__plot3,
          X200313_A00902_A_L001_APMG_8_2__plot3, 
          X200313_A00902_A_L001_APMG_8_1__plot3, 
          X191011_SND405_A_L007_APMG_6_1__plot3, 
          X191011_SND405_A_L007_APMG_6_8__plot3, 
          X200313_A00902_A_L001_APMG_8_7__plot3, 
          X200313_A00902_A_L002_APMG_9_12_plot3)
  


#Tobago

# First of all you have to read in following documents for TOBAGO

# biodiv/user/joweiss/Tobago/Age_Tobago.txt
# biodiv/user/joweiss/Tobago/Tobago_info.txt

#Tobago





#TOBAGo









####### LOOP FOR READ COUNTS ################
library(data.table)

folder_path_info <- "/Volumes/projects/biodiv/user/joweiss/Tobago/out.read_extract_Embryophyta"

# Get a list of all files in the folder
file_list_info_1 <- list.files(path = folder_path_info, full.names = TRUE)

# Filter the file list to include only files with "_merged_tax3193_extract_length" in their names
file_list_info_2 <- file_list_info_1[grep("_tax3193_extract_merged_count.txt", file_list_info_1)]

# Create an empty list to store the data tables for each file
table_list <- list()

# Loop through each file in the folder
table_list <- lapply(file_list_info_2, function(file_path) {
  # Read the data from the file (assuming it's a single value per file)
  value <- as.numeric(readLines(file_path))
  
  # Extract the file name without the path and extension
  file_name_info <- basename(file_path)
  
  # Create a data table with the file name and value
  data.table(file_name = file_name_info, value = value)
})

# Combine all the data tables into a single data table
result_table <- rbindlist(table_list, fill = TRUE)

# Print the resulting table
print(result_table)

# Write the table to a CSV file
output_file <- "/Volumes/projects/biodiv/user/joweiss/Tobago/out.read_extract_Embryophyta/result_table_Tobago.csv"
write.csv(result_table, file = output_file, row.names = FALSE)

#######LOOP FOR READ LENGTH ################
library(ggplot2)
library(tools)  
library(stringr)# Required for file_path_sans_ext function
library(gridExtra)



theme_set(theme(plot.title = element_text(size = 6)))

# Specify the folder path where the files are located
folder_path <- "/Volumes/projects/biodiv/user/joweiss/Tobago/output/out.kraken2/out.length_extract_Embryophyta"

output_path <- "/Volumes/projects/biodiv/user/joweiss/Tobago/output/Plots"
# Get a list of all files in the specified folder
file_list <- list.files(path = folder_path, full.names = TRUE)

# Filter the file list to include only files with "_merged_tax3193_extract_length" in their names
file_list_data <- file_list[grep("_merged_tax3193_extract_length", file_list)]

# Get a list of all files with "_merged_tax3193_extract_length" in the specified folder
#file_list <- list.files(path ="/Volumes/projects/biodiv/user/joweiss/KL12/output/out.kraken2/out.length_extract_Embryophyta", pattern = "_merged_tax3193_extract_length", full.names = TRUE)

# set NA to 0 in Age_KL12
Age_Tobago$age[is.na(Age_Tobago$age)] <- 0


# Set up an empty list to store the plots
plot_list <- list()

# empty data frame for later
result_data_for_plot_Tobago <- data.frame(age_value = numeric(),
                                          V2 = numeric(),
                                          file_name_partial = character(),
                                          stringsAsFactors = FALSE)# empty data frame for later

# Loop through each file in the list
for (file in file_list_data) {
  # Check if the file is empty or has no content
  if (is.na(file.size(file)) || file.size(file) == 0) {
    # Skip this file and move on to the next iteration
    next
  }

  # Read the data from the file
  data <- read.table(file, header=FALSE)
  
  #extract the file name
  file_name_data <- file_path_sans_ext(basename(file))
  
  # Extract the first 29 characters of the file name
  file_name_partial <- substr(file_name_data, 1, 13)
  
  # Find the matching row in the results table based on the partial file name
  matching_rows <- result_table[grepl(file_name_partial, result_table$file_name), ]
  
  # Check if there are matching rows
  if (nrow(matching_rows) > 0) {
    # Extract the value from the first matching row
    info_value <- matching_rows$value[1]
  } else {
    # Set info_value to a default value or handle the case when there are no matches
    info_value <- "N/A"
  }
  
  # Find the matching row in the result table for Age_KL12 file
  age_matching_row <- Age_Tobago[str_detect(Age_Tobago$samples, file_name_partial), ]
  
  # Check if age_matching_row is found
  if (nrow(age_matching_row) == 0) {
    # Skip this file and move on to the next iteration
    next
  }
  
  
  # Extract the Age value from the matching row
  age_value <- age_matching_row$age
  
  # Check if there are matching rows
  #  if (nrow(age_matching_row) > 0) {
  # Extract the value from the first matching row
  #   age_value <- age_matching_row$age[1]
  #  } else {
  # Set info_value to a default value or handle the case when there are no matches
  #   age_value <- "N/A"
  #}
  # Add the data to the result_data table
  result_data_for_plot_Tobago <- rbind(result_data_for_plot_Tobago, data.frame(age_value = age_value,
                                                                           V2 = data$V2,
                                                                           file_name_partial = file_name_partial,
                                                                           stringsAsFactors = FALSE))  
  
  # Create the ggplot objects
  plot1 <- ggplot(data, aes(x = V2)) +
    geom_histogram(binwidth = 0.5, colour = "black", fill = "white") +
    labs(title= paste("Age:", age_value),x= info_value)
  
  plot2 <- ggplot(data, aes(x = V2)) +
    geom_density() +
    labs(title= paste("Age:", age_value), x= info_value)
  
  # Generate a dynamic name for the plot
  plot_name <- paste0("X", gsub("-", "_", file_name_partial), "_plot3")
  
  # Create the plot with the dynamic name
  assign(plot_name, ggplot(data, aes(x = V2)) +
           geom_histogram(aes(y = after_stat(density)),
                          binwidth = 0.5,
                          colour = "black", fill = "white") +
           geom_density(alpha = 0.2, fill = "#FF6666") +
           labs(title = paste("Age:", age_value, gsub("-", "_", file_name_partial)), x = info_value))
  
  print(plot_name)
  
  # Print or save the plots here (e.g., using ggsave or print functions)
  #print(plot3)
  
  # Add the plot to the plot list
  
}


# Calculate the average read length over time
average_read_length_Tobago <- result_data_for_plot_Tobago %>%
  group_by(age_value) %>%
  summarize(average_length = mean(V2))

# Print the average read length over time
print(average_read_length_Tobago)

# Remove the first three rows from the data table
average_read_length_Tobago <- average_read_length_Tobago[-c(1:3), ]

# Print the updated data table
print(average_read_length_Tobago)

output_file <- "/Volumes/projects/biodiv/user/joweiss/Tobago/result_data_for_plot_Tobago_new.csv"
write.csv(result_data_for_plot_Tobago, file = output_file, row.names = FALSE)

# Create a geom_segment plot of read length over time (age_value)
plot_segments <- ggplot(average_read_length_KL12, aes(x = age_value, y = average_length, xend = age_value, yend = 0)) +
  geom_col_segs(color = "black") +
  scale_x_continuous(breaks = seq(0, 20000, by = 2000)) +
  labs(title = "Read Length over Time (KL12)", x = "Age (yrs)", y = "Mean Read Length per sample") +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(), axis.ticks = element_line())
print(plot_segments)

plot_grid(XJK0134L_1_S14_plot3,
          XJK0134L_11_S2_plot3,
          XJK0134L_12_S2_plot3,
          XJK0134L_13_S2_plot3,
          XJK0134L_14_S2_plot3,
          XJK0134L_15_S2_plot3,
          XJK0134L_2_S15_plot3,
          XJK0134L_3_S16_plot3,
          XJK0134L_4_S17_plot3,
          XJK0134L_5_S18_plot3,
          XJK0134L_6_S19_plot3,
          XJK0134L_7_S20_plot3,
          XJK0134L_8_S21_plot3,
          XJK0134L_9_S22_plot3,
          XJK0134L_10_S2_plot3,
          XJK0134L_16_S2_plot3,
          XJK0135L_1_S30_plot3,
          XJK0135L_2_S31_plot3,
          XJK0135L_3_S32_plot3,
          XJK0135L_4_S33_plot3,
          XJK0135L_5_bla_plot3,
          XJK0135L_6_ntc_plot3,
          XJK0136L_12_S1_plot3,
          XJK0136L_11_S1_plot3,
          XJK0136L_10_S1_plot3,
          XJK0136L_13_S1_plot3,
          XJK0136L_9_S9__plot3,
          XJoW007L_1_S1__plot3,
          XJoW007L_2_S2__plot3,
          XJoW007L_3_S3__plot3,
          XJoW007L_4_S4__plot3,
          XJoW007L_5_S5__plot3,
          XJoW007L_6_S6__plot3,
          XJoW007L_7_S7__plot3)



# First of all you have to read in following documents for MSM

# biodiv/user/joweiss/MSM/Age_MSM.txt
# biodiv/user/joweiss/MSM/INFO_Age_MSM.txt

#MSM
















####### LOOP FOR READ COUNTS ################
library(data.table)

folder_path_info <- "/Volumes/projects/biodiv/user/joweiss/MSM/out.read_extract_Embryophyta"

# Get a list of all files in the folder
file_list_info_1 <- list.files(path = folder_path_info, full.names = TRUE)

# Filter the file list to include only files with "_merged_tax3193_extract_length" in their names
file_list_info_2 <- file_list_info_1[grep("_tax3193_extract_merged_count.txt", file_list_info_1)]

# Create an empty list to store the data tables for each file
table_list <- list()

# Loop through each file in the folder
table_list <- lapply(file_list_info_2, function(file_path) {
  # Read the data from the file (assuming it's a single value per file)
  value <- as.numeric(readLines(file_path))
  
  # Extract the file name without the path and extension
  file_name_info <- basename(file_path)
  
  # Create a data table with the file name and value
  data.table(file_name = file_name_info, value = value)
})

# Combine all the data tables into a single data table
result_table <- rbindlist(table_list, fill = TRUE)

# Print the resulting table
print(result_table)

# Write the table to a CSV file
output_file <- "/Volumes/projects/biodiv/user/joweiss/MSM/out.read_extract_Embryophyta/result_table_MSM.csv"
write.csv(result_table, file = output_file, row.names = FALSE)

#######LOOP FOR READ LENGTH ################
library(ggplot2)
library(tools)  
library(stringr)# Required for file_path_sans_ext function
library(gridExtra)



theme_set(theme(plot.title = element_text(size = 6)))



# Specify the folder path where the files are located
folder_path <- "/Volumes/projects/biodiv/user/joweiss/MSM/out.length_extract_Embryophyta"

output_path <- "/Volumes/projects/biodiv/user/joweiss/MSM"
# Get a list of all files in the specified folder
file_list <- list.files(path = folder_path, full.names = TRUE)

# Filter the file list to include only files with "_merged_tax3193_extract_length" in their names
file_list_data <- file_list[grep("_merged_tax3193_extract_length", file_list)]

# Get a list of all files with "_merged_tax3193_extract_length" in the specified folder
#file_list <- list.files(path ="/Volumes/projects/biodiv/user/joweiss/KL12/output/out.kraken2/out.length_extract_Embryophyta", pattern = "_merged_tax3193_extract_length", full.names = TRUE)

# set NA to 0 in Age_KL12
Age_MSM$age[is.na(Age_MSM$age)] <- 0


# Set up an empty list to store the plots
plot_list <- list()

# MEDIAN Calculations 
medians <- list()

# Create an empty table to store the data
result_data_for_plot_MSM <- data.frame(age_value = numeric(),
                                        V2 = numeric(),
                                        file_name_partial = character(),
                                        stringsAsFactors = FALSE)


# Loop through each file in the list
for (file in file_list_data) {
  # Check if the file is empty or has no content
  if (is.na(file.size(file)) || file.size(file) == 0) {
    # Skip this file and move on to the next iteration
    next
  }
  # Read the data from the file
  data <- read.table(file, header=FALSE)
  
  
  #extract the file name
  file_name_data <- file_path_sans_ext(basename(file))
  
  # Extract the first 29 characters of the file name
  file_name_partial <- substr(file_name_data, 1, 12)
  
  # Find the matching row in the results table based on the partial file name
  matching_rows <- result_table[str_detect(result_table$file_name, file_name_partial), ]
  
  # Check if there are matching rows
  if (nrow(matching_rows) > 0) {
    # Extract the value from the first matching row
    info_value <- matching_rows$value[1]
  } else {
    # Set info_value to a default value or handle the case when there are no matches
    info_value <- "N/A"
  }
  
  # Find the matching row in the result table for Age_KL12 file
  age_matching_row <- Age_MSM[str_detect(Age_MSM$file_name, file_name_partial), ]
  
  
  # Extract the Age value from the matching row
  age_value <- age_matching_row$age
  
  # Check if there are matching rows
  #  if (nrow(age_matching_row) > 0) {
  # Extract the value from the first matching row
  #   age_value <- age_matching_row$age[1]
  #  } else {
  # Set info_value to a default value or handle the case when there are no matches
  #   age_value <- "N/A"
  #}
  
  result_data_for_plot_MSM <- rbind(result_data_for_plot_MSM, data.frame(age_value = age_value,
                                                                               V2 = data$V2,
                                                                               file_name_partial = file_name_partial,
                                                                               stringsAsFactors = FALSE))  
  
  # Create the ggplot objects
  plot1 <- ggplot(data, aes(x = V2)) +
    geom_histogram(binwidth = 0.5, colour = "black", fill = "white") +
    labs(title= paste("Age:", age_value),x= info_value)
  
  plot2 <- ggplot(data, aes(x = V2)) +
    geom_density() +
    labs(title= paste("Age:", age_value), x= info_value)
  
  # Generate a dynamic name for the plot
  plot_name <- paste0("X", gsub("-", "_", file_name_partial), "_plot3")
  
  # Create the plot with the dynamic name
  assign(plot_name, ggplot(data, aes(x = V2)) +
           geom_histogram(aes(y = after_stat(density)),
                          binwidth = 0.5,
                          colour = "black", fill = "white") +
           geom_density(alpha = 0.2, fill = "#FF6666") +
           labs(title = paste("Age:", age_value, gsub("-", "_", file_name_partial)), x = info_value))
  
  print(plot_name)
  
}

# Write the table to a CSV file
output_file <- "/Volumes/projects/biodiv/user/joweiss/MSM/result_data_for_plot_MSMpart2.csv"
write.csv(result_data_for_plot_MSM, file = output_file, row.names = FALSE)

# Calculate the average read length over time
#average_read_length_Tobago <- result_data_for_plot_Tobago %>%
 # group_by(age_value) %>%
  #summarize(average_length = mean(V2))

# Print the average read length over time
#print(average_read_length_Tobago)

# Remove the first three rows from the data table
#average_read_length_Tobago <- average_read_length_Tobago[-c(1:3), ]

# Print the updated data table
#print(average_read_length_Tobago)


# Convert age_value to numeric
result_data_for_plot_MSMgesamt$age_value <- as.numeric(as.character(result_data_for_plot_MSMgesamt$age_value))

filtered_data <- result_data_for_plot_MSMgesamt %>%
  filter(age_value != 0)


# Plot boxplots over time with filtered data
ggplot(filtered_data, aes(x = factor(age_value), y = V2)) +
  geom_boxplot() +
  labs(x = "Age (yrs)", y = "Read Length") +
  scale_x_discrete() +
  theme_minimal()







plot_grid(XJK0137L_1_S1_plot3, XJK0137L_2_S2_plot3, XJK0137L_3_S3_plot3, XJK0137L_4_S4_plot3, XJK0137L_5_S5_plot3, XJK0137L_6_S6_plot3, XJK0137L_7_S7_plot3, XJK0137L_8_S8_plot3)

# First of all you have to read in following documents for PS97

# biodiv/user/joweiss/PS97/Age_PS97.txt
# biodiv/user/joweiss/PS97/INFO_age_PS97.txt

#PS97








####### LOOP FOR READ COUNTS ################
library(data.table)

folder_path_info <- "/Volumes/projects/biodiv/user/joweiss/PS97/out.read_extract_Embryophyta"

# Get a list of all files in the folder
file_list_info_1 <- list.files(path = folder_path_info, full.names = TRUE)

# Filter the file list to include only files with "_merged_tax3193_extract_length" in their names
file_list_info_2 <- file_list_info_1[grep("_tax3193_extract_merged_count.txt", file_list_info_1)]

# Create an empty list to store the data tables for each file
table_list <- list()

# Loop through each file in the folder
table_list <- lapply(file_list_info_2, function(file_path) {
  # Read the data from the file (assuming it's a single value per file)
  value <- as.numeric(readLines(file_path))
  
  # Extract the file name without the path and extension
  file_name_info <- basename(file_path)
  
  # Create a data table with the file name and value
  data.table(file_name = file_name_info, value = value)
})

# Combine all the data tables into a single data table
result_table <- rbindlist(table_list, fill = TRUE)

# Print the resulting table
print(result_table)

# Write the table to a CSV file
output_file <- "/Volumes/projects/biodiv/user/joweiss/PS97/out.read_extract_Embryophyta/result_table_MSM.csv"
write.csv(result_table, file = output_file, row.names = FALSE)

#######LOOP FOR READ LENGTH ################
library(ggplot2)
library(tools)  
library(stringr)# Required for file_path_sans_ext function
library(gridExtra)



theme_set(theme(plot.title = element_text(size = 6)))



# Specify the folder path where the files are located
folder_path <- "/Volumes/projects/biodiv/user/joweiss/PS97/out.length_extract_Embryophyta"

output_path <- "/Volumes/projects/biodiv/user/joweiss/PS97"
# Get a list of all files in the specified folder
file_list <- list.files(path = folder_path, full.names = TRUE)

# Filter the file list to include only files with "_merged_tax3193_extract_length" in their names
file_list_data <- file_list[grep("_merged_tax3193_extract_length", file_list)]

# Get a list of all files with "_merged_tax3193_extract_length" in the specified folder
#file_list <- list.files(path ="/Volumes/projects/biodiv/user/joweiss/KL12/output/out.kraken2/out.length_extract_Embryophyta", pattern = "_merged_tax3193_extract_length", full.names = TRUE)

# set NA to 0 in Age_KL12
Age_PS97$age[is.na(Age_PS97$age)] <- 0


# Set up an empty list to store the plots
plot_list <- list()

# MEDIAN Calculations 
medians <- list()

# Create an empty table to store the data
result_data_for_plot_PS97 <- data.frame(age_value = numeric(),
                                       V2 = numeric(),
                                       file_name_partial = character(),
                                       stringsAsFactors = FALSE)


# Loop through each file in the list
for (file in file_list_data) {
  # Check if the file is empty or has no content
  if (is.na(file.size(file)) || file.size(file) == 0) {
    # Skip this file and move on to the next iteration
    next
  }
  # Read the data from the file
  data <- read.table(file, header=FALSE)
  
  
  #extract the file name
  file_name_data <- file_path_sans_ext(basename(file))
  
  # Extract the first 29 characters of the file name
  file_name_partial <- substr(file_name_data, 1, 12)
  
  # Find the matching row in the results table based on the partial file name
  matching_rows <- result_table[str_detect(result_table$file_name, file_name_partial), ]
  
  # Check if there are matching rows
  if (nrow(matching_rows) > 0) {
    # Extract the value from the first matching row
    info_value <- matching_rows$value[1]
  } else {
    # Set info_value to a default value or handle the case when there are no matches
    info_value <- "N/A"
  }
  
  # Find the matching row in the result table for Age_KL12 file
  age_matching_row <- Age_PS97[str_detect(Age_PS97$file_name, file_name_partial), ]
  
  
  # Extract the Age value from the matching row
  age_value <- age_matching_row$age
  
  # Check if there are matching rows
  #  if (nrow(age_matching_row) > 0) {
  # Extract the value from the first matching row
  #   age_value <- age_matching_row$age[1]
  #  } else {
  # Set info_value to a default value or handle the case when there are no matches
  #   age_value <- "N/A"
  #}
  
  result_data_for_plot_PS97 <- rbind(result_data_for_plot_PS97, data.frame(age_value = age_value,
                                                                         V2 = data$V2,
                                                                         file_name_partial = file_name_partial,
                                                                         stringsAsFactors = FALSE))  
  
  
  
  # Create the ggplot objects
  plot1 <- ggplot(data, aes(x = V2)) +
    geom_histogram(binwidth = 0.5, colour = "black", fill = "white") +
    labs(title= paste("Age:", age_value),x= info_value)
  
  plot2 <- ggplot(data, aes(x = V2)) +
    geom_density() +
    labs(title= paste("Age:", age_value), x= info_value)
  
  # Generate a dynamic name for the plot
  plot_name <- paste0("X", gsub("-", "_", file_name_partial), "_plot3")
  
  # Create the plot with the dynamic name
  assign(plot_name, ggplot(data, aes(x = V2)) +
           geom_histogram(aes(y = after_stat(density)),
                          binwidth = 0.5,
                          colour = "black", fill = "white") +
           geom_density(alpha = 0.2, fill = "#FF6666") +
           labs(title = paste("Age:", age_value, gsub("-", "_", file_name_partial)), x = info_value))
  
  print(plot_name)
  
}

# Write the table to a CSV file
output_file <- "/Volumes/projects/biodiv/user/joweiss/PS97/result_data_for_plot_PS97.csv"
write.csv(result_data_for_plot_PS97, file = output_file, row.names = FALSE)

# Calculate the average read length over time
#average_read_length_Tobago <- result_data_for_plot_Tobago %>%
# group_by(age_value) %>%
#summarize(average_length = mean(V2))

# Print the average read length over time
#print(average_read_length_Tobago)

# Remove the first three rows from the data table
#average_read_length_Tobago <- average_read_length_Tobago[-c(1:3), ]

# Print the updated data table
#print(average_read_length_Tobago)


# Convert age_value to numeric
result_data_for_plot_PS97gesamt$age_value <- as.numeric(as.character(result_data_for_plot_PS97gesamt$age_value))

filtered_data <- result_data_for_plot_PS97gesamt %>%
  filter(age_value != 0)


# Plot boxplots over time with filtered data
ggplot(filtered_data, aes(x = factor(age_value), y = V2)) +
  geom_boxplot() +
  labs(x = "Age (yrs)", y = "Read Length") +
  scale_x_discrete() +
  theme_minimal()


plot_grid(XJoW004L_6_S3_plot3,
          XJoW004L_5_S2_plot3,
          XJoW004L_4_S2_plot3,
          XJoW004L_3_S2_plot3,
          XJoW004L_2_S2_plot3,
          XJoW004L_1_S2_plot3,
          XJoW003L_6_S2_plot3,
          XJoW003L_5_S2_plot3,
          XJoW003L_4_S2_plot3,
          XJoW003L_3_S1_plot3,
          XJoW003L_2_S1_plot3,
          XJoW003L_1_S1_plot3,
          XJoW002L_6_S1_plot3,
          XJoW002L_5_S1_plot3,
          XJoW005L_4_S3_plot3,
          XJoW002L_4_S1_plot3,
          XJoW005L_3_S3_plot3,
          #XJoW002L_3_S1_plot3,
          XJoW002L_2_S1_plot3,
          XJoW002L_1_S9_plot3,
          XJoW001L_3_S3_plot3,
          XJoW001L_1_S1_plot3,
          XJoW001L_6_S6_plot3,
          XJoW001L_5_S5_plot3,
          XJoW001L_4_S4_plot3,
          XJoW005L_2_S3_plot3,
          XJoW001L_2_S2_plot3,
          XJoW005L_1_S3_plot3,
          #XJoW001L_7_S7_plot3,
          XJoW001L_8_S8_plot3,
          XJoW002L_7_S1_plot3,
          #XJoW002L_8_S1_plot3,
          #XJoW003L_7_S2_plot3,
          XJoW003L_8_S2_plot3,
          #XJoW004L_7_S3_plot3,
          #XJoW004L_8_S3_plot3,
          XJoW005L_5_S3_plot3)


# First of all you have to read in following documents for KL77

# biodiv/user/joweiss/KL77/Age_KL77.txt
# biodiv/user/joweiss/KL77/INFO_Age_KL77.txt

#KL77
####### LOOP FOR READ COUNTS ################
library(data.table)

folder_path_info <- "/Volumes/projects/biodiv/user/shotgun_data/05.2_shotgun_data_nt05/KL77_NextSeqBHV/read_extract_Embryophyta"

# Get a list of all files in the folder
file_list_info_1 <- list.files(path = folder_path_info, full.names = TRUE)

# Filter the file list to include only files with "_merged_tax3193_extract_length" in their names
file_list_info_2 <- file_list_info_1[grep("_tax3193_extract_merged_count.txt", file_list_info_1)]

# Create an empty list to store the data tables for each file
table_list <- list()

# Loop through each file in the folder
table_list <- lapply(file_list_info_2, function(file_path) {
  # Read the data from the file (assuming it's a single value per file)
  value <- as.numeric(readLines(file_path))
  
  # Extract the file name without the path and extension
  file_name_info <- basename(file_path)
  
  # Create a data table with the file name and value
  data.table(file_name = file_name_info, value = value)
})

# Combine all the data tables into a single data table
result_table <- rbindlist(table_list, fill = TRUE)

# Print the resulting table
print(result_table)

# Write the table to a CSV file
output_file <- "/Volumes/projects/biodiv/user/joweiss/result_table_KL77.csv"
write.csv(result_table, file = output_file, row.names = FALSE)

#######LOOP FOR READ LENGTH ################
library(ggplot2)
library(tools)  
library(stringr)# Required for file_path_sans_ext function
library(gridExtra)



theme_set(theme(plot.title = element_text(size = 5)))



# Specify the folder path where the files are located
folder_path <- "/Volumes/projects/biodiv/user/shotgun_data/05.2_shotgun_data_nt05/KL77_NextSeqBHV/read_length_Embryophyta"

output_path <- "/Volumes/projects/biodiv/user/joweiss/KL77"
# Get a list of all files in the specified folder
file_list <- list.files(path = folder_path, full.names = TRUE)

# Filter the file list to include only files with "_merged_tax3193_extract_length" in their names
file_list_data <- file_list[grep("_merged_tax3193_extract_length", file_list)]

# Get a list of all files with "_merged_tax3193_extract_length" in the specified folder
#file_list <- list.files(path ="/Volumes/projects/biodiv/user/joweiss/KL12/output/out.kraken2/out.length_extract_Embryophyta", pattern = "_merged_tax3193_extract_length", full.names = TRUE)

# set NA to 0 in Age_KL12
Age_KL77$age[is.na(Age_KL77$age)] <- 0


# Set up an empty list to store the plots
plot_list <- list()

# MEDIAN Calculations 
medians <- list()

# Create an empty table to store the data
result_data_for_plot_KL77 <- data.frame(age_value = numeric(),
                                        V2 = numeric(),
                                        file_name_partial = character(),
                                        stringsAsFactors = FALSE)

# Loop through each file in the list
for (file in file_list_data) {
  # Check if the file is empty or has no content
  if (is.na(file.size(file)) || file.size(file) == 0) {
    # Skip this file and move on to the next iteration
    next
  }
  # Read the data from the file
  data <- read.table(file, header=FALSE)
  
  
  #extract the file name
  file_name_data <- file_path_sans_ext(basename(file))
  
  # Extract the first 29 characters of the file name
  file_name_partial <- substr(file_name_data, 1, 15)
  
  # Find the matching row in the results table based on the partial file name
  matching_rows <- result_table[stringdist::stringdist(file_name_partial, file_name, method = "jaccard", q = 2) >= 0.98, ]
  
  # Check if there are matching rows
  if (nrow(matching_rows) > 0) {
    # Extract the value from the first matching row
    info_value <- matching_rows$value[1]
  } else {
    # Set info_value to a default value or handle the case when there are no matches
    info_value <- "N/A"
  }
  
  # Find the matching row in the result table for Age_KL12 file
  #age_matching_row <- Age_KL77[str_detect(Age_KL77$file_name, file_name_partial), ]
  age_matching_rows <- Age_KL77[stringdist::stringdist(file_name_partial, Age_KL77$file_name, method = "jaccard", q = 2) >= 0.99, ]
  
  
  # Extract the Age value from the matching row
  age_value <- age_matching_row$age
  
  # Check if there are matching rows
  #  if (nrow(age_matching_row) > 0) {
  # Extract the value from the first matching row
  #   age_value <- age_matching_row$age[1]
  #  } else {
  # Set info_value to a default value or handle the case when there are no matches
  #   age_value <- "N/A"
  #}
  
  result_data_for_plot_KL77 <- rbind(result_data_for_plot_KL77, data.frame(age_value = age_value,
                                                                           V2 = data$V2,
                                                                           file_name_partial = file_name_partial,
                                                                           stringsAsFactors = FALSE))  
  
  
  # Create the ggplot objects
  plot1 <- ggplot(data, aes(x = V2)) +
    geom_histogram(binwidth = 0.5, colour = "black", fill = "white") +
    labs(title= paste("Age:", age_value),x= info_value)
  
  plot2 <- ggplot(data, aes(x = V2)) +
    geom_density() +
    labs(title= paste("Age:", age_value), x= info_value)
  
  # Generate a dynamic name for the plot
  plot_name <- paste0("X", age_value, gsub("-", "_", file_name_partial), "_plot3")
  
  # Create the plot with the dynamic name
  assign(plot_name, ggplot(data, aes(x = V2)) +
           geom_histogram(aes(y = after_stat(density)),
                          binwidth = 0.5,
                          colour = "black", fill = "white") +
           geom_density(alpha = 0.2, fill = "#FF6666") +
           labs(title = paste("Age:", age_value, gsub("-", "_", file_name_partial)), x = info_value))
  
  print(plot_name)
  
}

# Write the table to a CSV file
output_file <- "/Volumes/projects/biodiv/user/joweiss/KL77/result_data_for_plot_KL77.csv"
write.csv(result_data_for_plot_KL77, file = output_file, row.names = FALSE)

# Calculate the average read length over time
#average_read_length_Tobago <- result_data_for_plot_Tobago %>%
# group_by(age_value) %>%
#summarize(average_length = mean(V2))

# Print the average read length over time
#print(average_read_length_Tobago)

# Remove the first three rows from the data table
#average_read_length_Tobago <- average_read_length_Tobago[-c(1:3), ]

# Print the updated data table
#print(average_read_length_Tobago)


# Convert age_value to numeric
result_data_for_plot_KL77$age_value <- as.numeric(as.character(result_data_for_plot_KL77$age_value))

filtered_data <- result_data_for_plot_KL77 %>%
  filter(age_value != 0)


# Plot boxplots over time with filtered data
ggplot(filtered_data, aes(x = factor(age_value), y = V2)) +
  geom_boxplot() +
  labs(x = "Age (yrs)", y = "Read Length") +
  scale_x_discrete() +
  theme_minimal()

ggplot(filtered_data, aes(x = factor(age_value), y = V2)) +
  geom_boxplot() +
  labs(x = "Age (yrs)", y = "Read Length") +
  scale_x_discrete(labels = function(x) sprintf("%.2f", as.numeric(x) / 1000)) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 6))

plot_grid(XStB106L_01_S9_plot3, 
          XJK136L_1_S1_m_plot3, 
          XStB106L_02_S10_plot3, 
          XStB106L_03_S11_plot3,
          XJK136L_2_S2_m_plot3, 
          XJK136L_3_S3_m_plot3,
          XStB106L_04_S12_plot3, 
          XStB106L_05_S13_plot3, 
          XJK136L_4_S4_m_plot3, 
          XStB097L_01_S1_plot3, 
          XStB097L_02_S2_plot3, 
          XJK136L_5_S5_m_plot3, 
          XStB097L_03_S3_plot3, 
          XJK136L_6_S6_m_plot3, 
          XStB097L_04_S4_plot3, 
          XStB097L_05_S5_plot3, 
          XStB097L_06_S6_plot3, 
          XStB106L_07_S15_plot3, 
          XJK136L_7_S7_m_plot3, 
          XStB108L_01_S17_plot3, 
          XJK137L_1_S9_m_plot3,
          XJK137L_2_S10__plot3, 
          #StB108L_02_S18_plot3, 
          XJK137L_3_S11__plot3, 
          XStB108L_03_S19_plot3, 
          #XJK137_4_S12__plot3, 
          XJK137L_5_S13__plot3, 
          XStB108L_05_S21_plot3, 
          XJK137L_6_S14__plot3, 
          XJK137L_7_S15__plot3, 
          XStB109L_01_S25_plot3, 
          XJK138L_1_S17__plot3, 
          #XJK138L_2_S18__plot3, 
          XJK138L_3_S19__plot3, 
          XStB109L_02_S26_plot3,
          XJK138L_4_S20__plot3, 
          XStB109L_04_S28_plot3,
          XStB109L_05_S29_plot3,
          XJK138L_6_S22__plot3, 
          XStB109L_06_S30_plot3,
          XJK138L_7_S23__plot3, 
          XStB108L_07_S23_plot3, 
          XStB097L_07_S7_plot3, 
          XStB097L_08_S8_plot3, 
          XStB106L_06_S14_plot3, 
          XStB106L_08_S16_plot3, 
          XStB108L_04_S20_plot3,
          XStB108L_06_S22_plot3, 
          XStB108L_08_S24_plot3,
          #XStB109L_03_S27_plot3, 
          XStB109L_07_S31_plot3, 
          XStB109L_08_S32_plot3)
          #XJK136L_8_S8__plot3, 
          #XJK137L_8_S16__plot3, 
          #XJK138L_8_S24__plot3)
