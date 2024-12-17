#Read Length 
setwd("/Volumes/projects/biodiv/user/joweiss")

# Import tax IDs, and Age files
tax_report_ID <- read.delim("~/Desktop/Weiss,Stoof-Leichsenring,Herzschuh (2022)/DATA/TAX_ID_REAL_FINAL_1.txt")
Age_KL12 <- read.delim("/Volumes/projects/biodiv/user/joweiss/KL12/output/out.kraken2/Age_KL12.txt")
Age_KL77 <- read.delim("/Volumes/projects/biodiv/user/joweiss/KL77/Age_KL77.txt")
Age_MSM <- read.delim("/Volumes/projects/biodiv/user/joweiss/MSM/Age_MSM.txt")
Age_PS97 <- read.delim("/Volumes/projects/biodiv/user/joweiss/PS97/Age_PS97.txt")
Age_Tobago <- read.delim("/Volumes/projects/biodiv/user/joweiss/Tobago/Age_Tobago.txt")

#KL12, PS97, KL77, Tobago, Fram Strait
####### LOOP FOR READ COUNTS KL12 ################
library(data.table)

folder_path_info <- "/Volumes/projects/biodiv/user/joweiss/KL12/out.read_extract_33090" #33090 sind embryophyta

# Get a list of all files in the folder
file_list_info_1 <- list.files(path = folder_path_info, full.names = TRUE)

# Filter the file list to include only files with "merged_tax33090_separated_extract_length" in their names
file_list_info_2 <- file_list_info_1[grep("_tax33090_extract_merged_count.txt", file_list_info_1)]

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
output_file <- "/Volumes/projects/biodiv/user/joweiss/KL12/out.read_extract_Embryophyta/result_table_KL12_with_taxID.csv"
write.csv(result_table, file = output_file, row.names = FALSE)



####### LOOP FOR READ LENGTH KL12 ################
library(ggplot2)
library(tools)  
library(stringr)# Required for file_path_sans_ext function
library(gridExtra)

theme_set(theme(plot.title = element_text(size = 6)))


# Specify the folder path where the files are located
folder_path <- "/Volumes/projects/biodiv/user/joweiss/KL12/out.length_extract_33090"

output_path <- "/Volumes/projects/biodiv/user/joweiss/KL12/output/out.kraken2/PLOTS"
# Get a list of all files in the specified folder
file_list <- list.files(path = folder_path, full.names = TRUE)

# Filter the file list to include only files with "merged_tax33090_separated_extract_length" in their names
file_list_data <- file_list[grep("merged_tax33090_separated_extract_length", file_list)]

# Get a list of all files with "merged_tax33090_separated_extract_length" in the specified folder
#file_list <- list.files(path ="/Volumes/projects/biodiv/user/joweiss/KL12/output/out.kraken2/out.length_extract_Embryophyta", pattern = "merged_tax33090_separated_extract_length", full.names = TRUE)

# set NA to 0 in Age_KL12
Age_KL12$age[is.na(Age_KL12$age)] <- 0

# Set up an empty list to store the plots
plot_list <- list()

# MEDIAN Calculations 
medians <- list()

# Create an empty table to store the data
result_data_for_plot_KL12 <- data.frame(age_value = numeric(),
                                        Read_Length = numeric(),
                                        taxid = character(),
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
                                                                           Read_Length = data$V3,
                                                                           taxid = data$V2,
                                                                           file_name_partial = file_name_partial,
                                                                           stringsAsFactors = FALSE))  
  # Create the ggplot objects
  plot1 <- ggplot(data, aes(x = Read_Length)) +
    geom_histogram(binwidth = 0.5, colour = "black", fill = "white") +
    labs(title= paste("Age:", age_value),x= info_value)
  
  plot2 <- ggplot(data, aes(x = Read_Length)) +
    geom_density() +
    labs(title= paste("Age:", age_value), x= info_value)
  
  # Generate a dynamic name for the plot
  plot_name <- paste0("X", gsub("-", "_", file_name_partial), "_plot3")
  
  # Create the plot with the dynamic name
  assign(plot_name, ggplot(data, aes(x = Read_Length)) +
           geom_histogram(aes(y = after_stat(density)),
                          binwidth = 0.5,
                          colour = "black", fill = "white") +
           geom_density(alpha = 0.2, fill = "#FF6666") +
           labs(title = paste("Age:", age_value, gsub("-", "_", file_name_partial)), x = info_value))
  
  print(plot_name)
  
}

print(result_data_for_plot_KL12)


# Load required packages
library(dplyr)

# Merge the tables by Tax_ID
merged_table <- merge(result_data_for_plot_KL12, tax_report_ID, by = "taxid")

# Remove rows with age_value equal to zero and NA
merged_table <- merged_table[!(merged_table$age_value == 0 | is.na(merged_table$age_value)), ]

# Print the merged table
print(merged_table)

unique(merged_table$family)

# Load required packages
library(ggplot2)

# Get all unique families
all_families <- unique(merged_table$family)

# Calculate the median read length for each family
family_median <- merged_table %>%
  group_by(family) %>%
  summarise(median_read_length = median(Read_Length)) %>%
  arrange(desc(median_read_length))

unique(tax_report_ID$taxid)

# Reorder the Family column based on median read length
merged_table$family <- factor(merged_table$family, levels = family_median$family
)

# Plotting the read length for each family as individual boxplots
ggplot(merged_table, aes(x = family, y = Read_Length)) +
  geom_boxplot() +
  xlab("Family") +
  ylab("Read Length") +
  ggtitle("Read Length for Each Family") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed