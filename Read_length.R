# Load necessary libraries
library(ggplot2)
library(gridExtra)

# Assuming "combined_cores_ts" is your data frame containing the data

# Step 1: Calculate the median read length for each family
median_read_length <- aggregate(Read_Length ~ family, data = combined_cores_ts, FUN = median)

# Step 2: Sort the families based on the median read length (highest to lowest)
sorted_families <- median_read_length[order(-median_read_length$Read_Length), "family"]

# Step 3: Create 6 subsets of data based on the ordered families
num_families <- length(sorted_families)
num_subplots <- 6
families_per_plot <- ceiling(num_families / num_subplots)

# Step 4: Create a list to store the ggplot objects for each subset
plot_list <- list()

# Step 5: Plot the data in 6 separate plots and store them in the list
for (i in 1:num_subplots) {
  start_idx <- (i - 1) * families_per_plot + 1
  end_idx <- min(i * families_per_plot, num_families)
  current_families <- sorted_families[start_idx:end_idx]
  
  # Filter the data to include only the current families
  current_data <- combined_cores_ts[combined_cores_ts$family %in% current_families, ]
  
  # Plot the boxplot for the current subset of data and store it in the list
  p <- ggplot(current_data, aes(x = family, y = Read_Length)) +
    geom_boxplot() +
    ylab("Read length (bp)") +
    #ggtitle(paste("Read Length for Families ", start_idx, " to ", end_idx, sep = "")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # Smaller font for family names
          axis.title.x = element_blank(),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"),
          # Remove x-axis title
          panel.background = element_rect(fill = "white"),  # Set background to white
          plot.margin = unit(c(1, 1, 1, 1), "cm"))  # Add margins to the plot to make space between plots
  
  plot_list[[i]] <- p
}

# Step 6: Assemble four plots together in one plot, ordered from top left to bottom right
final_plot <- grid.arrange(grobs = plot_list, ncol = 2)

# Step 7: Save the final plot
ggsave(filename = "final_plot.png", plot = final_plot, width = 12, height = 10, units = "in", dpi = 300)
