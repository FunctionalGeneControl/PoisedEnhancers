library(data.table)
setDTthreads(4)

# Get command-line arguments
args = commandArgs(trailingOnly = TRUE)

# Arguments: file with list of input files and output file
input_list_file <- args[1]
output_file <- args[2]

# Read the list of input files with aliases
input_list <- fread(input_list_file, header = FALSE, sep = "=", col.names = c("alias", "file_path"))
input_list[, file_path := trimws(file_path)]  # Remove extra spaces around file paths
input_list[, alias := trimws(alias)]          # Remove extra spaces around aliases

# Print the input files and aliases
cat("Input files with aliases are:\n")
print(input_list)

cat("Output read counts file is:\n", output_file, "\n")

# Prepare an empty list to store results
results_list <- vector("list", length(input_list$file_path))

for (i in seq_along(input_list$file_path)) {
  
  # Read the file
  input_data <- fread(input_list$file_path[i])
  
  # Filter rows where distSign is not empty
  filtered_data <- input_data[distSign != "" & !is.na(distSign)]
  
  # Summarize by baitID
  summarized_data <- filtered_data[, .(N = sum(N)), by = baitID]
  
  # Rename the "N" column to use the alias as the column name
  setnames(summarized_data, "N", input_list$alias[i])
  
  # Store the summarized result in the list
  results_list[[i]] <- summarized_data
}

# Merge all the data.tables on baitID, filling missing values with NA
wide_format_data <- Reduce(function(x, y) merge(x, y, by = "baitID", all = TRUE), results_list)

# Replace empty cells with NA explicitly (ensures compatibility)
wide_format_data[wide_format_data == ""] <- NA

# Write the result to the output file
fwrite(wide_format_data, output_file, sep = "\t", na = "NA")

cat("Wide-format table created. Output file is:", output_file, "\n")
