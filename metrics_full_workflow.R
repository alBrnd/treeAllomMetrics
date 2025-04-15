# full workflow for calculating metrics (for allometric tree volume/AGB models) from single tree point clouds

# install all required packages
install.packages(c("data.table", "ITSMe", "doParallel", "foreach", "dplyr", "patchwork",  "VoxR", "lidR", "Rvcg", "rgl", "TreeLS"))


# set paths
# main working directory
wd_path <- "F:/CRS-Allometries-Datasets"
setwd(wd_path)

# load code
path_to_src <- "./code"
source(file.path(path_to_src,"preprocess_filtering_functions.R"))
#source(file.path(path_to_src,"basic_metrics_pc_v2.R"))
source(file.path(path_to_src,"basic_metrics_pc_v2_ds.R"))
source(file.path(path_to_src,"summary_metrics_v2.R"))
source(file.path(path_to_src,"read_rct_qsm.R"))
source(file.path(path_to_src,"rctQSM_metrics.R"))
source(file.path(path_to_src,"plot_rctQSM.R"))
source(file.path(path_to_src,"tapering_metrics.R"))
source(file.path(path_to_src,"stem_taper_treels_function.R"))

#######################################################################################
# step 1: renaming files
# (adapt according to new folder structure!)

# rename point cloud files based on originalpclabel and newpointcloudname in database

# Read csv that contains info from database
df <- read.csv("database_refdata.csv", stringsAsFactors = FALSE)

# List of deliverers (adapt according to folder names)
deliverers <- c("Abegg") #, "Krucek", "Terryn"

# Base directory where the folders are located
base_dir <- getwd()

# Directory to save renamed files
output_dir <- "./renamed_pointclouds/"
if (!dir.exists(output_dir)) {dir.create(output_dir)}

# Function to search for a file recursively
find_file <- function(directory, filename) {
  files <- list.files(directory, pattern = filename, recursive = TRUE, full.names = TRUE)
  if (length(files) > 0) {
    return(files[1])  # Return the first match found
  } else {
    return(NA)
  }
}

# Loop through each deliverer
for (deliverer in deliverers) {
  
  # Filter dataframe for the current deliverer
  newdf <- df %>% filter(familyname == deliverer)
  
  # Loop through each row in the filtered dataframe
  for (i in 1:nrow(newdf)) {
    
    # Construct source folder name
    sourcefolder <- file.path(base_dir, paste0(newdf$name[i], newdf$familyname[i]))
    
    # Search for the old file in source folder and subfolders
    old_file <- find_file(sourcefolder, newdf$originalpclabel[i])
    
    # Define new file path
    new_file <- file.path(output_dir, newdf$newpointcloudname[i])
    
    # Check if the old file exists and copy it with a new name
    if (!is.na(old_file) && file.exists(old_file)) {
      if (!file.exists(new_file)) {
        file.copy(old_file, new_file)
        cat("Copied and renamed:", old_file, "->", new_file, "\n")
      } else {
        cat("File already exists:", new_file, "\n")
      }
    } else {
      cat("File not found:", newdf$originalpclabel[i], " in", sourcefolder, "\n")
    }
  }
}


#######################################################################################
# step 2: filter and preprocess point clouds

library(foreach)
library(doParallel)

inpath <- "./data/renamed_pointclouds"

outpath_laz <- "./data/renamed_pointclouds_filtered"  
if (!dir.exists(outpath_laz)) {dir.create(outpath_laz)}

outpath_ply <- "./data/rct_qsm/trees"  
if (!dir.exists(outpath_ply)) {dir.create(outpath_ply)}

infiles <- list.files(inpath, full.names = TRUE)  

# Existing filenames in outpath
existing_files <- list.files(outpath_laz)

# Input file base names without extension
input_basenames <- sub('\\..*$', '', basename(infiles))

# Filter: keep only those not found in existing_files
files_to_process <- infiles[
  !sapply(input_basenames, function(x) any(grepl(x, existing_files)))
]


# parallel processing with library doParallel

# set number of cores to use
nr_cores <- detectCores() - 1 # use as many as possible
#or:
nr_cores <- 6   # N of cores to use; do not use all cores available

cl <- makeCluster(spec = nr_cores)
registerDoParallel(cl)

n <- length(files_to_process)

foreach(i = 1:n, .packages = c("data.table", "ITSMe", "VoxR", "lidR", "Rvcg") ) %dopar%  {
  
  filter_by_cluster(pcpath = files_to_process[i], outpath_ply, outpath_laz)
}

stopCluster(cl)
registerDoSEQ()

print("filtering (step 2) finished")

#######################################################################################
# step 3: calculate basic metrics with ITSMe package

library(ITSMe)
library(patchwork)
library(foreach)
library(doParallel)
library(dplyr)


inpath <- "./data/renamed_pointclouds_filtered"
outpath <- "./data/renamed_pointclouds_metrics/"
if (!dir.exists(outpath)) {dir.create(outpath)}
outpath_csv <- paste0(outpath, "/metrics_per_tree")
if (!dir.exists(outpath_csv)) {dir.create(outpath_csv)}


infiles <- list.files(inpath, full.names = TRUE)  

# Existing filenames in outpath
existing_files <- list.files(outpath_csv)

# Input file base names without extension
input_basenames <- sub('\\..*$', '', basename(infiles))

# Filter: keep only those not found in existing_files
files_to_process <- infiles[
  !sapply(input_basenames, function(x) any(grepl(x, existing_files)))
]


# calculte metrics with summary_basic_pointcloud_metrics_pertree

# parallel with library doParallel
nr_cores <- detectCores() - 1
#or:
#nr_cores <- 6   # N of cores to use; do not use all cores available
cl <- makeCluster(spec = nr_cores)
registerDoParallel(cl)

n <- length(files_to_process)

foreach(i = 1:n, .packages = c("ITSMe", "patchwork") ) %dopar%  {
  
  tryCatch({
    
    pc_metrics <- summary_basic_pointcloud_metrics_pertree(files_to_process[i], metrics = c(
      "stem diameter",
      "tree height",
      "projected area",
      "alpha volume"), # "alpha volume"
      taper = FALSE,
      crown = TRUE,
      plot = TRUE,
      OUT_path = outpath)
    
    pc_metrics$treeid <- as.numeric(gsub("_", "", regmatches(pc_metrics$tree_id, regexpr("_\\d+", pc_metrics$tree_id))))
    
    utils::write.csv(pc_metrics,
                     paste0(outpath_csv, "/basic_metrics_", sub('\\..*$', '', basename(files_to_process[i])),".csv"),
                     row.names = FALSE)  
  })
}

stopCluster(cl)
registerDoSEQ()


## combine metrics into one dataframe

# List all CSV files in the folder
csv_files <- list.files(path = outpath_csv, pattern = "*.csv", full.names = TRUE)
# Initialize an empty list to store data frames
data_list <- list()
# Loop through each CSV file, read it, and add a filename column
for (file in csv_files) {
  data <- read.csv(file)
  data$filename <- basename(file)  # Add a column with the filename
  data_list[[length(data_list) + 1]] <- data
}
# Combine all data frames into one
combined_data <- bind_rows(data_list)
# Write the combined data into a new CSV file
write.csv(combined_data, file = "summary_basic_metrics.csv", row.names = FALSE)

#######################################################################################
# step 3.2: calculate diameters and stem volume with TreeLS package

library(foreach)
library(doParallel)

inpath <- "./data/renamed_pointclouds_filtered"
outpath <- "./data/renamed_pointclouds_stemmetrics/"
if (!dir.exists(outpath)) {dir.create(outpath)}

infiles <- list.files(inpath, full.names = TRUE)  

# Existing filenames in outpath
existing_files <- list.files(outpath)

# Input file base names without extension
input_basenames <- sub('\\..*$', '', basename(infiles))

# Filter: keep only those not found in existing_files
files_to_process <- infiles[
  !sapply(input_basenames, function(x) any(grepl(x, existing_files)))
]



# parallel with library doParallel
nr_cores <- detectCores() - 1
#or:
#nr_cores <- 6   # N of cores to use; do not use all cores available
cl <- makeCluster(spec = nr_cores)
registerDoParallel(cl)

n <- length(files_to_process)

foreach(i = 1:n, .packages = c("dplyr", "ggplot2", "lidR", "TreeLS") ) %dopar%  {
  
  tryCatch({
    
    dia_vol <- calculate_stem_volume(file=files_to_process[i],
      outpath = outpath, slice_length = 0.2)
    
  })
}

stopCluster(cl)
registerDoSEQ()


## combine metrics into one dataframe

# List all CSV files in the folder
csv_files <- list.files(path = outpath, pattern = "stemVolumes_circle_irls", full.names = TRUE)
# Initialize an empty list to store data frames
data_list <- list()
# Loop through each CSV file, read it, and add a filename column
for (file in csv_files) {
  data <- read.csv(file)
  data$filename <- basename(file)  # Add a column with the filename
  data_list[[length(data_list) + 1]] <- data
}
# Combine all data frames into one
combined_data <- bind_rows(data_list)
# Write the combined data into a new CSV file
write.csv(combined_data, file = "stem_dia_vol_metrics.csv", row.names = FALSE)

#######################################################################################
# step 4: generate QSMs with raycloudtools

# install raycloudtools (if not already done). instructions: https://github.com/qforestlab/rayextract-manual

# run batch script
# it assumes a folder called "trees" containing the .ply files to process
# cd ./data/rct_qsm
# bash batch_run_raycloudQSMs.sh


#######################################################################################
# step 5: extract metrics from raycloud QSMs

data <- extract_rct_qsm_metrics(
  inpath = "./data/rct_qsm/trees",
  outpath = "./")

# volume in m^3
# length in cm (?)


# optional: plot all QSMs
library(rgl)
library(Rvcg)
library(magick)

ply_grapher(file_path= "./data/rct_qsm/trees", color="#f6b251", out_path = "./data/rct_qsm/qsm_images")


#######################################################################################
# step 6: merge all metrics into one dataframe

basic_metrics <- read.csv("summary_basic_metrics.csv", stringsAsFactors = FALSE)
stem_metrics <- read.csv("stem_dia_vol_metrics.csv", stringsAsFactors = FALSE)
rct_metrics <- read.csv("RCT_QSM_metrics.csv", stringsAsFactors = FALSE)

joined_df <- left_join(basic_metrics, stem_metrics, by="treeid")
joined_df <- left_join(joined_df, rct_metrics, by="treeid")
write.csv(joined_df, file = "CRS_metrics.csv", row.names = FALSE)




