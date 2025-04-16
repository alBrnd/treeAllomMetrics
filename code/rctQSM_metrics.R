
library(ITSMe)


# Functions for volume calculation
trunk_volume_rct <- function(cylinder, cylindercutoff = 0) {
  sum(cylinder$volume[cylinder$BranchOrder == 0 & cylinder$diameter > cylindercutoff])
}

branch_volume_rct <- function(cylinder, cylindercutoff = 0) {
  sum(cylinder$volume[cylinder$BranchOrder > 0 & cylinder$diameter > cylindercutoff])
}

tree_volume_rct <- function(treedata, cylinder = NA, cylindercutoff = 0) {
  if (cylindercutoff > 0 & length(cylinder) > 1) {
    sum(cylinder$volume[cylinder$diameter > cylindercutoff])
  } else {
    treedata$TotalVolume[1]
  }
}



# Main function
extract_rct_qsm_metrics <- function(inpath, outpath) {
  
  # List all info.txt files
  qsmfiles <- list.files(inpath, pattern = "info\\.txt$", full.names = TRUE)
  
  # Process all files
  result_list <- lapply(qsmfiles, function(file) {
    
    qsm <- read_rct_qsm(file)
    treedata <- qsm$treedata
    cylinder <- qsm$cylinder
    
    data.frame(
      treeid = as.numeric(gsub("_", "", regmatches(file, regexpr("_\\d+", file)))),
      treelabel = sub(".*trees/(.*)_filtered.*", "\\1", file),
      filename = basename(file),
      DBH_rct = treedata$DBHqsm,
      height_rct = treedata$TreeHeight,
      crown_radius_rct = treedata$CrownRadius,
      total_vol_rct = treedata$TotalVolume/1000,
      vol_cut0.07_rct = tree_volume_rct(treedata, cylinder, cylindercutoff = 0.07),
      vol_cut0.02_rct = tree_volume_rct(treedata, cylinder, cylindercutoff = 0.02),
      vol_trunk_rct = trunk_volume_rct(cylinder, cylindercutoff = 0.07),
      vol_trunk_cut0.07_rct = trunk_volume_rct(cylinder),
      vol_branch_cut0.07_rct = branch_volume_rct(cylinder, cylindercutoff = 0.07),
      vol_branch_cut0.02_rct = branch_volume_rct(cylinder, cylindercutoff = 0.02),
      vol_branch_all_rct = branch_volume_rct(cylinder, cylindercutoff = 0),
      sum_length_branch_cut0.07_rct = sum(cylinder$length[cylinder$BranchOrder > 0 & cylinder$diameter > 0.07]),
      mean_length_branch_cut0.07_rct = mean(cylinder$length[cylinder$BranchOrder > 0 & cylinder$diameter > 0.07]),
      sum_length_branch_all_rct = sum(cylinder$length[cylinder$BranchOrder > 0]),
      mean_length_branch_all_rct = mean(cylinder$length[cylinder$BranchOrder > 0]),
      sum_length_totalcylinders_rct = sum(cylinder$length),
      mean_length_totalcylinders_rct = mean(cylinder$length),
      sum_length_totalcylinders_cut0.07_rct = sum(cylinder$length[cylinder$diameter > 0.07])
    )
  })
  
  # Combine results
  result_df <- do.call(rbind, result_list)
  
  # Write to file
  write.csv(result_df, file.path(outpath, "RCT_QSM_metrics.csv"), row.names = FALSE)
  
  cat("File saved to:", file.path(outpath, "RCT_QSM_metrics.csv"), "\n")
  
  return(result_df)
}


# data <- extract_rct_qsm_metrics(
#   inpath = "F:/CRS-Allometries-Datasets/forQSM/trees",
#   outpath = "F:/CRS-Allometries-Datasets/forQSM")
# 
# # volume in m^3
# # length in cm (?)


