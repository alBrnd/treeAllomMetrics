# full workflow for calculating metrics (for allometric tree volume/AGB models) from single tree point clouds
# part 2

source(file.path(path_to_src,"read_rct_qsm.R"))
source(file.path(path_to_src,"rctQSM_metrics.R"))
source(file.path(path_to_src,"plot_rctQSM.R"))

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
