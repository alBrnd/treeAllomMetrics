# preprocessing single tree point clouds: 
# - downsample to 0.02m
# - remove small point clusters
# - remove point clusters far away from largest cluster
# caveat: clustering can be very slow for larger point clouds

# if point cloud too big (> 100 000 points): 
# - remove points with negative Z
# - clip point cloud by bounding box of 40x40 m around center point (median)

library(lidR)
library(Rvcg)

# Function to remove small clusters
remove_small_clusters <- function(point_cloud, min_cluster_size = 10) {
  cluster_sizes <- table(point_cloud$cluster)
  valid_clusters <- names(cluster_sizes[cluster_sizes >= min_cluster_size])
  return(point_cloud[point_cloud$cluster %in% valid_clusters, ])
}

# Function to remove distant clusters
remove_distant_clusters <- function(point_cloud, max_distance = 5) {
  cluster_sizes <- table(point_cloud$cluster)
  largest_cluster <- names(which.max(cluster_sizes))
  largest_cluster_points <- point_cloud[point_cloud$cluster == largest_cluster, ]
  
  # Compute distance of smaller clusters from the largest cluster
  cluster_distances <- sapply(unique(point_cloud$cluster), function(cluster) {
    if (cluster == largest_cluster) return(0)
    
    cluster_points <- point_cloud[point_cloud$cluster == cluster, ]
    largest_cluster_center <- colMeans(largest_cluster_points[, c("X", "Y", "Z")])
    mean_dist <- mean(VoxR::point_distance(cluster_points, largest_cluster_center))
    return(mean_dist)
  })
  
  # Keep only clusters within max_distance from the largest cluster
  valid_clusters <- names(cluster_distances[cluster_distances <= max_distance])
  valid_clusters <- unique(point_cloud$cluster)[cluster_distances <= max_distance]
  return(point_cloud[point_cloud$cluster %in% valid_clusters, ])
}

# Function to filter point cloud by bounding box
box_filter <- function(pointcloud) {
  # Step 1: Remove points with negative Z values
  pointcloud <- pointcloud[ , .SD[Z >= 0], .SDcols = names(pointcloud)]
  
  # Step 2: Calculate the median center
  median_x <- median(pointcloud[[1]])
  median_y <- median(pointcloud[[2]])
  median_z <- median(pointcloud[[3]])
  
  # Step 3: Keep points within 20m in X/Y and 40m in Z from median center
  pointcloud <- pointcloud[
    abs(get(names(pointcloud)[1]) - median_x) <= 20 &
      abs(get(names(pointcloud)[2]) - median_y) <= 20 &
      abs(get(names(pointcloud)[3]) - median_z) <= 40
  ]
  
  return(pointcloud)
}


# function to save point clouds as .ply
save_as_ply <- function(pointcloud, output_file) {
  
  # Extract XYZ
  xyz <- as.matrix(pointcloud)
  
  # Create vcg object
  vcg_obj <- list(
    vb = t(rbind(xyz, rep(1, nrow(xyz)))),  # Homogeneous coords (4 x n)
    it = matrix(numeric(0), nrow = 3),      # No faces = point cloud only
    material = list()
  )
  class(vcg_obj) <- "mesh3d"
  
  # Write binary PLY
  Rvcg::vcgPlyWrite(vcg_obj, filename = output_file, binary = TRUE)
  
  #cat(sprintf("Converted %s to %s\n", input_file, output_file))
}



# Main processing function
filter_by_cluster <- function(pcpath, outpath_ply, outpath_laz) {
  point_cloud <- ITSMe::read_tree_pc(path = pcpath)
  point_cloud <- data.table::data.table(X=point_cloud$X, Y=point_cloud$Y, Z=point_cloud$Z)
  density_filtered <- VoxR::filter_point_density(point_cloud,0.02)
  sor_filtered_cloud <- VoxR::filter_noise(density_filtered, k=5, sigma=1.5, store_noise=FALSE)
  
  if(nrow(sor_filtered_cloud) > 100000){
    final_cloud <- box_filter(sor_filtered_cloud)  
  }else{
  
  clustering <- VoxR::distance_clustering(sor_filtered_cloud, d_clust=2) # d_clust = The distance required to consider two points as being part of two different clusters
  clust_filtered <- remove_small_clusters(clustering, min_cluster_size=5) # min_cluster_size= min number of points per cluster to keep it
  final_cloud <- if (length(unique(clust_filtered$cluster)) > 1) remove_distant_clusters(clust_filtered, max_distance = 5) else clust_filtered
  }
  
  # center point cloud by subtracting minimum coordinates
  min_coords <- apply(final_cloud, 2, min)
  centered_final_cloud <- sweep(final_cloud, 2, min_coords, FUN = "-")
  
  # save filtered pointcloud as txt
  #write.table(centered_final_cloud[,1:3], file = file.path(outpath, paste0(sub('\\..*$', '', basename(pcpath)), ".txt")), sep = " ", row.names = FALSE, col.names = FALSE)
  
  # save filtered point cloud as .ply
  save_as_ply(centered_final_cloud, file.path(outpath_ply, paste0(sub('\\..*$', '', basename(pcpath)), "_filtered.ply")))
  
  # save filtered pointcloud as laz
  laspc <- lidR::LAS(centered_final_cloud[,1:3])
  lidR::writeLAS(laspc, file = file.path(outpath_laz, paste0(sub('\\..*$', '', basename(pcpath)), "_filtered.laz")))
  
  #return(final_cloud)
}







# filepath <- "D:/STSMGent/data/ray/trees/Hajnowka_brz_44_m142.laz"
# filepath <- "D:/STSMGent/data/Bialowieza_sw_71_testclusters2.las"
# las <- lidR::readLAS(filepath)
# xyz <- data.table::data.table(x=las@data$X, y=las@data$Y, z=las@data$Z)
# #or:
# pc <- ITSMe::read_tree_pc(path = filepath)
# xyz <- data.table::data.table(x=pc$X, y=pc$Y, z=pc$Z)
# 
# pc_out <- filter_by_cluster(xyz) # or pc


