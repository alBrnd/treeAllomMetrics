# stem volume

# Load required libraries
  library(dplyr)
  library(ggplot2)
  library(lidR)
  library(TreeLS)

calculate_stem_volume <- function(file, outpath, slice_length = 0.2) {
  
  outpath_plots <- file.path(outpath, "images")
  if (!dir.exists(outpath_plots)) {dir.create(outpath_plots)}

  # Read the point cloud file and convert to LAS object
  message("Processing file: ", file)
  tree_cloud <- readTLS(file)
  
  treeid = as.numeric(gsub("_", "", regmatches(file, regexpr("_\\d+", file))))
  
  # Normalize Z-axis
  minZ <- min(tree_cloud@data[["Z"]])
  tree_cloud@data[["Z"]] <- tree_cloud@data[["Z"]] - minZ
  
  # Classify stem points
  tree_stem <- stemPoints(tree_cloud, stm.hough())
  stem_points <- filter_poi(tree_stem, Stem == TRUE)
  
  # Define sequence of heights for stem slices
  height_seq <- seq(from = min(stem_points@data[["Z"]]), 
                    by = slice_length, 
                    to = floor(max(stem_points@data[["Z"]])))
  
  # Pre-allocate a list for slices
  slice_list <- vector("list", length(height_seq))
  
  # Slice stem point cloud and fit circles
  for (j in seq_along(height_seq)) {
    segment <- filter_poi(stem_points, Z > height_seq[j] - slice_length * 0.75 & Z < height_seq[j] + slice_length * 0.75)
    
    if (length(segment@data[["Z"]]) < 50) {
      slice_list[[j]] <- data.frame(X = NA, Y = NA, Radius = NA, Error = NA)
    } else {
      slice_list[[j]] <- shapeFit(segment, shape = "circle", algorithm = "irls")
    }
  }
  
  # Combine slices into a data frame
  slices <- do.call(rbind, slice_list)
  segments <- cbind(slices, height = height_seq[1:nrow(slices)])
  segments$dia <- segments$Radius * 2
  
  # Calculate rough stem volume
  segments$vol.raw <- pi * (segments$dia / 2)^2 * slice_length
  seg.vol.raw <- sum(segments$vol.raw, na.rm = TRUE)
  
  # Clean and smooth diameter data
  segments$dia[segments$Radius > mean(segments$Radius[1:4]) + 0.1] <- NA
  segments <- segments[!is.na(segments$height), ]
  
  NonNAindex <- which(!is.na(segments$dia))
  firstNonNA <- min(NonNAindex)
  segments <- segments[firstNonNA:length(segments$height), ]
  
  if(length(segments$Radius)>4){
  
    for (k in 4:length(segments$height)) {
      if (is.na(segments$dia[k]) ||
          abs(mean(c(segments$dia[k - 3], segments$dia[k - 2], segments$dia[k - 1]), na.rm = TRUE) - segments$dia[k]) > 0.1 ||
          is.na(mean(c(segments$dia[k - 3], segments$dia[k - 2], segments$dia[k - 1]), na.rm = TRUE))) {
        segments$dia[k] <- segments$dia[k - 1] - 0.01
      }
    }
    
    for (k in 10:length(segments$height)) {
      if (is.na(segments$Radius[k]) && is.na(segments$Radius[k-1]) && is.na(segments$Radius[k-2]) && is.na(segments$Radius[k-3])) {
        segments$dia[k] <- NA
      }
    }
  }  
  
  # Fit a cubic smoothing spline
  segments.comp <- segments[!is.na(segments$dia), ]
  
  if(length(segments.comp$Radius) > 4){
  
    sm <- stats::smooth.spline(segments.comp$height, segments.comp$dia, spar = 0.7)
    
    dia.interp <- predict(sm, segments$height)
    segments$dia.interp <- dia.interp$y
    segments <- dplyr::filter(segments, dia.interp >= 0.07)
    
    # Calculate smoothed stem volume
    segments$vol <- pi * (segments$dia.interp / 2)^2 * slice_length
    seg.vol.smooth <- sum(segments$vol)
    
    # write segments dataframe
    csvsegments <- file.path(outpath, paste0("stemDia_circle_irls_", basename(file), ".csv"))
    write.csv(segments, csvsegments, row.names = FALSE)
    
    # Extract diameter at height 7 m
    if (any(abs(segments$height - 7) <= slice_length / 2)) {
      D7 <- segments$dia.interp[which.min(abs(segments$height - 7))]
    } else {
      D7 <- NA  # If no data point is close to 7 m, set D7 to NA
    }
    
    # Extract diameter at height 1.3 m
    if (any(abs(segments$height - 1.3) <= slice_length / 2)) {
      DBH <- segments$dia.interp[which.min(abs(segments$height - 1.3))]
    } else {
      DBH <- NA  # If no data point is close to 1.3 m, set D7 to NA
    }
    
    
    # Generate stem plot with a legend
    stem_plot <- ggplot() +
      # Add segment diameter points (black)
      geom_point(aes(x = segments$height, y = segments$dia, color = "Segment Diameter"), size = 0.8) +
      
      # Add interpolated diameter line (red)
      geom_line(aes(x = segments$height, y = segments$dia.interp, color = "Interp. Diameter"), linewidth = 0.9) +
      
      # Add vertical dashed lines for 1.3 m and 7 m heights (blue)
      geom_vline(aes(xintercept = 1.3, color = "Height = 1.3 m"), linetype = "dashed", linewidth = 0.8) +
      geom_vline(aes(xintercept = 7, color = "Height = 7 m"), linetype = "dashed", linewidth = 0.8) +
      
      # Customize labels
      ylab("Stem Diameter [m]") +
      xlab("Height [m]") +
      ggtitle(paste("Tree", treeid, "IRLS Circle Fit, spline TLS vol =", round(seg.vol.smooth, 3), "m^3")) +
      
      # Customize x-axis limits
      xlim(0, max(tree_cloud@data[["Z"]])) +
      
      # Minimal theme
      theme_minimal() +
      
      # Customize legend
      scale_color_manual(
        name = "Legend",  # Title for the legend
        values = c(
          "Height = 1.3 m" = "blue",
          "Height = 7 m" = "darkgreen",
          "Interp. Diameter" = "red",
          "Segment Diameter" = "black"
        )
      ) +
      
      # Adjust legend appearance
      theme(
        legend.position = c(0.8, 0.8),  # Top-right corner inside the plot area
        legend.background = element_rect(fill = "white", color = "white"),
        plot.title = element_text(size = 11)
      )
    
    plotname <- file.path(outpath_plots, paste0("stem_curve_", basename(file), ".png"))
    ggsave(plotname, stem_plot, bg = "white")
    
    # Save results
    vol.results.table <- data.frame(treeid = treeid,
                                    vol_stem_raw = seg.vol.raw, 
                                    vol_stem_smooth = seg.vol.smooth, 
                                    DBH_splineTreeLS = DBH,
                                    D7_splineTreeLS = D7,
                                    pointFile = basename(file))
    csvname <- file.path(outpath, paste0("stemVolumes_circle_irls_", basename(file), ".csv"))
    write.csv(vol.results.table, csvname, row.names = FALSE)
  }else{
    vol.results.table <- data.frame(treeid = treeid,
                                    vol_stem_raw = NA, 
                                    vol_stem_smooth = NA, 
                                    DBH_splineTreeLS = NA,
                                    D7_splineTreeLS = NA,
                                    pointFile = basename(file))
    csvname <- file.path(outpath, paste0("stemVolumes_circle_irls_", basename(file), ".csv"))
    write.csv(vol.results.table, csvname, row.names = FALSE)
  }
  
  message("Results saved to: ", csvname)
  return(vol.results.table)
}

# file <- "D:/STSMGent/data/renamed_pointclouds_filtered/TLS_PL_REMB_183_filtered.laz"
# results <- calculate_stem_volume(file = file, outpath = "D:/STSMGent/data/renamed_pointclouds_metrics", slice_length = 0.2)

