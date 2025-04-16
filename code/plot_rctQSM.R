
# plot rayextract QSMs

#library(png)
#library(stringr)
library(rgl)
library(Rvcg)
library(magick)


ply_grapher <- function(file_path, color, out_path) {
  ##function to create 2d prints of mesh ply files
  #file_path input is a directory containing meshes (.ply) and point cloud (.ply)
  #Creates a directory out_path 
  #color is a hexadecimal input that will monocolor the png print
  
  ##prepare input files
  #set wd to working figure folder
  if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
  #list obj files in input directory
  files <- list.files(path = file_path, pattern = "raycloud_trees_mesh.ply", full.names = TRUE, recursive = FALSE)
  
  #loop through all .obj files in data directory and generate png prints, saved in input directory
  for (i in 1:length(files)){ #length(files)
    #extract tree sample name (strip off file extension)
    str_temp <- sub(".*/","", files[i])
    str_temp <- str_extract(str_temp, "[^.]+")
    #read first ply file and store in mesh
    mesh <- vcgPlyRead(files[[i]])
    
    # Look for corresponding point cloud file
    pc_file <- list.files(path = file_path, pattern = paste0(sub("_raycloud_trees_mesh", "", str_temp), "_raycloud_segmented.ply"), full.names = TRUE)
    
    if (length(pc_file) > 0) {
      pc <- vcgPlyRead(pc_file)
      #print("reading point cloud")
    } else {
      pc <- NULL
    }
    #display model with shape4d
    open3d()
    shade3d(mesh, color = color, alpha = 0.5)
    if (!is.null(pc)) {
      points3d(pc$vb[1,], pc$vb[2,], pc$vb[3,], color = "black", size = 1)
    }
    par3d("windowRect" = c(500,0,2000,2000))
    view3d(theta = 0, phi = -90, zoom = 0.70)
    img1 <- tempfile(fileext = ".png")
    rgl.snapshot(img1)
    
    clear3d()
    
    shade3d(mesh, color = color)
    par3d("windowRect" = c(500,0,2000,2000))
    view3d(theta = 0, phi = -90, zoom = 0.70)
    img2 <- tempfile(fileext = ".png")
    rgl.snapshot(img2)
    
    clear3d()
    close3d()
    
    img_combined <- image_append(c(image_read(img1), image_read(img2)))
    image_write(img_combined, path = file.path(out_path, paste0(str_temp, ".png")))
    #clear device for loop
    clear3d()
    close3d()
  }
}



# file_path <- "F:/CRS-Allometries-Datasets/forQSM/trees"
# color <- "#f6b251"
# out_path <- "F:/CRS-Allometries-Datasets/forQSM/qsm_images"
# 
# ply_grapher(file_path, color, out_path)
