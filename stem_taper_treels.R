# stem volume

library(dplyr)
#library(tidyr)
library(ggplot2)
library(cowplot)
library(lidR)
library(TreeLS)
library(pracma)
library(readxl)
library(stringr)


# set up table for storing results
vol.results.table = data.frame(VolTLS= numeric(0), volTLSsmooth= numeric(0),  pointFile = character(0))

file <- paste0(pc_path[t], filename[i])
# read tree cloud and convert to LAS object
# baum <- fread(file[i])
# baum <- baum[,2:4]
# baum <- setTLS(baum)

# read las tree cloud
#baum <- readTLS("D:/TreeVolumeData/Waldlabor/Pointclouds/T19_Bu19_WL_cleaned.las")
disp(paste("processing:", filename[i]))
baum <- readTLS(file)

# normalize, as in subtract minimum z
minZ = min(baum@data[["Z"]]);
baum@data[["Z"]] = baum@data[["Z"]] - minZ;

# classify stem points
baum.stem = stemPoints(baum, stm.hough())

# store stem points separately 
stem = filter_poi(baum.stem, Stem == TRUE)

#plot(stem)


# define length of stem slices (in meters)
s = 0.2
# set up sequence of heights for stem slices, starting at lowest stem point
height = seq(from=min(stem@data[["Z"]]), by=s, to= floor(max(stem@data[["Z"]])))

#slices <- data.frame(rho= NA, theta= NA, phi = NA, alpha =NA, Radius= NA, Error= NA, PX= NA, PY=NA, PZ=NA)
slices <- data.frame(X= NA, Y= NA, Radius= NA, Error= NA)


# slice stem point cloud and fit a circle or cylinder if it contains enough points
for(j in 1:length(height)) {

  segment = filter_poi(stem, Z > height[j]-s*0.75 & Z < height[j]+s*0.75) 

    if(length(segment@data[["Z"]])< 50) {
      #slice = data.frame(rho= NA, theta= NA, phi = NA, alpha =NA, Radius= NA, Error= NA, PX= NA, PY=NA, PZ=NA) # set diameter to NA if the slice contains less than 10 points
      slice <- data.frame(X= NA, Y= NA, Radius= NA, Error= NA)
      } else {
      slice = shapeFit(segment, shape='circle', algorithm='irls')
      #slice = shapeFit(segment, shape='cylinder', algorithm='irls')
      }

  slices = rbind(slices, slice)
}

segments = cbind(slices, height = height[1:nrow(slices)])

# the volume of each stem section was then estimated according to Huber's formula
# and the total stem volume was then estimated as the sum over all sections
# cross-sectional area at the middle of a stem section multiplied with the length, i.e. 10 cm

segments$dia = segments$Radius*2


# rough total stem volume
segments$vol.raw = pi*(segments$dia/2)^2*0.1
seg.vol.raw = sum(segments$vol.raw, na.rm=TRUE)

#segments$dia[segments$Error > 0.3] = NA # set segments with very bad fit to NA-->Loch
segments$dia[segments$Radius > segments$Radius[1]+0.1] = NA # set segments with Radius bigger than 10cm than lowest to NA --> Loch

segments = segments[!is.na(segments$height),]

# start segments with first cylinder
NonNAindex <- which(!is.na(segments$dia))
firstNonNA <- min(NonNAindex)
segments <- segments[firstNonNA:length(segments$height),]

# 
for(k in 4:length(segments$height)) { # absolute difference zu mean of previous 3 diameters must not be more than 10 cm, (else set to Na) else set to 0.01 smaller than previous dia

  if(is.na(segments$dia[k])) {
    segments$dia[k] = NA
  } else if(abs(mean(c(segments$dia[k-3],segments$dia[k-2],segments$dia[k-1]), na.rm=TRUE) - segments$dia[k]) > 0.1 || is.na(mean(c(segments$dia[k-3],segments$dia[k-2],segments$dia[k-1]), na.rm=TRUE)) ) {
    segments$dia[k] = segments$dia[k-1] - 0.01
    #segments$dia[k] = NA
  #} else if(segments$dia[k] - mean(c(segments$dia[k-3],segments$dia[k-2],segments$dia[k-1]), na.rm=TRUE) > segments$dia[k]/10 || is.na(mean(c(segments$dia[k-3],segments$dia[k-2],segments$dia[k-1]), na.rm=TRUE)) ) {
  #  segments$dia[k] = segments$dia[k-1]    
  } else {segments$dia[k] = segments$dia[k]
    }
}

plot(segments$height, segments$dia)


# interpolate with spline to fill holes

# Fit a cubic smoothing spline

segments.comp = segments[!is.na(segments$dia),]
sm = stats::smooth.spline(segments.comp$height, segments.comp$dia, spar=0.7)

dia.interp = predict(sm, segments$height)

segments$dia.interp = dia.interp$y

segments = dplyr::filter(segments, dia.interp >= 0.07) # remove (tree top!) segments with diameter smaller than 7cm (nur Derbholz)

segments$vol = pi*(segments$dia.interp/2)^2*0.1


#plot(segments$height, segments$Radius*2)
#lines(sm)

#plot(segments$height, segments$dia_smooth)

# total stem volume from smoothed diameters
seg.vol.smooth = sum(pi*(segments$dia.interp/2)^2*0.1)


stem_plot <- ggplot()+
  geom_point(aes(segments$height, segments$dia), size=0.8, colour= 'black')+
  geom_line(aes(segments$height, segments$dia.interp), linewidth=0.9, colour= 'red')+
  ylab("Stem Diameter [m]")+
  xlab("Height [m]")+
  ggtitle(paste("irls circle fit,","TLS. vol.:", round(seg.vol.smooth,3)))+
  xlim(0, max(baum@data[["Z"]]))+
  #ylim(0, max(segments$dia, na.omit=TRUE)+0.2)+
  theme_minimal() +
  theme(plot.title = element_text(size=11))
stem_plot

plotname <- paste0(oupath, treeid, '_stem_curve.png') 
ggsave(plotname, stem_plot, bg="white")

vol.results = data.frame(VolTLS.raw= seg.vol.raw, volTLS.smooth= seg.vol.smooth, pointFile = file[i])

vol.results.table = rbind(vol.results.table, vol.results)


write.csv(vol.results.table, paste0(oupath, treeid, "_stemVolumes_circle_irls_s02.csv"))

rm(vol.results.table)





# Ahornholz: Rohdichte Mittelwert	Bergahorn 623 kg/m3. 
# gesamtgewicht derbholz: 3000 kg
#1790/700



# stem_plot <- ggplot()+
#   geom_point(aes(segments$dia, segments$height),shape=16, size=0.8, colour= 'black')+
#   #geom_line(aes(segments$dia.interp, segments$height),size=0.8, colour= 'red')+
#   geom_point(aes(baum.field.df$dA/100, baum.field.df$h),shape=17, size=2, colour='cyan4')+
#   geom_point(aes(baum.field.df$dB/100, baum.field.df$h),shape=17, size=2, colour='blue1')+
#   #geom_line(aes(ref.smooth$x, ref.smooth$y/100), size=0.8, colour= 'darkgreen')+
#   xlab("Stem Diameter [m]")+
#   ylab("Height [m]")+
#   #ggtitle(paste("irls circle fit,","TLS. vol.:", round(seg.vol.smooth,3),"ref. vol.:", round(ref.vol.interp, 3)))+
#   ylim(0, max(baum@data[["Z"]]))+
#   xlim(0,0.6)+
#   theme_minimal() +
#   theme(plot.title = element_text(size=11))
# stem_plot
# 
# 
# ggplot()+
#   geom_point(aes(segments$height, segments$dia), size=1, colour= 'black')+
#  # geom_point(aes(segments$height, segments$dia.interp), size=0.8, colour= 'red')+
#   geom_point(aes(baum.field.df$h, baum.field.df$dA/100),shape=17,colour='lightblue')+
#   geom_point(aes(baum.field.df$h, baum.field.df$dB/100),shape=17, colour='blue')+
#  # geom_point(aes(ref.smooth$x, ref.smooth$y/100), size=0.8, colour= 'darkgreen')+
#   ylab("Stem Diameter [m]")+
#   xlab("Height [m]")+
#  # ggtitle(paste("irls cylinder fit,","TLS. vol.:", round(seg.vol.smooth,3),"ref. vol.:", round(ref.vol.interp, 3)))+
#   xlim(0, max(baum@data[["Z"]]))+
#   ylim(0,0.6)+
#   theme_minimal() +
#   theme(plot.title = element_text(size=11))



#-------------------------------------------------------------------------
# extract stem measures
# seg = stemSegmentation(stem, sgt.ransac.circle(n = 20))
# #add_stemSegments(x, seg, color='white', fast=T)
# 
# seg.irls = stemSegmentation(stem, sgt.irls.cylinder())
# #add_stemSegments(x, seg.irls, color='white', fast=T)
# 
# seg.bf = stemSegmentation(stem, sgt.bf.cylinder())
# 
# library(DescTools)
# 
# segment = filter_poi(stem, Z > 1.3 & Z < 4)
# segment = filter_poi(stem, Z > height[3]-0.05 & Z < height[3]+0.05)
# pars = shapeFit(segment, shape='circle', algorithm='irls')
# segment@data %$% plot(Y ~ X, pch=20, asp=1)
# pars %$% points(X,Y,col='red', pch=3, cex=2)
# pars %$% lines(c(X,X+Radius),c(Y,Y), col='red',lwd=2,lty=2)
# pars %$% DrawCircle(X, Y, Radius)
# 
# #-------------------------------------------------------------------------
# ggplot()+
#   geom_point(aes(segments$height, segments$dia), size=0.8)+
#  # geom_point(aes(segments$height, segments$dia.interp), size=0.8, colour= 'red')+
#   geom_smooth(aes(segments$height, segments$dia), method=loess)+
#  # geom_point(aes(baum.field.df$h, baum.field.df$dA/100),shape=17,colour='green')+
#  # geom_point(aes(baum.field.df$h, baum.field.df$dB/100),shape=17, colour='blue')+
#  # geom_point(aes(ref.smooth$x, ref.smooth$y/100), size=0.8, colour= 'darkgreen')+
#   ylab("Stem Diameter [m]")+
#   xlab("Height [m]")+
#   ggtitle(paste("Stem curve,","stem volume =", round(seg.vol.smooth,3), "m^3"))+
#   xlim(0, max(baum@data[["Z"]]))+
#   ylim(0,0.6)+
#   theme_minimal() +
#   theme(plot.title = element_text(size=11))
# #-------------------------------------------------------------------------
# 
# ggplot(vol.results.table, aes(x=volTLS.smooth, y= volRef.smooth, label=BALABEL)) +
#   geom_point() +
#   geom_text(hjust = 0, nudge_x = 0.05) +
#   geom_smooth(method=lm, se=FALSE) +
#   geom_abline(slope=1, intercept= 0) +
#   stat_regline_equation( aes(label = ..rr.label..)) +
#   xlab("TLS stem volume [m^3]")+
#   ylab("reference stem volume [m^3]")+
#   theme_bw() + scale_x_continuous(expand = c(0, 0), limits = c(0, 8)) + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0,8))


