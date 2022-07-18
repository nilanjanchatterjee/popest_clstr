
pkgs <- c("factoextra",  "NbClust", "WeightedCluster")
install.packages("WeightedCluster", dependencies = T)
library(factoextra)
library(NbClust)
library(WeightedCluster)

#####################################################################################3
setwd("D:/Data_hub/TATR_CT/Nianjan_WD analysis")
###2014 data 
wddat <-read.csv("WD_locations_data.csv", header=TRUE)
head(wddat)
long <-wddat$northing
lat <-wddat$easting
x<- cbind(long, lat)
dst <-dist(x)
clstr <- hclust(dst, method = "complete", members = wddat$indv)

aggMvad <- wcAggregateCases(wddat[, 9:10])
uniqueMvad <- wddat[aggMvad$aggIndex, 9:10]
mvad.seq <- seqdef(uniqueMvad, weights = wddat$indv)
## Computing Hamming distance between sequence
diss <- seqdist(mvad.seq, method = "HAM")

averageTree <- as.seqtree(clstr, seqdata = mvad.seq, diss = dst,ncluster = 6)
seqtreedisplay(averageTree, type = "d", border = NA)
avgClustQual <- as.clustrange(clstr,  dst, weights = wddat$indv, ncluster = 10)
plot(avgClustQual)
write.csv(avgClustQual$stats, "2014_clustering_stat.csv")

### available methods are single,complete,avrage,wald,mcquitty ###
plot(clstr, main="Dendogram of wild-dog location",  xlab="location id",hang=-1)
rect.hclust(clstr, k=7, border="red")
groups<-cutree(clstr, k=7)
grpinfo <- cbind(wddat,groups)
write.csv(grpinfo, "D:/wdclust_2014_7grp.csv")

fviz_nbclust(x, FUN = hcut, method = "wss")
NbClust(data = x, distance = "euclidean",  min.nc = 2, max.nc = 15, method = "complete")


##################################################################################
##### 2016 analysis

wd16 <-read.csv("WD_loc_2016.csv")
head(wd16)
wd16$indv[wd16$indv > 8] <-8
loc<- cbind(wd16$Northing, wd16$Easting)
dst1 <-dist(loc)
clstr1 <- hclust(dst1, method = "complete", members = wd16$num)
plot(clstr1, hang=-1)
rect.hclust(clstr1, k=9, border="red")
groups<-cutree(clstr1, k=10)
grpinfo1 <- cbind(wd16,groups)
write.csv(grpinfo1, "wdclust_2016_10grp1.csv")


agg_wg <- wcAggregateCases(wd16[, 7:8])
uniqueMvad <- wddat[agg_wg$aggIndex, 7:8]
mvad.seq1 <- seqdef(uniqueMvad, weights = wd16$num)
## Computing Hamming distance between sequence
diss <- seqdist(mvad.seq, method = "HAM")

averageTree1 <- as.seqtree(clstr1, seqdata = mvad.seq1, diss = dst1,ncluster = 11)
seqtreedisplay(averageTree1, type = "d", border = NA)
avgClustQual1 <- as.clustrange(clstr1,  dst1, weights = wd16$indv, ncluster = 11)
avgClustQual1$stats
write.csv(avgClustQual1$stats, "2016_clustering_stat1.csv")

### Read the indicator values and plot
cls_16 <-read.csv("2016_clustering_stat1.csv")
head(cls_16)
jpeg("2014_clustrng1.jpeg", width = 12, height = 8, units = "in", res = 300)
plot(2:11, cls_16$HG, pch=19, ylim=c(0,1), type= "l", col= "blue", xlab="Number of clusters", ylab= "Indicator value", lwd=2)
par(new=T)
plot(2:11, cls_16$R2, pch=19, ylim=c(0,1), type= "l", col= "purple",xlab="Number of clusters", ylab= "Indicator value", lwd=2)
par(new=T)
plot(2:11, cls_16$HC, pch=19, ylim=c(0,1), type= "l", col= "green",xlab="Number of clusters", ylab= "Indicator value", lwd=2)
par(new=T)
plot(2:11, cls_16$ASWw, pch=19, ylim=c(0,1), type= "l", col= "brown",xlab="Number of clusters", ylab= "Indicator value", lwd=2)
par(new=T)
plot(2:11, cls_16$PBC, pch=19, ylim=c(0,1), type= "l", col= "black",xlab="Number of clusters", ylab= "Indicator value", lwd=2)
legend("bottomright", lwd=2,legend = c("Hubert's Gamma", "Pseudo R2", "Hubert's C", "Average Silhouette Width",
                                 "Point Biserial corrln"), bty="n", col= c("blue", "purple", "green", "brown", "black"))
dev.off()

######################################################################################################
##### Point pattern analysis
library(raster)
library(spatstat)
library(RColorBrewer)
library(plyr)
library(maptools)

ind <-shapefile("D:/Work/Ungulate Point process/India_Boundary/India_Boundary.shp")
tatr <-shapefile('D:/Work/Conferences/ISEC/ISEC_2018/ISEC workshop/DSM/new analysis17/tatr.shp')
wd_cls_16 <-read.csv("wdclust_2016_9grp1.csv", header =T)
wd_18 <-read.csv("dhole_idw_2018.csv")
head(wd_cls_16)
tatr
tatr <-spTransform(tatr, crs(ind))

wd_loc16 <-SpatialPoints(wd16[,8:9], proj4string = crs(tatr))
wd_loc_16 <-spTransform(wd_loc16, crs(tatr))
wd_loc14 <- SpatialPoints(wddat[,2:3], proj4string = crs(tatr)) 
wd_loc18 <-SpatialPoints(wd_18[,2:3], proj4string = crs(tatr))
  
plot(ind)
plot(trtl, add=T , pch=19, col="purple")

xx <-tatr@polygons[[1]]@Polygons[[4]]@coords[,1]
yy <-tatr@polygons[[1]]@Polygons[[4]]@coords[,2]
Z <- owin(poly=list(x=rev(xx), y=rev(yy)))

loc_16 <-ppp(wd_cls_16[,8], wd_cls_16[,9], window = Z)
plot(loc_16, add=T)
summary(loc_16)

loc_14 <-ppp(wddat[,2], wddat[,3], window = Z)
plot(loc_14, add=T)
summary(loc_14)

loc_18 <-ppp(wd_18[,2], wd_18[,3], window = Z)
plot(loc_18, add=T)
summary(loc_18)
plot(density(loc_18, 0.01, weights= wd_18$Dhole))

jpeg("WD_density_14_16.jpeg", width = 12, height = 8, units = "in", res = 300)
par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(density(loc_16, 0.01, weights= wd16$num), col=brewer.pal(9, "YlGnBu"), scalekernel=T, kernel= "gaussian", 
     main= "Density plot of 2016 with wild-dog captures")
plot(tatr, add=T)

plot(density(loc_14, 0.01, weights= wddat$indv), col=brewer.pal(9, "YlGnBu"), scalekernel=T, kernel= "gaussian", 
     main= "Density plot of 2014 with wild-dog captures")
plot(tatr, add=T)
dev.off()

 #plot(density(loc, 0.01, weights= wd16$indv), col=brewer.pal(9, "YlGnBu"),      scalekernel=T, kernel= "epanechnikov")
#plot(density(loc, 0.01, weights= wd16$indv), col=brewer.pal(9, "YlGnBu"),      scalekernel=T, kernel= "quadratic")
#plot(density(loc, 0.01, weights= wd16$indv), col=brewer.pal(9, "YlGnBu"),      scalekernel=T, kernel= "disc")

plot(den_16, col=brewer.pal(9,"PRGn"))

#### 2016 analysis and plot
den_16 <-density(loc_16, 0.01, weights= wd16$num)
den_16_sp <-as.SpatialGridDataFrame.im(den_16)
den_16_r <-as(den_16_sp, "RasterLayer")
writeRaster(den_16_r, "WD_den_16.tif", format= "GTiff", overwrite= T)

jpeg("WD_den_16_9grp.jpeg", width = 10, height = 8, units = "in", res = 300)
par(mar=c(5,5,1,2))
plot(den_16_r, xlab = "Longitude", ylab = "Latitude", ylim=c(19.8,20.6), cex.lab=2)
plot(tatr, add=T)
for(i in 1:9){
  grp1 <-subset(wd_cls_18, wd_cls_18$groups==i)
  dat1 <- cbind(grp1$easting, grp1$northing)
  dat <-cbind(grp1$northing, grp1$easting)
  cen <-apply(grp1[,3:4], 2, mean)
  euc.dist <- function(dat, cen) sqrt(sum((dat - cen) ^ 2))
  dist <- NULL
  for(j in 1:nrow(dat)) dist[j] <- euc.dist(dat[j,],cen)
  
  dsigma[i] <-mean(dist)
  
  ch1 <-chull(dat1)
  coords1 <- dat1[c(ch1, ch1[1]), ]
  #lines(coords)
  sp_poly1 <- SpatialPolygons(list(Polygons(list(Polygon(coords1)), ID=1)), proj4string = crs(tatr))
  sp_poly_df <- SpatialPolygonsDataFrame(sp_poly1, data=data.frame(ID=1))
  plot(sp_poly_df, add=T, lwd=2, border= "black")
}
dev.off()


#### 2014  analysis and plots
den_14 <-density(loc_14, 0.01, weights= wddat$indv)
den_14_sp <-as.SpatialGridDataFrame.im(den_14)
den_14_r <-as(den_14_sp, "RasterLayer")
writeRaster(den_14_r, "WD_den_14.tif", format= "GTiff", overwrite= T)

jpeg("WD_den_14.jpeg", width = 10, height = 8, units = "in", res = 300)
par(mar=c(5,5,1,2))
plot(den_14_r, xlab = "Longitude", ylab = "Latitude", ylim=c(19.8,20.6), cex.lab=2)
plot(tatr, add=T)
dsigma <- NULL
for(i in 1:7){
  grp1 <-subset(wddat, wddat$groups_7==i)
  dat1 <- cbind(grp1$DD.long, grp1$DD.lat)
  dat <-cbind(grp1$northing, grp1$easting)
  cen <-apply(grp1[,9:10], 2, mean)
  euc.dist <- function(dat, cen) sqrt(sum((dat - cen) ^ 2))
  dist <- NULL
  for(j in 1:nrow(dat)) dist[j] <- euc.dist(dat[j,],cen)
  
  dsigma[i] <-mean(dist)
  ch1 <-chull(dat1)
  coords1 <- dat1[c(ch1, ch1[1]), ]
  #lines(coords)
  sp_poly1 <- SpatialPolygons(list(Polygons(list(Polygon(coords1)), ID=1)), proj4string = crs(tatr))
  sp_poly_df <- SpatialPolygonsDataFrame(sp_poly1, data=data.frame(ID=1))
  plot(sp_poly_df, add=T, lwd=2, border= "black")
}
dev.off()

###########################################################################
#### Plotting the dhole specific cameras

dhl_cam <-read.csv("csv_dhole_cam.csv")
head(dhl_cam)
indv_cpt <-read.csv("indv_captures.csv")
head(indv_cpt)
indv <-unique(indv_cpt$indv)

camloc_18 <-read.csv("D:/Data_hub/TATR_CT/2018_ctpoints.csv")
head(camloc_18)
dhl_loc_18 <-read.csv("dhole_loc_18.csv")

cam_loc_utm <-SpatialPoints(dhl_cam[,3:2], proj4string = crs("+proj=utm +zone=44N +datum=WGS84")) 
cam_loc_ll <- spTransform(cam_loc_utm, crs(tatr))

plot(tatr)
plot(cam_loc, add=T, pch=19, col= "brown") ###Cameras deployed for dhole

dhl_loc_utm <-SpatialPoints(dhl_loc_18[,3:2], proj4string = CRS("+proj=utm +zone=44N +datum=WGS84"))
dhl_loc_ll<- spTransform(dhl_loc_utm, crs(tatr))
plot(dhl_loc_ll, add=T, cex=2)

for(i in 1:3){
  grp1 <-subset(indv_cpt, indv_cpt$indv== indv[i])
  dat1 <- cbind(grp1$long, grp1$lat)
  ch1 <-chull(dat1)
  coords1 <- dat1[c(ch1, ch1[1]), ]
  #lines(coords)
  sp_poly1 <- SpatialPolygons(list(Polygons(list(Polygon(coords1)), ID=1)), proj4string = crs(tatr))
  sp_poly_df <- SpatialPolygonsDataFrame(sp_poly1, data=data.frame(ID=1))
  plot(sp_poly_df, add=T, lwd=2, border= "blue")
}
indv[4]


#######################################################################
######### Function to estimate the pack size with the confidence interval


wd_clstr1 <-function(x){
  #x <-unique(x)
  max_grp <-max(x)
  detp_avg <-mean((max_grp-x)/max_grp)
  detp_sd <-sd((max_grp-x)/max_grp)
  est_grp <-max_grp/(1- detp_avg)
  lcl_grp <-max_grp/(1- (detp_avg-detp_sd))
  ucl_grp <-max_grp/(1- (detp_avg+detp_sd))
  return(c(est_grp, lcl_grp, ucl_grp, detp_sd, detp_avg))
}

############# 2014 groups
grp1 <-c(1,1,3,2,1)
grp2 <-c(1,2,1,1,1,2,3,2,2,2,2)
grp3 <-c(1,1,3,3,2,1,2,5,3,2,2,3,1,2)
grp4 <-c(2,2,2,4,5,5,6,2,2)
grp5 <-c(2,2,2,1,2,2,1,3,2,3,1,1,1,1,6)
grp6 <-c(6,2,1,1,1,2,1)
grp7 <-c(1,1,2,2,1,1,1,2,3,1,5,4)

#############  2016 groups
grp1 <-c(4,1,1,1,1,1,5,4,1,1,1,1,1,1,1,1)
grp9 <-c(5,5,7,1,1,4,1,1,2,1,2,2,1,2)
grp8 <-c(2,1,1,3,1,1,1,2,2,1,1,1,2,1,1,2,3,1,1)
grp7 <-c(1,2,5,1,3,1,1,5,4,1)
grp6 <-c(1,1,1,4)
grp5 <-c(1,1,2,2,2,2,1,1,1)
grp4 <-c(2,2,1,1,5,3,1,2,1,2,3,1,2,4,3,1)
grp3 <-c(6,1,1,2,3,1,1,1,3,4,1,1,5,5,6,3,1,1,2,1,1,5,4,5,5,1)
grp2 <-c(2,2,5,1,1,4,5)


wdsum <- matrix(data = NA, nrow = 9, ncol=5)
wdsum[7, ] <-wd_clstr1(grp7)

wdsum
colSums(wdsum[1:7,], na.rm = T)

wdsum1 <-rbind(wdsum[1:7,], colSums(wdsum[1:7,], na.rm=T))
wdsum1
write.csv(wdsum[,1:5],"WD_calc_14_edtd.csv")

############## 2018 groups
grp1 <-c(2,1,1,2,2,2,1,1,3,3,2,2,1,1,1,2,3,3,1,1)
grp9 <-c(1,1,1,1,1,1,1,1)
grp8 <-c(5,3,1,2,1,6,1,5,2,7,4,4,1,2,3,6,6,1,4,2,6,1,5,6,1,6,2,1,3,6,3,6,6,1,4,1,1,1,3,2,2,2,1,1,1,1,1,2,1,2,2,1)
grp7 <-c(1,2,3,1,1,2,3,1,2,3,1,4,2,1,1,1,1,5,1,5,1,5,2,1,5,2,3,2,1,4,1,2,1,3,3,4,5,3,3,3,3,1)
grp6 <-c(1,1,2,3,1,3,1,3,1,1,5,6,1,2,4,1)
grp5 <-c(5,1,1,2,4,2,1,4,1,5,1,5,2,2,1,1,1,4,4,4,1,1,2)
grp4 <-c(1,1,1,2,3,1,1,1)
grp3 <-c(4,4,5,2,1,1,1,1,1,1,2,2,1,1,1,2,2,2,1,1,2,1,3,1,2,2,1,2,2,2,1,1,1,2,2,2,2)
grp2 <-c(1,1,2,1,1,1,2,1,1,1,1,4,1,2,2,1,2,2,1,4,5,5,1,6,6,1,1,1,5,5)

wdsum <- matrix(data = NA, nrow = 9, ncol=3)
wdsum[8, ] <-wd_clstr(grp8)

wdsum
colSums(wdsum) #### Estimated population with confidence interval

wdsum <-rbind(wdsum, colSums(wdsum))
wdsum
write.csv(wdsum,"WD_calc_18.csv")

##########################################################################################################
###### 2018 analysis with marked individuals

library(camtrapR)
ct18 <-read.csv("D:/Data_hub/TATR_CT/2018_ctpoints.csv")
head(ct18)
colnames(ct18)<- c("station","Northing"  , "Easting"  ,  "Start_date" ,"End_date"  , "BLOCK")
camop18 <-cameraOperation(ct18, stationCol = "station", 
                          setupCol = "Start_date", retrievalCol = "End_date",
                          dateFormat = "%d-%m-%Y",
                          writecsv = T)

rectab18 <-read.csv("D:/Data_hub/TATR_CT/record_table_1min_deltaT_2018-12-22.csv")
head(rectab18)

wd_det_18 <-detectionHistory(recordTable = rectab18,
                             species = "dhole",
                             camOp = camop18,
                             stationCol = "station",
                             recordDateTimeCol = "Date",
                             recordDateTimeFormat = "%Y-%m-%d",
                             occasionLength=1,
                             speciesCol = "Species",
                             includeEffort= T, scaleEffort = F, 
                             timeZone = "Asia/Calcutta",
                             day1 = "station")
head(wd_det_18)

write.csv(wd_det_18$detection_history, "Wd_det_18.csv")

#######################################################################################
################  DENSITY COMPARISON
library(ggplot2)
library(tidyverse)

pop_est <-read.csv("D:/Papers/Dhole Population/Dhole_est_com.csv")
head(pop_est)

ggplot(pop_est, aes(x=Method, y=(Estimate*100/ETA), colour=as.factor(year))) + 
  geom_errorbar(aes(ymin=(lcl*100/ETA), ymax=(ucl*100/ETA), colour=as.factor(year)), 
                width=.2, size=1.5, position = position_dodge(0.35)) +
  geom_point(position = position_dodge(0.35), aes(size=1.5))+
  theme_bw()+
  labs(y="Density Estimate per 100 sq. km", colour= "Year", size="")+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=16), 
        legend.title = element_text(size=18),
        legend.text = element_text(size=14))
ggsave("Pop_est.jpeg", width = 12, height = 8, units = "in", dpi = 300)

ggplot(pop_est, aes(x=Method, y=(Estimate), colour=as.factor(year))) + 
  geom_errorbar(aes(ymin=(lcl), ymax=(ucl), colour=as.factor(year)), width=.2, size=1.5, position = position_dodge(0.2)) +
  geom_point(position = position_dodge(0.2), aes(size=1.5))+
  theme_bw()+
  labs(y="Population Estimate", colour= "Year", size="")+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size=14), 
        legend.title = element_text(size=16),
        legend.text = element_text(size=12))

pop_est %>% 
  drop_na(sigma) %>%
  ggplot(aes(x=Method, y=sigma, colour=as.factor(year)), na.rm=T) + 
  geom_errorbar(aes(ymin=(sig_lcl), ymax=(sig_ucl), colour=as.factor(year)), 
                width=.2, size=1.5, position = position_dodge(0.35)) +
  geom_point(position = position_dodge(0.35), aes(size=1.5))+
  theme_bw()+
  labs(y="Sigma Estimate (in meters)", colour= "Year", size="")+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=16), 
        legend.title = element_text(size=18),
        legend.text = element_text(size=14))

ggsave("Sig_est_final.jpeg", width = 12, height = 8, units = "in", dpi = 300)
