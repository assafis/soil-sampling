#list.of.packages <- c("sp", "geoR","gstat","cluster","colorspace","dplyr","ggplot2","raster","e1071","geosphere","rgdal","rgeos","RColorBrewer")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

library(sp)
library(geoR)
library(gstat)
library(cluster)
library(colorspace)
library(dplyr)
library(ggplot2)
library(raster)   
library(e1071)
library(geosphere)
library(rgdal)
library(rgeos)
library(reshape2)
library(RColorBrewer)

source('spatial_helper_functions.R')

#####################################


# Read csv file to data-frame
csv.v1 = read.csv("/data/em38neveyaar/Neve Ver Pt1.csv",stringsAsFactors=FALSE)
csv.v2 = read.csv("/data/em38neveyaar/Neve Ver Pt2.csv",stringsAsFactors=FALSE)

csv.h1 = read.csv("/data/em38neveyaar/Neve 29-11 Hor Pt1.csv",stringsAsFactors=FALSE)
csv.h2 = read.csv("/data/em38neveyaar/Neve 29-11 Hor Pt2.csv",stringsAsFactors=FALSE)
csv.h3 = read.csv("/data/em38neveyaar/Neve 29-11 Hor Pt3.csv",stringsAsFactors=FALSE)
csv.h4 = read.csv("/data/em38neveyaar/Neve 29-11 Hor Pt4.csv",stringsAsFactors=FALSE)

df.v1 = as.data.frame(csv.v1)[1:4]
df.v2 = as.data.frame(csv.v2)[1:4]

df.h1 = as.data.frame(csv.h1)[1:4]
df.h2 = as.data.frame(csv.h2)[1:4]
df.h3 = as.data.frame(csv.h3)[1:4]
df.h4 = as.data.frame(csv.h4)[1:4]


colnames(df.v1)<-colnames(df.v2)<-c("x","y","ECaV","MSaV")
colnames(df.h1)<-colnames(df.h2)<-colnames(df.h3)<-colnames(df.h4)<-c("x","y","ECaH","MSaH")

# remove comments
df.v1<-df.v1[!grepl("//", df.v1$x),]
df.v2<-df.v2[!grepl("//", df.v2$x),]
df.h1<-df.h1[!grepl("//", df.h1$x),]
df.h2<-df.h2[!grepl("//", df.h2$x),]
df.h3<-df.h3[!grepl("//", df.h3$x),]
df.h4<-df.h4[!grepl("//", df.h4$x),]


colnames(df.v1)<-colnames(df.v2)<-c("x","y","ECaV","MSaV")
colnames(df.h1)<-colnames(df.h2)<-colnames(df.h3)<-colnames(df.h4)<-c("x","y","ECaH","MSaH")

# remove comments
df.v1<-df.v1[!grepl("//", df.v1$x),]
df.v2<-df.v2[!grepl("//", df.v2$x),]
df.h1<-df.h1[!grepl("//", df.h1$x),]
df.h2<-df.h2[!grepl("//", df.h2$x),]
df.h3<-df.h3[!grepl("//", df.h3$x),]
df.h4<-df.h4[!grepl("//", df.h4$x),]



df.v<-rbind(df.v1,df.v2)
df.h<-rbind(df.h1,df.h2,df.h3,df.h4)

df.v <- as.data.frame(apply(df.v,2, as.numeric))
df.h <- as.data.frame(apply(df.h,2, as.numeric))

# fix GPS position: 3m behind, approx # of measurements
#aa=32476
#bb=5
#df.h[aa+bb,1]-df.h[aa,1]
#df.h[aa+bb,2]-df.h[aa,2]
df.v<-cbind(df.v[1:(nrow(df.v)-5),1:2],df.v[6:nrow(df.v),3:4])
df.h<-cbind(df.h[1:(nrow(df.h)-5),1:2],df.h[6:nrow(df.h),3:4])


# compact data
df.v<- meanNv(df.v, 20)
df.h<- meanNh(df.h, 20)


# explore eca
summary(df.v)
summary(df.h)


# Histograms - raw 
h.theme <- theme(plot.subtitle = element_text(hjust=1,size=11),panel.background = element_rect(NA),panel.grid = element_line(size=0.1,colour = "grey80"),panel.border = element_rect(size=0.1, fill = NA))

sk1<-round(skewness(df.v$ECaV),3)
sk2<-round(skewness(df.v$MSaV),3)
sk3<-round(skewness(df.h$ECaH),3)
sk4<-round(skewness(df.h$MSaH),3)
hist.ecv<-ggplot(df.v, aes(ECaV)) + geom_histogram(bins=100) + labs(subtitle=paste("skewness ",sk1)) + h.theme
hist.msv<-ggplot(df.v, aes(MSaV)) + geom_histogram(bins=100) + labs(subtitle=paste("skewness ",sk2)) + h.theme
hist.ech<-ggplot(df.h, aes(ECaH)) + geom_histogram(bins=100) + labs(subtitle=paste("skewness ",sk3)) + h.theme
hist.msh<-ggplot(df.h, aes(MSaH)) + geom_histogram(bins=100) + labs(subtitle=paste("skewness ",sk4)) + h.theme

multiplot(hist.ecv,hist.msv,hist.ech,hist.msh,cols=2)



# trim data
df.v <- df.v[df.v$ECaV<quantile(df.v$ECaV,0.99) & df.v$ECaV>quantile(df.v$ECaV,0.01),]
df.h <- df.h[df.h$ECaH<quantile(df.h$ECaH,0.99) & df.h$ECaH>quantile(df.h$ECaH,0.01),]

sk1<-round(skewness(df.v$ECaV),3)
sk2<-round(skewness(df.v$MSaV),3)
sk3<-round(skewness(df.h$ECaH),3)
sk4<-round(skewness(df.h$MSaH),3)
hist.ecv<-ggplot(df.v, aes(ECaV)) + geom_histogram(bins=100) + labs(subtitle=paste("Trimmed | skewness ",sk1)) + h.theme
hist.msv<-ggplot(df.v, aes(MSaV)) + geom_histogram(bins=100) + labs(subtitle=paste("Trimmed | skewness ",sk2)) + h.theme
hist.ech<-ggplot(df.h, aes(ECaH)) + geom_histogram(bins=100) + labs(subtitle=paste("Trimmed | skewness ",sk3)) + h.theme
hist.msh<-ggplot(df.h, aes(MSaH)) + geom_histogram(bins=100) + labs(subtitle=paste("Trimmed | skewness ",sk4)) + h.theme

multiplot(hist.ecv,hist.msv,hist.ech,hist.msh,cols=2)



# convert to log
df.v$ECaV <- log(df.v$ECaV)
df.v$MSaV <- log(df.v$MSaV)
df.h$ECaH <- log(df.h$ECaH)
df.h$MSaH <- log(df.h$MSaH)

sk1<-round(skewness(df.v$ECaV),3)
sk2<-round(skewness(df.v$MSaV),3)
sk3<-round(skewness(df.h$ECaH),3)
sk4<-round(skewness(df.h$MSaH),3)
hist.ecv<-ggplot(df.v, aes(ECaV)) + geom_histogram(bins=100) + labs(subtitle=paste("Trimmed, log | skewness ",sk1)) + h.theme
hist.msv<-ggplot(df.v, aes(MSaV)) + geom_histogram(bins=100) + labs(subtitle=paste("Trimmed, log | skewness ",sk2)) + h.theme
hist.ech<-ggplot(df.h, aes(ECaH)) + geom_histogram(bins=100) + labs(subtitle=paste("Trimmed, log | skewness ",sk3)) + h.theme
hist.msh<-ggplot(df.h, aes(MSaH)) + geom_histogram(bins=100) + labs(subtitle=paste("Trimmed, log | skewness ",sk4)) + h.theme

multiplot(hist.ecv,hist.msv,hist.ech,hist.msh,cols=2)


N=22
br=quantile(df.v$ECaV, probs = seq(0, 1, length.out = N + 1), na.rm = TRUE) 

hist(df.v$ECaV, breaks=br,freq=TRUE,col="grey70",xlab="ECa V")
hist(df.v$ECaV, bins=20,freq=TRUE,col="grey70",xlab="ECa V")


sp.v<-df.v
sp.h<-df.h
# Set spatial coordinates to create a Spatial object:
coordinates(sp.v) = ~x + y
coordinates(sp.h) = ~x + y



wgs84 <- "+proj=longlat + ellps=WGS84"
utm <- "+proj=utm +zone=36N +datum=WGS84"
proj4string(sp.v) <- CRS(utm)
proj4string(sp.h) <- CRS(utm)
#dfcv1.7 <-spTransform(dfcv1.7,CRS=utm)   # re-project WGS84 to UTM

# Define the grid extent (measurements + 5m)
x.range <- as.numeric(c(min(min(df.v$x)-5,min(df.h$x)-5), max(max(df.v$x)+5,max(df.h$x)+5)))
y.range <- as.numeric(c(min(min(df.v$y)-5,min(df.h$y)-5), max(max(df.v$y)+5,max(df.h$y)+5)))

# Expand points to grid
grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 1), y = seq(from = y.range[1], to = y.range[2], by = 1))  
grd10 <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 10), y = seq(from = y.range[1], to = y.range[2], by = 10))  

coordinates(grd) <- coordinates(grd10) <- ~x + y
gridded(grd) <- gridded(grd10) <- TRUE
proj4string(grd) <- proj4string(grd10) <- CRS(utm)

par(mfrow=c(1,2),bty = 'n',mar=c(2, 3, 2, 3))

#plot(grd, cex = 1.5, col = "grey")
plot(grd10, cex = 1.5, col = "grey",axes=T)
title(main="ECa V measurements",sub="")
points(df.v, pch = 19, col = df.v$ECaV, cex = 0.5)
plot(grd, cex = 1.5, col = "grey",axes=T)
points(df.h, pch = 19, col = df.h$ECaH, cex = 0.5)
legend( x="topleft",legend=c("Vertical","Horizontal"),col=c("red","blue"), lwd=1, lty=c(0,0),pch=19)


#saveRDS(grd, 'ny_1x1_grid.RData')

# Semi Variogram - Meters
v.ecv = variogram(ECaV~1, sp.v)
v.msv = variogram(MSaV~1, sp.v)
v.ech = variogram(ECaH~1, sp.h)
v.msh = variogram(MSaH~1, sp.h)

# Fit function to variogram
fit.ecv<-fit.variogram(v.ecv, vgm(c("Exp", "Mat", "Sph")))
fit.msv<-fit.variogram(v.msv, vgm(c("Exp", "Mat", "Sph")))
fit.ech<-fit.variogram(v.ech, vgm(c("Exp", "Mat", "Sph")))
fit.msh<-fit.variogram(v.msh, vgm(c("Exp", "Mat", "Sph")))

fit.ecv
fit.msv
fit.ech
fit.msh

range.ecv<-round(fit.ecv[2,3],2)
range.msv<-round(fit.msv[2,3],2)
range.ech<-round(fit.ech[2,3],2)
range.msh<-round(fit.msh[2,3],2)

# Plot each variogram
plot(v.ecv, fit.ecv,xlab="Distance (m)",main="ECa V",sub=paste0("range: ",range.ecv,"m"))
plot(v.msv, fit.msv,xlab="Distance (m)",main="MSa V",sub=paste0("range: ",range.msv,"m"))
plot(v.ech, fit.ech,xlab="Distance (m)",main="ECa H",sub=paste0("range: ",range.ech,"m"))
plot(v.msh, fit.msh,xlab="Distance (m)",main="MSa H",sub=paste0("range: ",range.msh,"m"))


# Kriging
print(paste("Start kriging ECaV 1x1 ~~~",Sys.time()))
kriged.ecv = krige(ECaV~1, sp.v, grd, model = fit.ecv)
#saveRDS(kriged.ecv, 'ny_ecv_1x1_fullgrid_fixgps.RData')

print(paste("Start kriging MSaV 1x1 ~~~",Sys.time()))
kriged.msv = krige(MSaV~1, sp.v, grd, model = fit.msv)
#saveRDS(kriged.msv, 'ny_msv_1x1_fullgrid_fixgps.RData')

print(paste("Start kriging ECaH 1x1 ~~~",Sys.time()))
kriged.ech = krige(ECaH~1, sp.h, grd, model = fit.ech)
#saveRDS(kriged.ech, 'ny_ech_1x1_fullgrid_fixgps.RData')

print(paste("Start kriging MSaH 1x1 ~~~",Sys.time()))
kriged.msh = krige(MSaH~1, sp.h, grd, model = fit.msh)
#saveRDS(kriged.msh, 'ny_msh_1x1_fullgrid_fixgps.RData')
print(paste("Done kriging ~~~",Sys.time()))




################################################


kriged.ecv<-readRDS('/data/ny_ecv_1x1_fullgrid_fixgps.RData')
kriged.msv<-readRDS('/data/ny_msv_1x1_fullgrid_fixgps.RData')
kriged.ech<-readRDS('/data/ny_ech_1x1_fullgrid_fixgps.RData')
kriged.msh<-readRDS('/data/ny_msh_1x1_fullgrid_fixgps.RData')

r.kriged.ecv <- raster(kriged.ecv) 
r.kriged.msv <- raster(kriged.msv) 
r.kriged.ech <- raster(kriged.ech) 
r.kriged.msh <- raster(kriged.msh) 

palc <- rev(brewer.pal(6,"Spectral"))
plot(r.kriged.ecv)#,col=palc)
#plot(kriged.ecv)

# get field perimeter area from shapefile
perimeter <- readOGR(dsn = "/data/", layer = "perimeter")
perimeter.df <- SpatialPolygonsDataFrame(perimeter,data=data.frame(row.names=row.names(perimeter)))
proj4string(perimeter)
proj4string(r.kriged.ecv)

perimeter <-spTransform(perimeter,CRS=utm) 
perimeter.df <-spTransform(perimeter.df,CRS=utm) 

r.kriged.ecv <- projectRaster(r.kriged.ecv, crs=utm)
r.kriged.msv <- projectRaster(r.kriged.msv, crs=utm)
r.kriged.ech <- projectRaster(r.kriged.ech, crs=utm)
r.kriged.msh <- projectRaster(r.kriged.msh, crs=utm)


plot(r.kriged.ecv,col=palc)
title("ECa V (log)", line = 0.5,cex.main=1)


plot(perimeter,add=T)
# extent(perimeter)
# extent(r.kriged.ecv)

# Crop using extent, rasterize polygon and finally, create poly-raster
cr.eca.v <- crop(r.kriged.ecv, extent(perimeter))
fr.eca.v <- rasterize(perimeter.df, cr.eca.v)   
lr.eca.v <- mask(x=cr.eca.v, mask=fr.eca.v)

cr.msa.v <- crop(r.kriged.msv, extent(perimeter))
fr.msa.v <- rasterize(perimeter.df, cr.msa.v)   
lr.msa.v <- mask(x=cr.msa.v, mask=fr.msa.v)

cr.eca.h <- crop(r.kriged.ech, extent(perimeter))
fr.eca.h <- rasterize(perimeter.df, cr.eca.h)   
lr.eca.h <- mask(x=cr.eca.h, mask=fr.eca.h)

cr.msa.h <- crop(r.kriged.msh, extent(perimeter))
fr.msa.h <- rasterize(perimeter.df, cr.msa.h)   
lr.msa.h <- mask(x=cr.msa.h, mask=fr.msa.h)

par(mfrow=c(2,2))
plot(lr.eca.v,col=palc,main="ECa V")
plot(lr.eca.h,col=palc,main="ECa H")
plot(lr.msa.v,col=palc,main="MSa V")
plot(lr.msa.h,col=palc,main="MSa H")

eca.v<-lr.eca.v
msa.v<-lr.msa.v
eca.h<-lr.eca.h
msa.h<-lr.msa.v


# field area 
sqm<-gArea(perimeter)
dunam<-sqm/1000
print(paste("Area of field:",round(dunam, digits=3),"dunam = ", round(sqm, digits=1),"m2"))



# Correlation - ECa 
rstack<-stack(lr.eca.v,lr.eca.h,lr.msa.v,lr.msa.h)
jnk=layerStats(rstack, 'pearson', na.rm=T)
corr_matrix=jnk$'pearson correlation coefficient'
colnames(corr_matrix)<-c('ECa V','ECa H','MSaV','MSa H')
rownames(corr_matrix)<-c('ECa V','ECa H','MSaV','MSa H')
corr_matrix

melted_cormat <- melt(corr_matrix)
melted_cormat$value<-round(melted_cormat$value,3)

par(mfrow=c(1,1),bty = 'n',mar=c(3, 2, 3, 2))
pald <- rev(sequential_hcl(7,h=c(20,60),c=120,l=c(30,90)))
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=factor(value))) + 
  geom_tile() + scale_fill_manual(values=pald) + theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=12)) + labs(title="Correlation Matrix ECa", caption = "Research Center 18-19/7/18") +geom_text(aes(label=value))





field.df <- readRDS('field_df.RData') # the data frame
feasible <- subset(field.df,feasible==1)
palg <- rev(colorRampPalette(brewer.pal(4,"Paired"))(4))
ggplot(field.df,aes(x,y)) + geom_tile(aes(fill=factor(fc)),alpha = 0.9) + coord_equal() + theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=12),plot.subtitle = element_text(size=10,hjust = 0.5),axis.text=element_text(size=10),panel.background = element_rect(fill = NA),panel.grid.major = element_line(size=0.1,colour = "grey30"),panel.grid.minor = element_line(size=0.1,colour = "grey30"),panel.border = element_rect(size=0.2, fill = NA),axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank()) + labs( y="", x="", title="Feasible search area") + scale_fill_manual(values=palg)

r.feasible<-feasible
coordinates(r.feasible) = ~x + y
proj4string(r.feasible) <- CRS(utm)
gridded(r.feasible) <- TRUE
field.r=raster(r.feasible,layer=8)
#writeRaster(field.r , filename="field_feasible.tif", format="GTiff", overwrite=TRUE)

mzs <- raster("/home/ai/Documents/Data/ny_MZs_smooth.tif")
proj4string(mzs) <- CRS(utm)
mzs <- projectRaster(mzs, crs=utm)

paln <- rev(colorRampPalette(brewer.pal(9,"Purples"))(10))[3:6]
par(bty="n")
plot(mzs,col=paln,ylab="Northing",xlab="Easting",bty="n",legend=FALSE,axes=F)
plot(field.r,legend=FALSE,add=T,col = alpha("#FFFFFF", 0.5),bty="n")
grid(lty = 1,col="grey70",lw=0.5)

title("Feasible area", line = 0.5,cex.main=1)
raster::scalebar(100, xy=c(705160, 3620840), type='bar', divs=2,below='m')









############### Get NDVI data #####################

### NDVI 

ndvi <- raster("/data/NDVI/Multi_Merge_index_ndvi.tif")



#origin(ndvi)
proj4string(ndvi)

res(eca.v)
res(eca.h)
res(ndvi)
# adjust the resolution 
ndvi.adj <- projectRaster(ndvi, crs=crs(utm),res=1)
#eca.v.adj <- projectRaster(lr.eca.v, crs=crs(utm),res=1)
res(ndvi.adj)

extent(ndvi.adj)
extent(eca.v)
extent(perimeter)

#pal <- heat_hcl(12)
par(mfrow=c(1,1))
#plot(ndvi,col=pal,main="NDVI - raw",line = 1,cex.main=1)
plot(ndvi.adj,col=palc,main="NDVI - rescaled",line = 1,cex.main=1)
plot(eca.v,col=palc,add=T,legend=F)
plot(eca.h,col=palc,add=T,legend=F)


### Crop
# Crop using extent, rasterize polygon and finally, create poly-raster
cr.ndvi <- crop(ndvi.adj, extent(perimeter))
fr.ndvi <- rasterize(perimeter.df, cr.ndvi)   
ndvi.c <- mask(x=cr.ndvi, mask=fr.ndvi)

# Plot - CROPPED 
par(mfrow=c(1,1),bty = 'n',mar=c(3, 4, 3, 4))
plot(ndvi.c,col=palc,main="NDVI", legend=TRUE)
plot(perimeter,add=T)




ndvi.min = cellStats(ndvi.c, "min")
ndvi.max = cellStats(ndvi.c, "max")
ndvi.n <- (ndvi.c - ndvi.min) / (ndvi.max - ndvi.min)

plot(ndvi.n,col=palc)
title("NDVI - normalized", line = 1,cex.main=1)
#writeRaster(ndvi.n , filename="ndvi_n_short.tif", format="GTiff", overwrite=TRUE)




############
# NDVI - 26/2/19
ndvi26 <- raster("/data/NDVI/26.2.19/ndvi/MultiFull_index_ndvi.tif")

#origin(ndvi)
proj4string(ndvi26)

res(eca.v)
res(eca.h)
res(ndvi26)
# adjust the resolution 
ndvi.adj26 <- projectRaster(ndvi26, crs=crs(utm),res=1)
#eca.v.adj <- projectRaster(lr.eca.v, crs=crs(utm),res=1)
res(ndvi.adj26)

extent(ndvi.adj26)
extent(eca.v)
extent(perimeter)

#pal <- heat_hcl(12)
par(mfrow=c(1,1))
#plot(ndvi,col=pal,main="NDVI - raw",line = 1,cex.main=1)
plot(ndvi.adj26,col=palc,main="NDVI - rescaled",line = 1,cex.main=1)
plot(eca.v,col=palc,add=T,legend=F)

### Crop
# Crop using extent, rasterize polygon and finally, create poly-raster
cr.ndvi <- crop(ndvi.adj26, extent(perimeter))
fr.ndvi <- rasterize(perimeter.df, cr.ndvi)   
ndvi.c26 <- mask(x=cr.ndvi, mask=fr.ndvi)
#writeRaster(ndvi.c26 , filename="ndvi.c26_short.tif", format="GTiff", overwrite=TRUE)

palgr <- brewer.pal(10,"Greens")
# Plot - CROPPED 
par(mfrow=c(1,2),bty = 'n',mar=c(3, 4, 3, 4))
plot(ndvi.c26,col=palgr,main="NDVI 26/2/19", legend=TRUE)
#plot(perimeter,add=T)
plot(ndvi.c,col=palgr,main="NDVI 10/12/18", legend=TRUE)






### LIDAR 

lidar <- raster("/data/LiDAR/Extract_dtm.tif")

proj4string(lidar)

res(eca.v)
res(eca.h)
res(lidar)
# adjust the resolution 
lidar.adj <- projectRaster(lidar, crs=crs(utm),res=1)
res(lidar.adj)

extent(lidar.adj)
extent(perimeter)

palid <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(60)
par(mfrow=c(1,1))
plot(lidar,col=palid,main="LIDAR - raw",cex.main=1)
plot(lidar.adj,col=palid,main="LIDAR - rescaled to resolution 1x1m ",cex.main=1)

### Crop
# Crop using extent, rasterize polygon and finally, create poly-raster
cr.lidar <- crop(lidar.adj, extent(perimeter))
fr.lidar <- rasterize(perimeter.df, cr.lidar)   
lidar.c <- mask(x=cr.lidar, mask=fr.lidar)

# Plot - CROPPED 
par(mfrow=c(1,1),bty = 'n',mar=c(3, 4, 3, 4))
plot(lidar.c,col=palid,main="LIDAR", legend=TRUE)
plot(perimeter,add=T)

#writeRaster(lidar.c , filename="lidar_c_short.tif", format="GTiff", overwrite=TRUE)



# flow accumulation

flowacc <- raster("/data/Flow/Flow_Acc_Neve_Yaar.tif")

proj4string(flowacc)

res(eca.v)
res(flowacc)
# adjust the resolution 
flowacc.adj <- projectRaster(flowacc, crs=crs(utm),res=1)
res(flowacc.adj)

extent(flowacc.adj)
extent(perimeter)

palid <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(60)
par(mfrow=c(1,1))
plot(flowacc.adj,col=palid,main="Flow Accumulation - rescaled to resolution 1x1m ",cex.main=1)

### Crop
# Crop using extent, rasterize polygon and finally, create poly-raster
cr.flowacc <- crop(flowacc.adj, extent(perimeter))
fr.flowacc <- rasterize(perimeter.df, cr.flowacc)   
flowacc.c <- mask(x=cr.flowacc, mask=fr.flowacc)

# Plot - CROPPED 
par(mfrow=c(1,1),bty = 'n',mar=c(3, 4, 3, 4))
plot(flowacc.c,col=palid,main="Flow accumulation", legend=TRUE)
plot(perimeter,add=T)
#writeRaster(flowacc.c , filename="flowacc_short.tif", format="GTiff", overwrite=TRUE)



### Thermal IR 

tir <- raster("/data/TIR images/Calibrated_projected_mosaics/Oct31_n_CPGR1.tif")

#origin(ndvi)
proj4string(tir)

res(eca.v)
res(tir)
# adjust the resolution 
tir.adj <- projectRaster(tir, crs=crs(utm),res=1)
res(tir.adj)

extent(tir.adj)
extent(eca.v)
extent(perimeter)

#pal <- heat_hcl(12)
par(mfrow=c(1,1))
#plot(ndvi,col=pal,main="NDVI - raw",line = 1,cex.main=1)
plot(tir.adj,col=palc,main="Thermal IR - rescaled",line = 1,cex.main=1)
plot(eca.v,col=palc,add=T,legend=F)


### Crop
# Crop using extent, rasterize polygon and finally, create poly-raster
cr.tir <- raster::crop(tir.adj, extent(perimeter))
fr.tir <- rasterize(perimeter.df, cr.tir)   
tir.c <- mask(x=cr.tir, mask=fr.tir)

# Plot - CROPPED 
par(mfrow=c(1,1),bty = 'n',mar=c(3, 4, 3, 4))
plot(tir.c,col=palc,main="Thermal IR 31/10/19", legend=TRUE)
plot(perimeter,add=T)



tir.min = cellStats(tir.c, "min")
tir.max = cellStats(tir.c, "max")
tir.n <- (tir.c - tir.min) / (tir.max - tir.min)

palgr <- brewer.pal(9,"Reds")
plot(tir.n,col=palgr)
title("TIR - normalized", line = 1,cex.main=1)
#writeRaster(tir.n , filename="tir_n_311019_short.tif", format="GTiff", overwrite=TRUE)




### raster to matrix

coo<-data.frame(rasterToPoints(eca.v)[,1:2])
ecv<-rasterToPoints(eca.v)[,3]
msv<-rasterToPoints(msa.v)[,3]
ech<-rasterToPoints(eca.h)[,3]
msh<-rasterToPoints(msa.h)[,3]

# Adjust resolution to reach same number of cells
ndvi2 <- projectRaster(ndvi.c, crs=crs(utm),res=1.001526)
length(ecv)
endvi<-rasterToPoints(ndvi.c)[,3]
length(endvi)
ndvi.e<-endvi[1:length(ecv)]

# Adjust resolution to reach same number of cells
ndvi26.2 <- projectRaster(ndvi.c26, crs=crs(utm),res=1.001507)
length(ecv)
endvi26<-rasterToPoints(ndvi26.2)[,3]
length(endvi26)
ndvi26.e<-endvi26[1:length(ecv)]

lidar2 <- projectRaster(lidar.c, crs=crs(utm),res=1.001505)
length(ecv)
elidar<-rasterToPoints(lidar2)[,3]
length(elidar)
lidar.e<-elidar[1:length(ecv)]

fa2 <- projectRaster(flowacc.c, crs=crs(utm),res=1.001515)
length(ecv)
efa<-rasterToPoints(fa2)[,3]
length(efa)
flowacc.e<-efa[1:length(ecv)]


tir2 <- projectRaster(tir.c, crs=crs(utm),res=1.001515)
length(ecv)
etir<-rasterToPoints(tir2)[,3]
length(etir)
tir.e<-etir[1:length(ecv)]

summary(as.data.frame(ndvi.c))
summary(as.data.frame(tir.c))
res(tir.c)
res(ndvi.c)
# Dataframe with Coordinates and all variables columns
data.all<-data.frame(coo[,1],coo[,2],ecv,msv,ech,msh,ndvi.e,ndvi26.e,lidar.e,flowacc.e,tir.e)
data.all<-na.omit(data.all)
colnames(data.all)<-c('x','y','ecv','msv','ech','msh','ndvi','ndvi26','lidar','flowacc','tir3110')

### normalize - Min-Max
data.n <-data.frame(coo[,1],coo[,2],lapply(data.all[3:11], normal))
colnames(data.n)<-c('x','y','ecv','msv','ech','msh','ndvi','ndvi26','lidar','flowacc','tir3110')
summary(data.n)





# # # # # # # # # # # # # # # # # # # # # # # # ######################
### Clustering 



##
# Fuzzy clusters validitiy indices
x <- data.n[c(3:7)] # data columns for clustering
caption="ECa, MSa, NDVI 12/18"


max.k=12
fuzzy.indices<-data.frame(k=2:max.k)
pe.v<-pc.v<-fs.v<-chc.v<-list()

# Partition entropy
for (k in 2:max.k){
  fc <- cmeans(x,k,1000,verbose=F,m=2,method="cmeans")
  pe.v <- append(pe.v,fclustIndex(fc,x, index="partition.entropy"))
  pc.v <- append(pc.v,fclustIndex(fc,x, index="partition.coefficient"))
  fs.v <- append(fs.v,fclustIndex(fc,x, index="fukuyama.sugeno"))  
  chc.v <- append(chc.v,calinhara(x,fc$cluster))
  print(k)
}
fuzzy.indices<-cbind(fuzzy.indices,pe=unlist(pe.v),pc=unlist(pc.v),fs=unlist(fs.v),chc=unlist(chc.v))


###

theme.fvi <- theme(plot.subtitle = element_text(hjust=1,size=11),panel.background = element_rect(NA),panel.grid = element_line(size=0.1,colour = "grey30"),panel.border = element_rect(size=0.1, fill = NA))

fvi.1 <- ggplot(fuzzy.indices,aes(k,pe)) + geom_point(size=1) + geom_line(size=1) + theme.fvi + labs(title="Partition entropy (max)",caption=caption) + scale_x_discrete("k", 1:max.k, 1:max.k, 1:max.k)
fvi.2 <- ggplot(fuzzy.indices,aes(k,pc)) + geom_point(size=1) + geom_line(size=1) + theme.fvi + labs(title="Partition coefficient (min)",caption=caption) + scale_x_discrete("k", 1:max.k, 1:max.k, 1:max.k)
fvi.3 <- ggplot(fuzzy.indices,aes(k,fs)) + geom_point(size=1) + geom_line(size=1) + theme.fvi + labs(title="Fukuyama-Sugeno (min)",caption=caption) + scale_x_discrete("k", 1:max.k, 1:max.k, 1:max.k)
fvi.4 <- ggplot(fuzzy.indices,aes(k,chc)) + geom_point(size=1) + geom_line(size=1) + theme.fvi + labs(title="Calinski-Harbaz criterion (max)",caption=caption) + scale_x_discrete("k", 1:max.k, 1:max.k, 1:max.k)

multiplot(fvi.1,fvi.2,fvi.3,fvi.4,cols = 2)
multiplot(fvi.3,fvi.4)
##################


# Fuzzy clusters validitiy indices
x <- data.n[c(3:7)] # data columns for clustering
caption="TIR 3/20, ECa MSa 12/19"
plot(tir.c)

### Fuzzy c-means clustering
fc <- cmeans(x,4,1000,verbose=T,m=2,method="cmeans")
data.n$fc <- fc$cluster
data.n$fc <- mapvalues(data.n$fc, from=c(unique(data.n$fc)), to=c(index(unique(data.n$fc))))

palcm<-rev(brewer.pal(9,"Spectral")[c(1,5,9,7)])
ggplot(data.n,aes(x,y)) + geom_tile(aes(fill=factor(fc))) + coord_equal() + theme_bw() + theme(plot.title = element_text(hjust = 0.5,size=12)) + labs( y="Northing", x="Easting", title=paste0("Clusters by ",caption), caption = "") + scale_fill_manual(values=palcm,name="cluster",labels=c("1", "2", "3", "4", "5", "6","7","8"))




### Save to raster
# Set spatial coordinates to create a Spatial object:
field <- data.n
coordinates(field) = ~x + y
proj4string(field) <- CRS(utm)
field.r=raster(field)


gg1 <- 
  ggplot(data.n,aes(x,y)) + geom_tile(aes(fill=ndvi26)) + coord_equal() + theme_bw() + theme(plot.title = element_text(hjust = 0.5,size=12)) + labs( y="", x="", title="ECa V") + scale_fill_gradientn(colours = rev(brewer.pal(20,"Spectral")))

gg2 <- ggplot(data.n,aes(x,y)) + geom_tile(aes(fill=ech)) + coord_equal() + theme_bw() + theme(plot.title = element_text(hjust = 0.5,size=12)) + labs( y="", x="", title="ECa H") + scale_fill_gradientn(colours = rev(brewer.pal(20,"Spectral")))

gg3 <- ggplot(data.n,aes(x,y)) + geom_tile(aes(fill=msv)) + coord_equal() + theme_bw() + theme(plot.title = element_text(hjust = 0.5,size=12)) + labs( y="", x="", title="MSa V") + scale_fill_gradientn(colours = rev(brewer.pal(20,"Spectral")))

gg4 <- ggplot(data.n,aes(x,y)) + geom_tile(aes(fill=msh)) + coord_equal() + theme_bw() + theme(plot.title = element_text(hjust = 0.5,size=12)) + labs( y="", x="", title="MSa H") + scale_fill_gradientn(colours = rev(brewer.pal(20,"Spectral")))

gg5 <- ggplot(data.n,aes(x,y)) + geom_tile(aes(fill=ndvi)) + coord_equal() + theme_bw() + theme(plot.title = element_text(hjust = 0.5,size=12)) + labs( y="", x="", title="NDVI 12/18") + scale_fill_gradientn(colours = rev(brewer.pal(20,"Spectral")))


gg6 <- ggplot(data.n,aes(x,y)) + geom_tile(aes(fill=ndvi26)) + coord_equal() + theme_bw() + theme(plot.title = element_text(hjust = 0.5,size=12)) + labs( y="", x="", title="NDVI 2/19") + scale_fill_gradientn(colours = rev(brewer.pal(20,"Spectral")))

gg7 <- ggplot(data.n,aes(x,y)) + geom_tile(aes(fill=lidar)) + coord_equal() + theme_bw() + theme(plot.title = element_text(hjust = 0.5,size=12)) + labs( y="", x="", title="LIDAR") + scale_fill_gradientn(colours = rev(brewer.pal(20,"Spectral")))

gg8 <- ggplot(data.n,aes(x,y)) + geom_tile(aes(fill=flowacc)) + coord_equal() + theme_bw() + theme(plot.title = element_text(hjust = 0.5,size=12)) + labs( y="", x="", title="Flow accumulation") + scale_fill_gradientn(colours = rev(brewer.pal(20,"Spectral")))

gg9 <- 
  ggplot(data.n,aes(x,y)) + geom_tile(aes(fill=ndvi26)) + coord_equal() + theme_bw() + theme(plot.title = element_text(hjust = 0.5,size=12)) + labs( y="", x="", title="Thermal IR 31/10/19") + scale_fill_gradientn(colours = rev(brewer.pal(20,"Spectral")))



multiplot(gg1,gg2,gg3,gg4,gg5,gg6,gg7,gg8,gg9,cols=2)
