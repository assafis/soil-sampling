library(entropy)
#' @importFrom graphics hist
#' 
# Objective functions for ECR: d-kl, clhs

.dkl_obj <- function(
  pop,
  df,
  d.cols,
  continuous_strata,
  cor_mat,
  N
) {
  #data_continuous_sampled<-x
  #n_cont_variables <- ncol(continuous_data)
  #print(continuous_data)
  data_sampled<-df[pop,]
  # print(data_continuous_sampled)
  
  
  full.ecv<-df$ecv
  h.full <- hist(full.ecv, breaks = seq(0,1,by=1/(N/2-1)),plot=FALSE)
  h.full$counts[h.full$counts==0]<-0.01
  h.full$counts=h.full$counts/sum(h.full$counts)
  d.full<-h.full$counts
  
  samp.ecv<-data_sampled$ecv
  h.samp <- hist(samp.ecv, breaks = seq(0,1,by=1/(N/2-1)), plot=FALSE)
  h.samp$counts[h.samp$counts==0]<-0.01
  h.samp$counts=h.samp$counts/sum(h.samp$counts)
  d.samples<-h.samp$counts
  
  
  ## DKL(P???Q) ##
  # P typically represents the "true" distribution of data.
  # Q typically represents a model / approximation of P
  
  kl <- signif(KL.plugin(d.full,d.samples),3)
  # plot(h.full, col="orange",main="Histograms")
  # plot(h.samp, col="orange",main=paste("D-KL:",kl))
  
  return(kl)
  
}




#' Objective function for cLHS ECR
#'
#' adapted from Pierre Roudier clhs package
#' 

.lhs_obj <- function(
  pop,
  df,
  d.cols,
  continuous_strata,
  cor_mat,
  N
) {
  #pop=population[[1]]
  #df=x
  continuous_data<-df[d.cols]
  n_cont_variables <- ncol(continuous_data)
  data_continuous_sampled<-continuous_data[pop,]
  
  cont_data_strata <- lapply(1:n_cont_variables, function(i) list(data_continuous_sampled[, i], continuous_strata[, i]) )
  cont_obj_sampled <- lapply(cont_data_strata, function(x) hist(x[[1]], breaks = x[[2]], plot = FALSE)$counts)
  cont_obj_sampled <- matrix(unlist(cont_obj_sampled), ncol = n_cont_variables, byrow = FALSE)
  
  delta_obj_continuous <- rowSums(abs(cont_obj_sampled - 1))
  
  # Correlation of continuous data
  cor_sampled <- suppressWarnings(cor(data_continuous_sampled))
  cor_sampled[is.na(cor_sampled)] <- 1 
  delta_obj_cor <- sum(abs(cor_mat - cor_sampled))
  
  obj <- sum(delta_obj_continuous) + delta_obj_cor # Objective function
  #obj <- sum(delta_obj_continuous)
  obj.n <- signif(obj/(N*n_cont_variables),3) # normalized to N
  
  return(obj.n)
}





library(gstat)
library(raster)
RMSE <- function(x, y) { sqrt(mean((x - y)^2)) } 

.rmse_obj <- function(
  pop,
  df,
  d.cols,
  continuous_strata,
  cor_mat,
  N
) {
  
  data_sampled<-df[pop,]
  
  grd10 <- readRDS('grid10.Rdata') # grid
  hull <- readRDS('hull.Rdata') # polygon perimeter
  hull.df <- readRDS('hull_df.Rdata') # polygon perimeter df
  full.ecv.n <- readRDS('ecv_cr_n.Rdata') # ECa V raster
  
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  # ecvr.min = cellStats(lr.eca.v, stat='min')
  # ecvr.max = cellStats(lr.eca.v, stat='max')
  # full.ecv.n = (lr.eca.v-ecvr.min)/(ecvr.max-ecvr.min)
  # saveRDS(full.ecv.n,'ecv_cr_n.Rdata')
  # 
  coordinates(data_sampled) = ~x + y
  utm <- "+proj=utm +zone=36S +datum=WGS84"
  proj4string(data_sampled) <- CRS(utm)
  
  v1 = variogram(ecv~1, data_sampled)
  fitv<-fit.variogram(v1, vgm("Sph"))
  
  ecv.kr = krige(ecv~1, data_sampled, grd10, model = fitv,debug.level = 0)
  
  r.ecv.kr <- raster(ecv.kr)
  
  cr.eca.v <- crop(r.ecv.kr, extent(hull))
  fr.eca.v <- rasterize(hull.df, cr.eca.v)   
  samp.ecv <- mask(x=cr.eca.v, mask=fr.eca.v)
  
  #  plot(samp.ecv)
  #  plot(full.ecv.n)
  #  plot(full.ecv.n-samp.ecv)
  
  options(warn = oldw)
  
  diff=full.ecv.n-samp.ecv
  diff.sq=diff^2
  diff.avg = cellStats(diff.sq, stat='mean')
  diff.rt = sqrt(diff.avg)
  
  return(diff.rt)
  
  
  
}




.kvar_obj <- function(
  pop,
  df,
  d.cols,
  continuous_strata,
  cor_mat,
  N
) {
  
  data_sampled<-df[pop,]
  
  grd10 <- readRDS('grid10.Rdata') # grid
  hull <- readRDS('hull.Rdata') # polygon perimeter
  hull.df <- readRDS('hull_df.Rdata') # polygon perimeter df
  
  coordinates(data_sampled) = ~x + y
  utm <- "+proj=utm +zone=36S +datum=WGS84"
  proj4string(data_sampled) <- CRS(utm)
  
  oldw <- getOption("warn") 
  options(warn = -1) # hide variogram warnings
  v1 = variogram(ecv~1, data_sampled)
  #class(v1)
  fitv<-fit.variogram(v1, vgm("Sph"))
  #plot(v1,fitv,cloud=TRUE)
  options(warn = oldw)
  
  ecv.kr = krige(ecv~1, data_sampled, grd10, model = fitv,debug.level = 0)
  r.ecv.krv <- raster(ecv.kr["var1.var"])
  
  cr.eca.var <- crop(r.ecv.krv, extent(hull))
  fr.eca.var <- rasterize(hull.df, cr.eca.var)   
  samp.ecvar <- mask(x=cr.eca.var, mask=fr.eca.var)
  
  mean.var = cellStats(samp.ecvar, stat='mean')
  #print(mean.var)
  if (length(mean.var)<1 | !is.numeric(mean.var)){mean.var=1}
  return(mean.var)
}



.ks_obj <- function(
  pop,
  df,
  d.cols,
  continuous_strata,
  cor_mat,
  N
) {
  continuous_data<-df[d.cols]
  n_cont_variables <- ncol(continuous_data)
  data_continuous_sampled<-continuous_data[pop,]
  
  ks.mean=mean(unlist(lapply(1:n_cont_variables, function(i) {
    ks.test(data_continuous_sampled[,i],continuous_data[,i])$statistic[[1]]
  })))
  
  return(ks.mean)
}




.chisq_obj <- function(pop,df,d.cols,continuous_strata,cor_mat,N) {
  data_sampled<-df[pop,]
  full.ecv<-df$ecv
  h.full <- hist(full.ecv, breaks = seq(0,1,by=1/(N/2-1)),plot=FALSE)
  h.full$counts[h.full$counts==0]<-0.01
  h.full$counts=h.full$counts/sum(h.full$counts)
  d.full<-h.full$counts
  
  samp.ecv<-data_sampled$ecv
  h.samp <- hist(samp.ecv, breaks = seq(0,1,by=1/(N/2-1)), plot=FALSE)
  h.samp$counts[h.samp$counts==0]<-0.01
  h.samp$counts=h.samp$counts/sum(h.samp$counts)
  d.samples<-h.samp$counts
  
  chisq <- signif(chi2.plugin(d.full,d.samples),3)
  #plot(h.full, col="orange",main="Full field")
  #plot(h.samp, col="orange",main="Samples", caption=paste("D-KL:",chisq))
  
  return(chisq)
  
}


.dist_obj <- function(pop,df,d.cols,continuous_strata,cor_mat,N,theta) {
  popxy<-df[pop,1:2]
  dst<-1/(min(dist(popxy)))
  return(dst)
}


# pop=population[[1]]
# df=x
# j=1
# i=4

.SPdiv.xy_obj <- function(pop,df,d.cols,continuous_strata,cor_mat,N,theta) {
  popxy<-df[pop,1:2]
  D=as.matrix(dist(popxy))
  lend=dim(D)[1]
  M=diag(lend)
  for (i in 1:lend-1){
    for (j in (i+1):lend){
      M[i,j]=exp(theta*D[i,j])
      M[j,i]=M[i,j]
    }
  }
  SPdiv=sum(solve(M))
  return(1/SPdiv)
}


.SPdiv.anc_obj <- function(pop,df,d.cols,continuous_strata,cor_mat,N,theta) {
  popxy<-df[pop,d.cols]
  D=as.matrix(dist(popxy))
  lend=dim(D)[1]
  M=diag(lend)
  for (i in 1:lend-1){
    for (j in (i+1):lend){
      #print(paste("i",i," , j",j))
      M[i,j]=exp(theta*D[i,j])
      M[j,i]=M[i,j]
    }
  }
  SPdiv=sum(solve(M))
  return(1/SPdiv)
}





#function SPdiv = SP_diversityMeasure(D)
# SP_diversityMeasure calculates the Solow-Polasky diversity measure,
# based upon A.R. Solow and S. Polasky, Measuring biological diversity.
# Environmental and Ecological Statistics, 1:95-103, 1994
#
# SPdiv = SP_diversityMeasure(D)
#
#  Input and output arguments:
#   D     (matrix) data: the distances matrix of the input data
#   SPdiv   (scalar) the Solow-Polasky diversity measure for the data
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%% EOF %%%
# theta=-10.0;
# M=eye(size(D));
# for i=1:length(D),
# for j=i+1:length(D),
# M(i,j)=exp(theta*D(i,j));
# M(j,i)=M(i,j);
# end
# end
# SPdiv=sum(sum(inv(M)));
#%%%%%%%%%%%%%%%%%%



