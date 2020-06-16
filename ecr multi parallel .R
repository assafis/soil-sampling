library(parallel)
library(parallelMap)
library(ecr)
library(ggplot2)
library(RColorBrewer)
library(pracma)

source('ecr multi runs inner function.R')
source('multiplot_func.R' ) 
source("ranking_functions.R")
source('spatial_helper_functions.R')


MU = 10
LAMBDA = 10
nsamp=10
MAX.ITER=50000
d.cols = c(3:7) # ancillary data columns
theta=-20
fitness.names = c("cLHS","max-min-diversity")
#fitness.names = c("SP xy","SP ancillary")
threads=30


parallelStop()
### run in parallel

npoints<-c(10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)
thetas<-c(-0.01,-0.05,-0.1,-0.2,-0.5,-1,-2,-5,-10)

nsamp=10
i=1
for(theta in thetas){
  
  #nsamp=np
  print(paste("Start ~",Sys.time(),'~~~',nsamp,"points /",MAX.ITER,"iterations/",threads,"threads"))
  parallelStart(mode = "multicore", cpus=threads)
 # assign(paste0("sp.theta.",i),parallelMap(inner, 1:threads, more.args=list(MU,LAMBDA,MAX.ITER,nsamp,d.cols,theta)))
  parallelStop()
  #    saveRDS(eval(parse(text=paste0("fr.",np))), paste0('fr_clhs_mmd_',np,'.RData'))
#  saveRDS(eval(parse(text=paste0("sp.theta.",i))), paste0('fr_sp_theta_',i,'.RData'))
  print(paste("End ~~~",Sys.time()))
  i<-i+1
}




###################

field.df <- readRDS('field_df.RData') # the data frame
feasible <- subset(field.df,feasible==1)


# read runs results from file
npoints<-c(10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)
for(np in npoints){
  assign(paste0("fr.",np),readRDS(paste0("fr_clhs_mmd_",np,".RData")))
}


# stage all runs in dataframe, ranking, fronts
fid=1
for(np in npoints){
  fr<-eval(parse(text=paste0("fr.",np)))
  frind<-lapply(fr, `[[`, 2)
  
  pop.fit<-data.frame(double(),double(),integer())
  colnames(pop.fit)<-c("obj1","obj2","pop")
  
  for(i in 1:length(frind)){
    for(j in 1:length(frind[[1]][[1]])){
      pop.fit[(i-1)*10+j,1]=frind[[i]][[2]][1,j]
      pop.fit[(i-1)*10+j,2]=frind[[i]][[2]][2,j]
      pop.fit[(i-1)*10+j,3]=list(frind[[i]][[1]][j])
    }
  }
  
  pop.fit$ID <- seq.int(nrow(pop.fit))
  pop.fit$rank<-as.factor(ranking(pop.fit[1:2]))
  
  front<-pop.fit[pop.fit$rank==0,]
  front$frontID <- seq(fid,(fid+nrow(front)-1),by=1)
  fid<-fid+nrow(front)
  front<-front[!duplicated(front[, -which(names(front) %in% c("ID"))]), ] # remove duplicate solutions
  
  assign(paste0("front.",np),front)
  assign(paste0("pop.fit.",np),pop.fit)
  print(paste("done",np))
}

# bind all fronts
df.names<-paste0("front.",npoints)
all.fronts<-do.call("rbind", as.list(parse(text=paste0(df.names))))
all.fronts<-all.fronts[!duplicated(all.fronts[, -which(names(all.fronts) %in% c("ID","rank","frontID"))]), ] # remove duplicates

all.fronts$n<-NA # add n column
for(i in 1:nrow(all.fronts)){ 
  all.fronts[i,"n"]<-length(all.fronts$pop[[i]])
  all.fronts$frontID[i]<-i
}
all.fronts$n<-factor(all.fronts$n)

#saveRDS(all.fronts, 'ny_all.fronts_clhs_mmd_10-50.RData')

maxx<-max(all.fronts$obj1)
maxy<-max(all.fronts$obj2)

steps<-length(unique(all.fronts$n))
palt <- colorRampPalette(brewer.pal(11,"Set1"))(steps)

#gg.fronts.sp<-
  ggplot(all.fronts)+geom_point(aes(obj1,obj2,col=n),size=0.5) + 
  theme(legend.position="right",legend.key = element_rect(fill = "transparent", colour = "transparent"),plot.title = element_text(hjust = 0,size=10),plot.subtitle = element_text(hjust = 0.5),panel.background = element_rect(NA),panel.grid.major = element_line(size=0.1,colour = "grey90"),panel.grid.minor = element_line(size=0.1,colour = "grey80"),panel.border = element_rect(size=0.1, fill = NA)) +
  labs( x=paste("(<-min)",fitness.names[1]), y=paste("(<-min)",fitness.names[2]), title=paste0("(",MU,"+",LAMBDA,")ES ",threads,"x",MAX.ITER," iterations"), subtitle=paste0("Pareto fronts by sample size (n)")) + geom_line(aes(obj1,obj2,col=n),size=0.5) + scale_x_continuous(expand = c(0, 0),limits=c(0, maxx+(0.1*maxx))) + scale_y_continuous(expand = c(0, 0),limits=c(0, maxy+(0.1*maxy))) + scale_color_manual(values=palt)



#ggsave("sp_fronts", plot = gg.fronts.sp, device = "eps", scale = 1,dpi = 300)


# plot all results of one n
pop.fit<-pop.fit.26
nsamp<-length(pop.fit$pop[[1]])
steps<-length(unique(pop.fit$rank))
palf <- colorRampPalette(brewer.pal(10,"Paired"))(steps)

maxx.p<-max(pop.fit$obj1)
maxy.p<-max(pop.fit$obj2)
ggplot(pop.fit)+geom_point(aes(obj1,obj2,col=rank),size=0.8) + 
  labs( x=paste(fitness.names[1]), y=paste(fitness.names[2]), subtitle=paste0(nsamp, " points"), title=paste0("(",MU,"+",LAMBDA,")NSGA-II ",threads,"x",MAX.ITER," iterations")) + 
  geom_line(data=front.26,aes(obj1,obj2),col="grey15") +
  #geom_text(data=front.26,aes(obj1,obj2,label=ID,hjust=1.4, vjust=2)) +
  theme(legend.position="right",legend.key = element_rect(fill = "transparent", colour = "transparent"),plot.title = element_text(hjust = 0,size=13),plot.subtitle = element_text(hjust = 0.5,size=14),panel.background = element_rect(NA),axis.title=element_text(size=14),axis.text=element_text(size=12),panel.grid.major = element_line(size=0.1,colour = "grey90"),panel.grid.minor = element_line(size=0.1,colour = "grey80"),panel.border = element_rect(size=0.1, fill = NA)) + scale_x_continuous(expand = c(0, 0),limits=c(0, maxx.p+(0.1*maxx.p))) + scale_y_continuous(expand = c(0, 0),limits=c(0, maxy.p+(0.1*maxy.p))) + scale_color_manual(values=palf)




pl.front<-ggplot(front.26)+geom_point(aes(obj1,obj2),size=2) + 
  theme(legend.position="top",legend.key = element_rect(fill = "transparent", colour = "transparent"),plot.title = element_text(hjust = 0.5,size=12),plot.subtitle = element_text(hjust = 0.5),panel.background = element_rect(NA),panel.grid.major = element_line(size=0.1,colour = "grey80"),panel.grid.minor = element_line(size=0.1,colour = "grey80"),panel.border = element_rect(size=0.2, fill = NA)) +
  labs( x=paste("(<-min)",fitness.names[1]), y=paste("(<-min)",fitness.names[2]), title=paste0("Pareto front - ",nsamp, " points"), subtitle=paste0("(",MU,"+",LAMBDA,")ES ",threads,"x",MAX.ITER," iterations")) + 
  geom_line(data=front.26,aes(obj1,obj2)) + scale_color_brewer(palette = "Set1") +  geom_text(aes(obj1,obj2,label=front.26$frontID,hjust=-0.2, vjust=-0.3))
pl.front



pop.fit<-pop.fit.26
pop.fit$dist<-pop.fit$obj2
pop.front<-front.26
pop.front$dist<-pop.front$obj2
maxy.d<-max(pop.fit$dist)

nsamp<-length(pop.fit$pop[[1]])
steps<-length(unique(pop.fit$rank))
palf <- colorRampPalette(brewer.pal(10,"Paired"))(steps)



palgr <- colorRampPalette(rev(brewer.pal(9,"Greys")[3:9]))(steps)
#palgr <- colorRampPalette(rev(brewer.pal(9,"Greys")[3:9]))(steps)

ggplot(pop.fit)+geom_point(aes(obj1,dist,col=rank),size=0.9) + 
  labs( x=paste(fitness.names[1]), y=paste(fitness.names[2]), subtitle=paste0(nsamp, " points"), title=paste0("(",MU,"+",LAMBDA,")NSGA-II ",threads,"x",MAX.ITER," iterations")) + 
  geom_line(data=pop.front,aes(obj1,dist),col="black") +  #c(0.89, 0.6)
  theme(legend.key = element_rect(fill = "transparent", colour = "transparent"),legend.position = "right",plot.title = element_text(hjust = 0,size=10),legend.background = element_rect(fill="white",size=0.1, linetype="solid", colour ="transparent"),plot.subtitle = element_text(hjust = 0.5),panel.background = element_rect(NA),panel.grid.major = element_line(size=0.1,colour = "grey90"),panel.grid.minor = element_line(size=0.1,colour = "grey80"),panel.border = element_rect(size=0.1, fill = NA)) + scale_x_continuous(expand = c(0, 0),limits=c(0, maxx.p+(0.1*maxx.p))) + scale_y_continuous(expand = c(0, 0),limits=c(0, maxy.d+(0.1*maxy.d))) + scale_color_manual(values=palf) 
# + annotate("text", x = 0.08, y = 0.01, label = "front",col="brown")

ggsave("solutions26", plot = pp26, device = "eps", scale = 1,dpi = 300)



### plot Hyper Volume Indicator progress
fr<-fr.26
nsamp=26
frhv<-lapply(fr, `[[`, 1)
hvruns=do.call("cbind", frhv)
gen<-hvruns[,1]
hvruns<-hvruns[c(FALSE,TRUE)]
n.runs<-ncol(hvruns)
hvall<-cbind(gen,hvruns)

hvall$hv.min<-apply(hvall[2:(n.runs+1)],1,min)
hvall$hv.max<-apply(hvall[2:(n.runs+1)],1,max)
hvall$hv.mean<-apply(hvall[2:(n.runs+1)],1,mean)


ggplot(hvall,aes(gen,hv.max)) + geom_line(aes(lty="a")) + 
  labs(x='generation',y='HV',title=paste("Hyper-Volume Indicator -",nsamp,"points")) + 
  geom_line(aes(gen,hv.mean,lty="b")) + 
  geom_line(aes(gen,hv.min,lty="c")) + 
  theme(legend.position=c(0.91,0.2),legend.key = element_rect(fill = "transparent", colour = "transparent"),plot.title = element_text(hjust = 0.5,size=12),plot.subtitle = element_text(hjust = 0.5),panel.background = element_rect(NA),panel.grid.major = element_line(size=0.2,colour = "grey80"), panel.grid.minor = element_line(size=0.2,colour = "grey80"),panel.border = element_rect(size=0.1, fill = NA)) + 
  scale_linetype_manual(name="", values=c("a"=1,"b"=4,"c"=3), labels = c("Max","Mean","Min")) 



# assign samples of front members
for (k in 1:nrow(all.fronts)){
  sample<-feasible[all.fronts[[k,3]],]  # sample points with full data
  sample$sampleID<-all.fronts[k,"ID"]
  sample$frontID<-k
  assign(paste0("sample", k),sample)
}


#csv.pts<-read.csv("/home/ai/Documents/Data/neve_yaar_22_sampling_coords_21.csv")
#df.pts = as.data.frame(csv.pts)
r.mz.smooth<-raster("/home/assafi/Documents/Data/Neve Yaar/ny_MZs_smooth.tif")
plot(r.mz.smooth)
mz.smooth<-as.data.frame(rasterToPoints(r.mz.smooth))
mz.smooth<-mz.smooth[1:nrow(field_df),]
colnames(mz.smooth)[3]<-"fc"
palg <- rev(brewer.pal(9,"PRGn"))
palc <- rev(brewer.pal(4,"Spectral"))
pale <- colorRampPalette(brewer.pal(11,"Spectral"))(20)

sample<-df.pts
sid<-sample[1,"sampleID"]
fid<-21#sample[1,"frontID"]
np=nrow(sample)

palg <- rev(colorRampPalette(brewer.pal(4,"PRGn"))(4))
palg <- brewer.pal(9,"PRGn")[c(3,1,4,2)]
ggplot(mz.smooth,aes(x,y)) + geom_tile(aes(fill=factor(fc)),alpha = 0.7) + coord_equal() + theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=12),plot.subtitle = element_text(size=11,hjust = 0.5),axis.text=element_text(size=10),panel.background = element_rect(fill = NA),panel.grid.major = element_line(size=0.1,colour = "grey30"),panel.grid.minor = element_line(size=0.1,colour = "grey30"),panel.border = element_rect(size=0.2, fill = NA),axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank()) + labs( y="Northing", x="Easting", title="", subtitle=paste0(np," points / Sample #",fid)) + geom_point(data=sample, aes(x, y),size=1) + scale_fill_manual(values=palg)



#################

### Scan indices by sample size

#########
### AIC - by sample size


aic.res <- data.frame(aic.ecv=numeric(),aic.msv=numeric(),aic.ech=numeric(),aic.msh=numeric(),aic.ndvi=numeric(),sid=character(),fid=integer(),n=integer())


for (i in 1:nrow(all.fronts)){
  samp<-eval(parse(text=paste0("sample",i)))
  
  lm <- lm(ecv ~ . , samp[d.cols])
  aic.ecv <- AIC(lm)
  
  lm <- lm(msv ~ . , samp[d.cols])
  aic.msv <- AIC(lm)
  
  lm <- lm(ech ~ . , samp[d.cols])
  aic.ech <- AIC(lm)
  
  lm <- lm(ecv ~ . , samp[d.cols])
  aic.ecv <- AIC(lm)
  
  lm <- lm(msh ~ . , samp[d.cols])
  aic.msh <- AIC(lm)
  
  lm <- lm(ndvi ~ . , samp[d.cols])
  aic.ndvi <- AIC(lm)
  
  sid<-samp[1,"sampleID"]
  fid<-samp[1,"frontID"]
  n=nrow(samp)
  newrow<-data.frame(aic.ecv,aic.msv,aic.ech,aic.msh,aic.ndvi,sid,fid,n)
  aic.res<-rbind(aic.res,newrow)
}

aic.res
aic.res$aic.mean<-apply(aic.res[1:5],1,mean)
aic.res$n<-as.factor(aic.res$n)
aic.agg <-aggregate(aic.res[,c(1:5,9)], list(n=aic.res$n),mean)
aic.agg.min <-aggregate(aic.res[,c(1:5,9)], list(n=aic.res$n),min)
aic.agg.max <-aggregate(aic.res[,c(1:5,9)], list(n=aic.res$n),max)
aic.agg$aic.mean.min<-aic.agg.min$aic.mean
aic.agg$aic.mean.max<-aic.agg.max$aic.mean

#pl.aic.all<-
ggplot(aic.res) + 
  theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=12),plot.subtitle = element_text(size=10,hjust = 0.5),axis.text=element_text(size=12),panel.background = element_rect(fill = NA),panel.grid.major = element_line(size=0.1,colour = "grey50"),panel.border = element_rect(size=0.3, fill = NA,colour = "black")) +
  scale_x_discrete() +
  labs( x="n",y="AIC", title='AIC by sample size (n)') +
  geom_point(aes(n,aic.ecv),col="midnightblue",alpha=0.5,size=2) + 
  geom_point(aes(n,aic.msv),col="navyblue",alpha=0.5,size=2) + 
  geom_point(aes(n,aic.ech),col="royalblue",alpha=0.5,size=2) + 
  geom_point(aes(n,aic.msh),col="skyblue",alpha=0.5,size=2) + 
  geom_point(aes(n,aic.ndvi),col="red4",alpha=0.5,size=2) + 
  geom_point(data=aic.agg,aes(n,aic.mean),size=1) + 
  geom_line(data=aic.agg,aes(n,aic.mean,group=1),size=1.2,lty=1,alpha=0.9) + 
  geom_line(data=aic.agg,aes(n,aic.ecv,group=1),size=0.5,lty=1,col="midnightblue",alpha=0.7)+ 
  geom_line(data=aic.agg,aes(n,aic.msv,group=1),size=0.5,lty=1,col="navyblue",alpha=0.7)+ 
  geom_line(data=aic.agg,aes(n,aic.ech,group=1),size=0.5,lty=1,col="royalblue",alpha=0.7)+ 
  geom_line(data=aic.agg,aes(n,aic.msh,group=1),size=0.5,lty=1,col="skyblue",alpha=0.7)+ 
  geom_line(data=aic.agg,aes(n,aic.ndvi,group=1),size=0.5,lty=1,col="red4",alpha=0.7)


# plot mean AIC by sample size
ggplot(aic.agg,aes(n,aic.mean)) + 
  theme(legend.position="top",plot.title = element_text(hjust = 0.5,size=12),plot.subtitle = element_text(size=10,hjust = 0.5),axis.text=element_text(size=12),panel.background = element_rect(fill = NA),panel.grid.major = element_line(size=0.1,colour = "grey50"),panel.border = element_rect(size=0.3, fill = NA,colour = "black")) +
  scale_x_discrete() +
  labs( x="n",y="AIC", title='AIC mean by sample size (n)') +
  geom_point(aes(n,aic.mean),size=2) + 
  geom_line(aes(n,aic.mean,group=1,lty="a"),size=1) +
  geom_point(aes(n,aic.mean.min),size=1) + 
  geom_line(aes(n,aic.mean.min,group=1,lty="b"),size=0.5) +
  geom_point(aes(n,aic.mean.max),size=1) + 
  geom_line(aes(n,aic.mean.max,group=1,lty="c"),size=0.5) +
  scale_linetype_manual(name = "", values = c("a" = 1, "b" = 3, "c" = 3), labels = c("Mean","Min", "Max"))


#########
### D-KL - by sample size
dkl.res <- data.frame(dkl.ecv=numeric(),dkl.aic.msv=numeric(),dkl.ech=numeric(),dkl.msh=numeric(),dkl.ndvi=numeric(),sid=character(),fid=integer(),n=integer())
bins=20

if(TRUE){ # run block
  h.full <- hist(field.df$ecv, breaks = seq(0,1,by=1/bins),plot=FALSE)
  h.full$counts[h.full$counts==0]<-0.01
  d.full.ecv=h.full$counts/sum(h.full$counts)
  
  h.full <- hist(field.df$msv, breaks = seq(0,1,by=1/bins),plot=FALSE)
  h.full$counts[h.full$counts==0]<-0.01
  d.full.msv=h.full$counts/sum(h.full$counts)
  
  h.full <- hist(field.df$ech, breaks = seq(0,1,by=1/bins),plot=FALSE)
  h.full$counts[h.full$counts==0]<-0.01
  d.full.ech=h.full$counts/sum(h.full$counts)
  
  h.full <- hist(field.df$msh, breaks = seq(0,1,by=1/bins),plot=FALSE)
  h.full$counts[h.full$counts==0]<-0.01
  d.full.msh=h.full$counts/sum(h.full$counts)
  
  h.full <- hist(field.df$ndvi, breaks = seq(0,1,by=1/bins),plot=FALSE)
  h.full$counts[h.full$counts==0]<-0.01
  d.full.ndvi=h.full$counts/sum(h.full$counts)
}



for (i in 1:nrow(all.fronts)){
  samp<-eval(parse(text=paste0("sample",i)))
  
  samp.ecv<-samp$ecv
  h.samp <- hist(samp.ecv, breaks = seq(0,1,by=1/bins), plot=FALSE)
  h.samp$counts[h.samp$counts==0]<-0.01
  d.samp.ecv=h.samp$counts/sum(h.samp$counts)
  
  samp.msv<-samp$msv
  h.samp <- hist(samp.msv, breaks = seq(0,1,by=1/bins), plot=FALSE)
  h.samp$counts[h.samp$counts==0]<-0.01
  d.samp.msv=h.samp$counts/sum(h.samp$counts)
  
  samp.ech<-samp$ech
  h.samp <- hist(samp.ech, breaks = seq(0,1,by=1/bins), plot=FALSE)
  h.samp$counts[h.samp$counts==0]<-0.01
  d.samp.ech=h.samp$counts/sum(h.samp$counts)
  
  samp.msh<-samp$msh
  h.samp <- hist(samp.msh, breaks = seq(0,1,by=1/bins), plot=FALSE)
  h.samp$counts[h.samp$counts==0]<-0.01
  d.samp.msh=h.samp$counts/sum(h.samp$counts)
  
  samp.ndvi<-samp$ndvi
  h.samp <- hist(samp.ndvi, breaks = seq(0,1,by=1/bins), plot=FALSE)
  h.samp$counts[h.samp$counts==0]<-0.01
  d.samp.ndvi=h.samp$counts/sum(h.samp$counts)
  
  
  kl.ecv <- signif(KL.plugin(d.full.ecv,d.samp.ecv),3)
  kl.msv <- signif(KL.plugin(d.full.msv,d.samp.msv),3)
  kl.ech <- signif(KL.plugin(d.full.ech,d.samp.ech),3)
  kl.msh <- signif(KL.plugin(d.full.msh,d.samp.msh),3)
  kl.ndvi <- signif(KL.plugin(d.full.ndvi,d.samp.ndvi),3)
  
  sid<-samp[1,"sampleID"]
  fid<-samp[1,"frontID"]
  n=nrow(samp)
  newrow<-data.frame(kl.ecv,kl.msv,kl.ech,kl.msh,kl.ndvi,sid,fid,n)
  dkl.res<-rbind(dkl.res,newrow)
}

dkl.res$kl.mean<-apply(dkl.res[1:5],1,mean)
dkl.res$n<-as.factor(dkl.res$n)
dkl.agg <-aggregate(dkl.res[,c(1:5,9)], list(n=dkl.res$n),mean)
dkl.agg.min <-aggregate(dkl.res[,c(1:5,9)], list(n=dkl.res$n),min)
dkl.agg$kl.mean.min<-dkl.agg.min$kl.mean
dkl.agg.max <-aggregate(dkl.res[,c(1:5,9)], list(n=dkl.res$n),max)
dkl.agg$kl.mean.max<-dkl.agg.max$kl.mean


#pl.dkl.all<-
ggplot(dkl.res,aes(n,kl.ecv)) + 
  theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=12),plot.subtitle = element_text(size=10,hjust = 0.5),axis.text=element_text(size=12),panel.background = element_rect(fill = NA),panel.grid.major = element_line(size=0.1,colour = "grey50"),panel.border = element_rect(size=0.3, fill = NA,colour = "black")) +
  scale_x_discrete() +
  labs( x="n",y="D-KL", title='D-KL by sample size (n)') +
  geom_point(aes(n,kl.ecv),col="purple4",alpha=0.5,size=2) + 
  geom_point(aes(n,kl.msv),col="navyblue",alpha=0.5,size=2) + 
  geom_point(aes(n,kl.ech),col="royalblue",alpha=0.5,size=2) + 
  geom_point(aes(n,kl.msh),col="seagreen4",alpha=0.5,size=2) + 
  geom_point(aes(n,kl.ndvi),col="red4",alpha=0.5,size=2) + 
  geom_point(data=dkl.agg,aes(n,kl.mean),size=4,shape=1) + 
  geom_line(data=dkl.agg,aes(n,kl.mean,group=1),size=1.2,lty=1) + 
  geom_line(data=dkl.agg,aes(n,kl.ecv,group=1),size=0.7,lty=1,col="purple4",alpha=0.7)+ 
  geom_line(data=dkl.agg,aes(n,kl.msv,group=1),size=0.7,lty=1,col="navyblue",alpha=0.7)+ 
  geom_line(data=dkl.agg,aes(n,kl.ech,group=1),size=0.7,lty=1,col="royalblue",alpha=0.7)+ 
  geom_line(data=dkl.agg,aes(n,kl.msh,group=1),size=0.7,lty=1,col="seagreen4",alpha=0.7)+ 
  geom_line(data=dkl.agg,aes(n,kl.ndvi,group=1),size=0.7,lty=1,col="red4",alpha=0.7)


# plot mean D-KL by sample size
ggplot(dkl.agg,aes(n,kl.mean)) + 
  theme(legend.position="top",plot.title = element_text(hjust = 0.5,size=12),plot.subtitle = element_text(size=10,hjust = 0.5),axis.text=element_text(size=12),panel.background = element_rect(fill = NA),panel.grid.major = element_line(size=0.1,colour = "grey50"),panel.border = element_rect(size=0.3, fill = NA,colour = "black")) +
  scale_x_discrete() +
  labs( x="n",y="D-KL", title='D-KL mean by sample size (n)') +
  geom_point(aes(n,kl.mean),size=2) + 
  geom_line(aes(n,kl.mean,group=1,lty="a"),size=1) +
  geom_point(aes(n,kl.mean.min),size=1) + 
  geom_line(aes(n,kl.mean.min,group=1,lty="b"),size=0.5) +
  geom_point(aes(n,kl.mean.max),size=1) + 
  geom_line(aes(n,kl.mean.max,group=1,lty="c"),size=0.5) +
  scale_linetype_manual(name = "", values = c("a" = 1, "b" = 3, "c" = 3), labels = c("Mean","Min", "Max"))





### Kriging Variance - by sample size
library(gstat)
library(raster)
library(rgdal)

utm <- "+proj=utm +zone=36N +datum=WGS84"
grd <- readRDS('ny_1x1_grid.RData') # grid
perimeter <- readOGR(dsn = ".", layer = "perimeter")
perimeter.df <- SpatialPolygonsDataFrame(perimeter,data=data.frame(row.names=row.names(perimeter)))
proj4string(perimeter)
perimeter <-spTransform(perimeter,CRS=utm) 
perimeter.df <-spTransform(perimeter.df,CRS=utm) 

mokv.res <- data.frame(mokv.ecv=numeric(),mokv.ech=numeric(),mokv.mean=numeric(),sid=character(),fid=integer(),n=integer())

#all.fronts.30.42<-all.fronts[as.numeric(as.character(all.fronts$n))>29,]
for (i in 1:nrow(all.fronts)){
  samp<-eval(parse(text=paste0("sample",i)))
  sid<-samp[1,"sampleID"]
  fid<-samp[1,"frontID"]
  n=nrow(samp)
  
  coordinates(samp) = ~x + y
  utm <- "+proj=utm +zone=36N +datum=WGS84"
  proj4string(samp) <- CRS(utm)
  
  oldw <- getOption("warn") 
  options(warn = -1) # hide variogram warnings
  v1 = variogram(ecv~1, samp)
  h1 = variogram(ech~1, samp)
  fitv<-fit.variogram(v1, vgm(c("Exp", "Mat", "Sph")))
  fith<-fit.variogram(h1, vgm(c("Exp", "Mat", "Sph")))
  #plot(v1,fitv,cloud=TRUE)
  options(warn = oldw)
  
  ecv.kr = krige(ecv~1, samp, grd, model = fitv,debug.level = 0)
  ech.kr = krige(ech~1, samp, grd, model = fith,debug.level = 0)
  
  r.ecv.kv <- raster(ecv.kr["var1.var"])
  r.ech.kv <- raster(ech.kr["var1.var"])
  
  # Crop using extent, rasterize polygon and finally, create poly-raster
  cr.eca.v <- crop(r.ecv.kv, extent(perimeter))
  fr.eca.v <- rasterize(perimeter.df, cr.eca.v)   
  lr.eca.v <- mask(x=cr.eca.v, mask=fr.eca.v)
  
  cr.eca.h <- crop(r.ech.kv, extent(perimeter))
  fr.eca.h <- rasterize(perimeter.df, cr.eca.h)   
  lr.eca.h <- mask(x=cr.eca.h, mask=fr.eca.h)
  
  mean.var.v = cellStats(lr.eca.v, stat='mean')
  mean.var.h = cellStats(lr.eca.h, stat='mean')
  mean.var<-mean(mean.var.v,mean.var.h)
  newrow<-data.frame(mean.var.v,mean.var.h,mean.var,sid,fid,n)
  mokv.res<-rbind(mokv.res,newrow)
  
  print(paste("mokv -",n,"points /done"))
}

#saveRDS(mokv.res, 'mokv_res_10-50a.RData')
mokv.res<-readRDS('mokv_res_10-50.RData')

mokv.res$n<-as.factor(mokv.res$n)
mokv.agg <-aggregate(mokv.res[,c(1:3)], list(n=mokv.res$n),mean)
mokv.agg.min <-aggregate(mokv.res[,c(1:3)], list(n=mokv.res$n),min)
mokv.agg$mokv.mean.min<-mokv.agg.min$mean.var
mokv.agg.max <-aggregate(mokv.res[,c(1:3)], list(n=mokv.res$n),max)
mokv.agg$mokv.mean.max<-mokv.agg.max$mean.var

#pl.mokv.all<-
ggplot(mokv.res,aes(n,mean.var.v)) + 
  theme(legend.position=c(0.8,0.8),legend.key = element_rect(colour = "transparent", fill = "white"),plot.title = element_text(hjust = 0.5,size=12),plot.subtitle = element_text(size=10,hjust = 0.5),axis.text=element_text(size=12),panel.background = element_rect(fill = NA),panel.grid.major = element_line(size=0.1,colour = "grey50"),panel.border = element_rect(size=0.3, fill = NA,colour = "black")) +
  scale_x_discrete() +
  labs( x="n",y="MOKV", title='MOKV by sample size (n)') +
  geom_point(aes(n,mean.var.v,col="a"),alpha=0.8,size=1) + 
  geom_point(aes(n,mean.var.h,col="b"),alpha=0.8,size=1) +
  geom_point(data=mokv.agg,aes(n,mean.var,col="c"),alpha=1,size=3)+
  geom_line(data=mokv.agg,aes(n,mean.var,group=1),size=1.2,lty=1) +
  scale_color_manual(name = "", values = c("a" = "navyblue", "b" = "red", "c" = "black"), labels = c("ECa V","ECa H", "Mean"))

# plot mean MOKV by sample size
ggplot(mokv.agg,aes(n,mean.var)) + 
  theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=12),plot.subtitle = element_text(size=10,hjust = 0.5),axis.text=element_text(size=12),panel.background = element_rect(fill = NA),panel.grid.major = element_line(size=0.1,colour = "grey50"),panel.border = element_rect(size=0.3, fill = NA,colour = "black")) +
  scale_x_discrete() +
  labs( x="n",y="MOKV", title='mean MOKV by sample size (n)') +
  geom_point(aes(n,mean.var),size=2) + 
  geom_line(aes(n,mean.var,group=1),size=1,lty=1) +
  geom_point(aes(n,mokv.mean.min),size=1) + 
  geom_line(aes(n,mokv.mean.min,group=1,lty="b"),size=0.5) +
  geom_point(aes(n,mokv.mean.max),size=1) + 
  geom_line(aes(n,mokv.mean.max,group=1,lty="c"),size=0.5) +
  scale_linetype_manual(name = "", values = c("a" = 1, "b" = 3, "c" = 3), labels = c("Mean","Min", "Max"))



### merge
all.indices<-cbind(aic.res,dkl.res[c(1:5,9)],mokv.res[1:3])
all.indices$aic.n <- normal(all.indices$aic.mean)
all.indices$dkl.n <- normal(all.indices$kl.mean)
all.indices$nokv.n <- normal(all.indices$mean.var)
#all.indices$sum.n <- y.aic+y.dkl+y.var


all.indices.agg<-cbind(aic.agg,dkl.agg[2:9],mokv.agg[2:6])
#saveRDS(all.indices,"all.indices.RData")
#saveRDS(all.indices.agg,"all.indices.agg.RData")

#all.indices<-readRDS("all.indices.RData")
#all.indices.agg<-readRDS("all.indices.agg.RData")


indices.theme <-theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=12),plot.subtitle = element_text(size=10,hjust = 0.5),axis.text=element_text(size=12),panel.background = element_rect(fill = NA),panel.grid.major = element_line(size=0.1,colour = "grey50"),panel.border = element_rect(size=0.3, fill = NA,colour = "black")) 



pl.aic<-
  ggplot(all.indices) + geom_point(aes(n,aic.mean,col="a"),size=0.5) + indices.theme + scale_x_discrete() + labs(x="n",y="AIC", title='AIC mean by sample size (n)') + 
  geom_line(data=all.indices.agg,aes(n,aic.mean.max,group=1,lty="a"),size=0.5) + 
  geom_line(data=all.indices.agg,aes(n,aic.mean,group=1,lty="b"),size=0.7) +
  geom_line(data=all.indices.agg,aes(n,aic.mean.min,group=1,lty="c"),size=0.5) + 
  scale_linetype_manual(name = "", values = c("a" = 3, "b" = 1, "c" = 3), labels = c("max","mean", "min")) + scale_color_manual(name = "", values = c("a" = "midnightblue", "b" = "purple4", "c" = "steelblue"), labels = c("AIC")) +
  theme(legend.position=c(0.94,0.66),panel.grid.minor = element_line(size=0.1,colour = "grey50"),panel.grid.major = element_line(size=0.1,colour = "grey50"),legend.background = element_blank(), legend.title = element_blank(),legend.box.background = element_rect(colour = "grey50",size=0.2),legend.key = element_rect(colour = "transparent", fill = "white"),legend.text = element_text(size=11)) + 
    # + theme(legend.position="none") + 
    guides(color=FALSE) 

pl.dkl<-
  ggplot(all.indices) + geom_point(aes(n,kl.mean,col="a"),size=1) + indices.theme + scale_x_discrete() + labs(x="n",y="D-KL", title='D-KL mean by sample size (n)') + 
  geom_line(data=all.indices.agg,aes(n,kl.mean.max,group=1,lty="a"),size=0.5) + 
  geom_line(data=all.indices.agg,aes(n,kl.mean,group=1,lty="b"),size=0.7) +
  geom_line(data=all.indices.agg,aes(n,kl.mean.min,group=1,lty="c"),size=0.5) + 
  scale_linetype_manual(name = "", values = c("a" = 3, "b" = 1, "c" = 3), labels = c("max","mean", "min")) + scale_color_manual(name = "", values = c("a" = "midnightblue", "b" = "purple4", "c" = "steelblue"), labels = c("D-KL")) +
  theme(legend.position=c(0.8,0.82),panel.grid.minor = element_line(size=0.1,colour = "grey50"),panel.grid.major = element_line(size=0.1,colour = "grey50"),legend.background = element_blank(), legend.title = element_blank(),legend.box.background = element_rect(colour = "grey50",size=0.2),legend.key = element_rect(colour = "transparent", fill = "white"),legend.text = element_text(size=11)) + theme(legend.position="none")

pl.mokv<-
  ggplot(all.indices) + geom_point(aes(n,mean.var,col="a"),size=1) + indices.theme + scale_x_discrete() + labs(x="n",y="MOKV", title='MOKV mean by sample size (n)') + 
    geom_line(data=all.indices.agg,aes(n,mokv.mean.max,group=1,lty="a"),size=0.5) + 
    geom_line(data=all.indices.agg,aes(n,mean.var,group=1,lty="b"),size=0.7) +
    geom_line(data=all.indices.agg,aes(n,mokv.mean.min,group=1,lty="c"),size=0.5) + 
    scale_linetype_manual(name = "", values = c("a" = 3, "b" = 1, "c" = 3), labels = c("max","mean", "min")) + scale_color_manual(name = "", values = c("a" = "midnightblue", "b" = "purple4", "c" = "steelblue"), labels = c("MOKV")) +
    theme(legend.position=c(0.8,0.82),panel.grid.minor = element_line(size=0.1,colour = "grey50"),panel.grid.major = element_line(size=0.1,colour = "grey50"),legend.background = element_blank(), legend.title = element_blank(),legend.box.background = element_rect(colour = "grey50",size=0.2),legend.key = element_rect(colour = "transparent", fill = "white"),legend.text = element_text(size=11)) + theme(legend.position="none")

  multiplot(pl.aic,pl.dkl,pl.mokv) 
  multiplot(pl.aic,pl.dkl,pl.mokv,pl.grad) 

multiplot(pl.aic.all,pl.dkl.all,pl.mokv.all)

multiplot(pl.aic,pl.dkl,pl.mokv,pl.grad.1st.2nd,pl.grad.2nd)




### normalize
x <- npoints

# Individual schemes
# y.aic <- normal(all.indices$aic.mean)
# y.dkl <- normal(all.indices$kl.mean)
# y.var <- normal(all.indices$mean.var)
# y.sum <- y.aic+y.dkl+y.var

# f.aic <- gradient(y.aic,x)
# f.dkl <- gradient(y.dkl,x)
# f.var <- gradient(y.var,x)
# f.sum<-f.aic+f.dkl+f.var
# ff.sum <- gradient(f.sum,x)

# grads<-data.frame(n=as.factor(npoints),f.aic=f.aic,f.dkl=f.dkl,f.var=f.var,f.sum=f.sum,ff.sum=ff.sum,y.sum=y.sum)
# 
# # Individual schemes - top 3
# library(dplyr)
# grads.top <- grads %>%
#   group_by(n) %>%
#   arrange(desc(y.sum)) %>%
#   slice(1:3)
# 
# grads.top.agg <-aggregate(grads.top[,c(1:3)], list(n=grads.top$n),mean)
# 
# 
# f.aic <- gradient(y.aic,x)
# f.dkl <- gradient(y.dkl,x)
# f.var <- gradient(y.var,x)
# f.sum<-f.aic+f.dkl+f.var
# ff.sum <- gradient(f.sum,x)


#<-data.frame(n=as.factor(npoints),f.aic=f.aic,f.dkl=f.dkl,f.var=f.var,f.sum=f.sum,ff.sum=ff.sum,y.sum=y.sum)


# Aggregated by sample size
y.aic <- normal(all.indices.agg$aic.mean)
y.dkl <- normal(all.indices.agg$kl.mean)
y.var <- normal(all.indices.agg$mean.var)
y.sum <- y.aic+y.dkl+y.var

f.aic <- gradient(y.aic,x)
f.dkl <- gradient(y.dkl,x)
f.var <- gradient(y.var,x)
f.sum<-f.aic+f.dkl+f.var
ff.sum <- gradient(f.sum,x)

all.indices.agg$y.sum<-y.sum
all.indices.agg$f.sum<-f.sum
all.indices.agg$ff.sum<-ff.sum


### get best value
y.aic.min <- normal(all.indices.agg$aic.mean.min)
y.dkl.min <- normal(all.indices.agg$kl.mean.min)
y.var.min <- normal(all.indices.agg$mokv.mean.min)
y.sum.min <- y.aic.min+y.dkl.min+y.var.min

y.aic.max <- normal(all.indices.agg$aic.mean.max)
y.dkl.max <- normal(all.indices.agg$kl.mean.max)
y.var.max <- normal(all.indices.agg$mokv.mean.max)
y.sum.max <- y.aic.max+y.dkl.max+y.var.max

f.aic.min <- gradient(y.aic.min,x)
f.dkl.min <- gradient(y.dkl.min,x)
f.var.min <- gradient(y.var.min,x)
f.sum.min<-f.aic.min+f.dkl.min+f.var.min
ff.sum.min <- gradient(f.sum.min,x)

all.indices.agg$y.sum.min<-y.sum.min
all.indices.agg$y.sum.max<-y.sum.max
all.indices.agg$f.sum.min<-f.sum.min
all.indices.agg$ff.sum.min<-ff.sum.min


#saveRDS(all.indices.agg,"all.indices.agg.RData")

#grads.agg<-data.frame(n=as.factor(npoints),f.aic=f.aic,f.dkl=f.dkl,f.var=f.var,f.sum=f.sum,ff.sum=ff.sum,y.sum=y.sum)

pl.index<-
ggplot(all.indices.agg) + 
  geom_point(aes(n,y.sum,col="a"),size=1.3) + 
  geom_line(aes(n,y.sum,group=1,col="a"),size=1.2,lty=1) + 
  #geom_point(aes(n,y.sum.min,col="b"),size=1) + 
  #geom_line(aes(n,y.sum.min,group=1,col="b"),size=1,lty=3) + 
  #geom_point(aes(n,y.sum.max,col="c"),size=1) + 
  #geom_line(aes(n,y.sum.max,group=1,col="c"),size=1,lty=3) + 
  
  scale_color_manual(name = "", values = c("a" = "midnightblue", "b" = "purple4", "c" = "steelblue"), labels = c("mean","min", "f(x)\'\'")) +
  indices.theme + scale_x_discrete() + 
  labs(x="n",y="", title='Information index by sample size', subtitle='f(x) = n(aic) + n(dkl) + n(mokv)') +
  theme(legend.position=c(0.94,0.75),panel.grid.minor = element_line(size=0.2,colour = "grey50"),panel.grid.major = element_line(size=0.2,colour = "grey50"),legend.background = element_blank(), legend.title = element_blank(),legend.box.background = element_rect(colour = "grey50",size=0.2),legend.key = element_rect(colour = "transparent", fill = "white"),legend.text = element_text(size=11)) + theme(legend.position="none")

multiplot(pl.aic,pl.dkl,pl.mokv,pl.grad) 

#pl.grad.agg<-
  ggplot(all.indices.agg) + 
    geom_point(aes(n,y.sum,col="a"),size=1.3) + 
    geom_point(aes(n,f.sum,col="b"),size=1) + 
    geom_point(aes(n,ff.sum,col="c"),size=1) + 
    geom_line(aes(n,y.sum,group=1,col="a"),size=1.2,lty=1) + 
    geom_line(aes(n,f.sum,group=1,col="b"),size=1,lty=1) + 
    geom_line(aes(n,ff.sum,group=1,col="c"),size=1,lty=1) + 
    scale_color_manual(name = "", values = c("a" = "midnightblue", "b" = "purple4", "c" = "steelblue"), labels = c("f(x)","f(x)\'", "f(x)\'\'")) +
    geom_hline(yintercept = 0) +
    indices.theme + scale_x_discrete() + 
    labs(x="n",y="", title='Aggregated indices - mean', subtitle='f(x) = n(aic) + n(dkl) + n(mokv)') +
    theme(legend.position=c(0.8,0.82),panel.grid.minor = element_line(size=0.1,colour = "grey50"),panel.grid.major = element_line(size=0.1,colour = "grey50"),legend.background = element_blank(), legend.title = element_blank(),legend.box.background = element_rect(colour = "grey50",size=0.2),legend.key = element_rect(colour = "transparent", fill = "white"),legend.text = element_text(size=11))# + theme(legend.position="none")
  
pl.grad.agg.2nd<-
  ggplot(all.indices.agg) + 
    geom_point(aes(n,ff.sum,col="c"),size=1) + 
    geom_line(aes(n,ff.sum,group=1,col="c"),size=1.1,lty=1) + 
    scale_color_manual(name = "", values = c("a" = "midnightblue", "b" = "purple4", "c" = "steelblue"), labels = c("f(x)\'\'")) +
    geom_hline(yintercept = 0) +
    indices.theme + scale_x_discrete() + 
    labs(x="n",y="", title='Information acceleration') +
    theme(legend.position=c(0.93,0.8),panel.grid.minor = element_line(size=0.2,colour = "grey50"),panel.grid.major = element_line(size=0.2,colour = "grey50"),legend.background = element_blank(), legend.title = element_blank(),legend.box.background = element_rect(colour = "grey50",size=0.2),legend.key = element_rect(colour = "transparent", fill = "white"),legend.text = element_text(size=11)) + ylim(-0.1,0.1)

  

multiplot(pl.index+theme(plot.margin = unit(c(0.2,1,0,1), "cm")),pl.grad.agg.2nd+theme(plot.margin = unit(c(0.3,1,0,0.3), "cm"))) 
multiplot(pl.grad.agg,pl.grad.agg.2nd)


 # pl.grad.agg.top<-  
 #  ggplot(all.indices.agg) + 
 #    geom_point(aes(n,y.sum.min,col="a"),size=1) + 
 #    geom_point(aes(n,f.sum.min,col="b"),size=1) + 
 #    geom_point(aes(n,ff.sum.min,col="c"),size=1) + 
 #    geom_line(aes(n,y.sum.min,group=1,col="a"),size=1.1,lty=1) + 
 #    geom_line(aes(n,f.sum.min,group=1,col="b"),size=1.1,lty=1) + 
 #    geom_line(aes(n,ff.sum.min,group=1,col="c"),size=1.1,lty=1) + 
 #    scale_color_manual(name = "", values = c("a" = "midnightblue", "b" = "purple4", "c" = "steelblue"), labels = c("f(x)","f(x)\'", "f(x)\'\'")) +
 #    geom_hline(yintercept = 0) +
 #    indices.theme + scale_x_discrete() + 
 #    labs(x="n",y="", title='Aggregated indices - by top value', subtitle='f(x) = min(aic.min) + min(dkl) + min(mokv)') +
 #    theme(legend.position=c(0.8,0.82),panel.grid.minor = element_line(size=0.2,colour = "grey50"),panel.grid.major = element_line(size=0.2,colour = "grey50"),legend.background = element_blank(), legend.title = element_blank(),legend.box.background = element_rect(colour = "grey50",size=0.2),legend.key = element_rect(colour = "transparent", fill = "white"),legend.text = element_text(size=11))# + theme(legend.position="none")  
  
  

# pl.grad.agg.top<-
# ggplot(grads.agg) + 
#   geom_point(aes(n,y.sum,color="a"),size=2) + geom_line(aes(n,y.sum,group=1,color="a"),size=1,lty=1) + 
#   geom_point(aes(n,f.sum,color="b"),size=2) + geom_line(aes(n,f.sum,group=1,color="b"),size=1,lty=1) +
#   geom_point(aes(n,ff.sum,color="c"),size=2)+geom_line(aes(n,ff.sum,group=1,color="c"),size=1,lty=1) +
#   scale_color_manual(name = "", values = c("a" = "midnightblue", "b" = "purple4", "c" = "steelblue"), labels = c("f(x)","f(x)\'", "f(x)\'\'")) +
#   labs(x="n",y="", title='Aggregated indices (by top value)',subtitle="f(x) = min(aic) + min(dkl) + min(mokv)") + geom_hline(yintercept = 0) +
#   indices.theme + theme(legend.position=c(0.9,0.7),panel.grid.minor = element_line(size=0.2,colour = "grey50"),panel.grid.major = element_line(size=0.2,colour = "grey50"),legend.background = element_blank(), legend.title = element_blank(),legend.box.background = element_rect(colour = "black"),        legend.key = element_rect(colour = "transparent", fill = "white"),legend.text = element_text(size=12)) 


pl.grad.agg.top.2nd<-
ggplot(all.indices.agg) + 
    geom_point(aes(n,ff.sum.min,col="c"),size=1) + 
    geom_line(aes(n,ff.sum.min,group=1,col="c"),size=1.1,lty=1) + 
    scale_color_manual(name = "", values = c("a" = "midnightblue", "b" = "purple4", "c" = "steelblue"), labels = c("f(x)\'\'")) +
    geom_hline(yintercept = 0) +
    indices.theme + scale_x_discrete() + 
    labs(x="n",y="", title='Information acceleration (by top value)', subtitle='f(x) = n(aic) + n(dkl) + n(mokv)') +
    theme(legend.position=c(0.8,0.82),panel.grid.minor = element_line(size=0.2,colour = "grey50"),panel.grid.major = element_line(size=0.2,colour = "grey50"),legend.background = element_blank(), legend.title = element_blank(),legend.box.background = element_rect(colour = "grey50",size=0.2),legend.key = element_rect(colour = "transparent", fill = "white"),legend.text = element_text(size=11))

multiplot(pl.grad.agg.top,pl.grad.agg.top.2nd)

#pl.grad.2nd.agg<-
# ggplot(grads.agg) + 
#   geom_point(aes(n,ff.sum,color="c"),size=2)+geom_line(aes(n,ff.sum,group=1,color="c"),size=1,lty=1) +
#   scale_color_manual(name = "", values = c("a" = "midnightblue", "b" = "purple4", "c" = "steelblue"), labels = c("f(x)\'\'"))+
#   labs(x="n",y="f(x)\'\'", title="f(x) = aic.mean [0..1] + dkl.mean [0..1] + mokv.mean [0..1]") + geom_hline(yintercept = 0) +
#   indices.theme + theme(legend.position=c(0.9,0.7),panel.grid.minor = element_line(size=0.1,colour = "grey50"),panel.grid.major = element_line(size=0.1,colour = "grey50"),legend.background = element_blank(), legend.title = element_blank(),legend.box.background = element_rect(colour = "black"),        legend.key = element_rect(colour = "transparent", fill = "white"),legend.text = element_text(size=12)) 


multiplot(pl.grad.1st.2nd.agg,pl.grad.2nd.agg)

multiplot(pl.aic,pl.dkl,pl.mokv,pl.grad.1st.2nd.agg,pl.grad.2nd.agg)









###################### Individual schemes evaluation

#################

### AIC - Individual scheme

front<-front.26
aic.res <- data.frame(aic.ecv=numeric(),aic.msv=numeric(),aic.ech=numeric(),aic.msh=numeric(),aic.ndvi=numeric(),sid=character(),fid=integer())


for (i in 1:nrow(front)){
  samp<-eval(parse(text=paste0("sample",front[i,"frontID"])))
  
  lm <- lm(ecv ~ . , samp[d.cols])
  aic.ecv <- AIC(lm)
  
  lm <- lm(msv ~ . , samp[d.cols])
  aic.msv <- AIC(lm)
  
  lm <- lm(ech ~ . , samp[d.cols])
  aic.ech <- AIC(lm)
  
  lm <- lm(ecv ~ . , samp[d.cols])
  aic.ecv <- AIC(lm)
  
  lm <- lm(msh ~ . , samp[d.cols])
  aic.msh <- AIC(lm)
  
  lm <- lm(ndvi ~ . , samp[d.cols])
  aic.ndvi <- AIC(lm)
  
  sid<-samp[1,"sampleID"]
  fid<-samp[1,"frontID"]
  newrow<-data.frame(aic.ecv,aic.msv,aic.ech,aic.msh,aic.ndvi,sid,fid)
  aic.res<-rbind(aic.res,newrow)
}

aic.res
aic.res$mean<-apply(aic.res[1:5],1,mean)

#pl.front

ggplot(aic.res,aes(factor(fid),aic.ecv)) + 
  theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=12),plot.subtitle = element_text(size=10,hjust = 0.5),axis.text=element_text(size=12),panel.background = element_rect(fill = NA),panel.grid.major = element_line(size=0.1,colour = "grey50"),panel.grid.minor = element_line(size=0.1,colour = "grey50"),panel.border = element_rect(size=0.2, fill = NA)) + labs(x="") +   scale_x_discrete() + labs( x="Sample #",y="AIC", title=paste('AIC by sample -',nsamp,"points")) + geom_point(aes(factor(fid),aic.ecv),col="blue") + geom_point(aes(factor(fid),aic.msv),col="orange") + geom_point(aes(factor(fid),aic.ech),col="green") + geom_point(aes(factor(fid),aic.msh),col="purple") + geom_point(aes(factor(fid),aic.ndvi),col="red") + geom_point(aes(factor(fid),mean),size=2) + geom_line(aes(factor(fid),mean,group=1),size=1)






### D-KL - individual scheme
dkl.res <- data.frame(dkl.ecv=numeric(),dkl.aic.msv=numeric(),dkl.ech=numeric(),dkl.msh=numeric(),dkl.ndvi=numeric(),sid=character(),id=integer())
bins=20

if(TRUE){ # run block
  h.full <- hist(field.df$ecv, breaks = seq(0,1,by=1/bins),plot=FALSE)
  h.full$counts[h.full$counts==0]<-0.01
  d.full.ecv=h.full$counts/sum(h.full$counts)
  
  h.full <- hist(field.df$msv, breaks = seq(0,1,by=1/bins),plot=FALSE)
  h.full$counts[h.full$counts==0]<-0.01
  d.full.msv=h.full$counts/sum(h.full$counts)
  
  h.full <- hist(field.df$ech, breaks = seq(0,1,by=1/bins),plot=FALSE)
  h.full$counts[h.full$counts==0]<-0.01
  d.full.ech=h.full$counts/sum(h.full$counts)
  
  h.full <- hist(field.df$msh, breaks = seq(0,1,by=1/bins),plot=FALSE)
  h.full$counts[h.full$counts==0]<-0.01
  d.full.msh=h.full$counts/sum(h.full$counts)
  
  h.full <- hist(field.df$ndvi, breaks = seq(0,1,by=1/bins),plot=FALSE)
  h.full$counts[h.full$counts==0]<-0.01
  d.full.ndvi=h.full$counts/sum(h.full$counts)
}


for (i in 1:nrow(front)){
  samp<-eval(parse(text=paste0("sample",front[i,"frontID"])))
  
  samp.ecv<-samp$ecv
  h.samp <- hist(samp.ecv, breaks = seq(0,1,by=1/bins), plot=FALSE)
  h.samp$counts[h.samp$counts==0]<-0.01
  d.samp.ecv=h.samp$counts/sum(h.samp$counts)
  
  samp.msv<-samp$msv
  h.samp <- hist(samp.msv, breaks = seq(0,1,by=1/bins), plot=FALSE)
  h.samp$counts[h.samp$counts==0]<-0.01
  d.samp.msv=h.samp$counts/sum(h.samp$counts)
  
  samp.ech<-samp$ech
  h.samp <- hist(samp.ech, breaks = seq(0,1,by=1/bins), plot=FALSE)
  h.samp$counts[h.samp$counts==0]<-0.01
  d.samp.ech=h.samp$counts/sum(h.samp$counts)
  
  samp.msh<-samp$msh
  h.samp <- hist(samp.msh, breaks = seq(0,1,by=1/bins), plot=FALSE)
  h.samp$counts[h.samp$counts==0]<-0.01
  d.samp.msh=h.samp$counts/sum(h.samp$counts)
  
  samp.ndvi<-samp$ndvi
  h.samp <- hist(samp.ndvi, breaks = seq(0,1,by=1/bins), plot=FALSE)
  h.samp$counts[h.samp$counts==0]<-0.01
  d.samp.ndvi=h.samp$counts/sum(h.samp$counts)
  
  
  kl.ecv <- signif(KL.plugin(d.full.ecv,d.samp.ecv),3)
  kl.msv <- signif(KL.plugin(d.full.msv,d.samp.msv),3)
  kl.ech <- signif(KL.plugin(d.full.ech,d.samp.ech),3)
  kl.msh <- signif(KL.plugin(d.full.msh,d.samp.msh),3)
  kl.ndvi <- signif(KL.plugin(d.full.ndvi,d.samp.ndvi),3)
  
  sid<-samp[1,"sampleID"]
  fid<-samp[1,"frontID"]
  newrow<-data.frame(kl.ecv,kl.msv,kl.ech,kl.msh,kl.ndvi,sid,fid)
  dkl.res<-rbind(dkl.res,newrow)
}

dkl.res$mean<-apply(dkl.res[1:5],1,mean)

#pl.front

ggplot(dkl.res,aes(factor(fid),kl.ecv)) + 
  theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=12),plot.subtitle = element_text(size=10,hjust = 0.5),axis.text=element_text(size=12),panel.background = element_rect(fill = NA),panel.grid.major = element_line(size=0.1,colour = "grey50"),panel.grid.minor = element_line(size=0.1,colour = "grey50"),panel.border = element_rect(size=0.2, fill = NA)) + labs(x="") +   scale_x_discrete() + labs( x="Sample #",y="D-KL", title=paste('D-KL by sample -',nsamp,"points")) + geom_point(aes(factor(fid),kl.ecv),col="blue") + geom_point(aes(factor(fid),kl.msv),col="darkblue") + geom_point(aes(factor(fid),kl.ech),col="green") + geom_point(aes(factor(fid),kl.msh),col="darkgreen") + geom_point(aes(factor(fid),kl.ndvi),col="red") + geom_point(aes(factor(fid),mean),size=2) + geom_line(aes(factor(fid),mean,group=1),size=1)



### MOKV - individual schemes

mokv.scheme<-mokv.res[mokv.res$n==26,]

ggplot(mokv.scheme,aes(factor(fid),mean.var)) + 
  theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=12),plot.subtitle = element_text(size=10,hjust = 0.5),axis.text=element_text(size=12),panel.background = element_rect(fill = NA),panel.grid.major = element_line(size=0.1,colour = "grey50"),panel.border = element_rect(size=0.3, fill = NA,colour = "black")) +
  scale_x_discrete() +
  labs( x="Sample #",y="MOKV", title='MOKV by sample -',nsamp,'points') +
  geom_point(aes(factor(fid),mean.var.v),col="navyblue",alpha=0.8,size=1) + 
  geom_point(aes(factor(fid),mean.var.h),col="red",alpha=0.8,size=1) +
  geom_point(aes(factor(fid),mean.var),col="black",alpha=0.8,size=1.1)+ 
  geom_line(aes(factor(fid),mean.var,group=1),size=1)
  #geom_point(data=mokv.agg,aes(factor(fid),mean.var),col="black",alpha=1,size=3)+
  #geom_line(data=mokv.agg,aes(factor(fid),mean.var,group=1),size=1.2,lty=1) 


### mokv by scheme
mokv.res <- data.frame(mokv.ecv=numeric(),mokv.ech=numeric(),mokv.mean=numeric(),sid=character(),fid=integer(),n=integer())


for (i in 1:nrow(front)){
  samp<-eval(parse(text=paste0("sample",front[i,"frontID"])))
  sid<-samp[1,"sampleID"]
  fid<-samp[1,"frontID"]
  n=nrow(samp)
  
  coordinates(samp) = ~x + y
  utm <- "+proj=utm +zone=36N +datum=WGS84"
  proj4string(samp) <- CRS(utm)
  
  oldw <- getOption("warn") 
  options(warn = -1) # hide variogram warnings
  v1 = variogram(ecv~1, samp)
  h1 = variogram(ech~1, samp)
  fitv<-fit.variogram(v1, vgm(c("Exp", "Mat", "Sph")))
  fith<-fit.variogram(h1, vgm(c("Exp", "Mat", "Sph")))
  #plot(v1,fitv,cloud=TRUE)
  options(warn = oldw)
  
  ecv.kr = krige(ecv~1, samp, grd, model = fitv,debug.level = 0)
  ech.kr = krige(ech~1, samp, grd, model = fith,debug.level = 0)
  
  r.ecv.kv <- raster(ecv.kr["var1.var"])
  r.ech.kv <- raster(ech.kr["var1.var"])
  
  # Crop using extent, rasterize polygon and finally, create poly-raster
  cr.eca.v <- crop(r.ecv.kv, extent(perimeter))
  fr.eca.v <- rasterize(perimeter.df, cr.eca.v)   
  lr.eca.v <- mask(x=cr.eca.v, mask=fr.eca.v)
  
  cr.eca.h <- crop(r.ech.kv, extent(perimeter))
  fr.eca.h <- rasterize(perimeter.df, cr.eca.h)   
  lr.eca.h <- mask(x=cr.eca.h, mask=fr.eca.h)
  
  mean.var.v = cellStats(lr.eca.v, stat='mean')
  mean.var.h = cellStats(lr.eca.h, stat='mean')
  mean.var<-mean(mean.var.v,mean.var.h)
  newrow<-data.frame(mean.var.v,mean.var.h,mean.var,sid,fid,n)
  mokv.res<-rbind(mokv.res,newrow)
  
  print(paste("mokv -",n,"points /done"))
}

#saveRDS(mokv.res, 'mokv_res_10-50a.RData')
mokv.res<-readRDS('mokv_res_10-50a.RData')

mokv.res$n<-as.factor(mokv.res$n)

mokv.scheme<-mokv.res[mokv.res$n==26,]

ggplot(mokv.scheme,aes(factor(fid),mean.var)) + 
  theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=12),plot.subtitle = element_text(size=10,hjust = 0.5),axis.text=element_text(size=12),panel.background = element_rect(fill = NA),panel.grid.major = element_line(size=0.1,colour = "grey50"),panel.grid.minor = element_line(size=0.1,colour = "grey50"),panel.border = element_rect(size=0.2, fill = NA)) + labs(x="") +   
  scale_x_discrete() + 
  labs( x="Sample #",y="MOKV", title='MOKV by sample -',nsamp,'points') +
  geom_point(aes(factor(fid),mean.var.v),col="navyblue",alpha=0.8,size=1) + 
  geom_point(aes(factor(fid),mean.var.h),col="red",alpha=0.8,size=1) +
  #geom_point(aes(factor(fid),mean.var),col="black",alpha=0.8,size=1.1)+ 
  geom_line(aes(factor(fid),mean.var,group=1),size=1)


