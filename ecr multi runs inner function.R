library(ecr)
source('ecr multi internal-obj.R')  # obj functions
source('spatial_helper_functions.R')


pointsReplace = makeMutator(
  mutator = function(i_sampled,n_data){
    i_unsampled <- setdiff(1:n_data, i_sampled)
    idx_removed <- sample(1:length(i_sampled), size = 1, replace = FALSE)
    idx_added <- sample(1:length(i_unsampled), size = 1, replace = FALSE)
    i_sampled <- i_sampled[-idx_removed]
    i_sampled <- c(i_sampled, i_unsampled[idx_added])
    
    return(i_sampled)
  },
  supported = "custom")




inner <- function(k,mu,lambda,max.iter,nsamp,d.cols,theta){
  
  fitness.fun = function(pop,x,d.cols,continuous_strata,cor_mat,N,theta){
    c(
      #.lhs_obj(pop,x,d.cols,continuous_strata,cor_mat,N,theta),
      #.dist_obj(pop,x,d.cols,continuous_strata,cor_mat,N,theta)
      .SPdiv.xy_obj(pop,x,d.cols,continuous_strata,cor_mat,N,theta),
      .SPdiv.anc_obj(pop,x,d.cols,continuous_strata,cor_mat,N,theta)
      
      #.dkl_obj(pop,x,d.cols,continuous_strata,N)
      #.rmse_obj(pop,x,d.cols,continuous_strata,N)
      #.kvar_obj(pop,x,d.cols,continuous_strata,N)
      #.chisq_obj(pop,x,d.cols,continuous_strata,N)
      #.ks_obj(pop,x,d.cols,continuous_strata,N)    
    )  
  }
  #fitness.names = c("cLHS","maxMinDist")
  fitness.names = c("SP_xy","SP_ancillary")
  
  field_df <- field.xy.n<- readRDS("field_df.RData") # the data frame
  field.xy.n$x<-normal(field.xy.n$x)  # normalize for SP
  field.xy.n$y<-normal(field.xy.n$y)
  
  #x <- subset(field_df,feasible==1) # utm coords
  x <- subset(field.xy.n,feasible==1) # SP - normalized x,y
  n_data <- nrow(x)
  #d.cols = c(3:7) # ancillary data columns 
  cor_mat <- cor(x[d.cols], use = "complete.obs") # Data correlation - for clhs
  fronts<-data.frame(double(),double(),integer())
  colnames(fronts)<-c(fitness.names,"N")
  hv<-data.frame()[1:(MAX.ITER+1), ]
  popBest<-list()
  
  
  print(paste("Start ~~~",Sys.time()))
  
  N=nsamp
  
  # Edge of the strata - for clhs
  continuous_strata <- apply(x[d.cols], 2, function(x) { quantile(x, probs = seq(0, 1, length.out = N + 1), na.rm = TRUE) })
  
  ctrl = initECRControl(fitness.fun, n.objectives = 2L, minimize = c(TRUE,TRUE))
  ctrl = registerECROperator(ctrl, "mutate", pointsReplace)
  ctrl = registerECROperator(ctrl, "selectForSurvival", selNondom)
  
  population <- gen(sample(1:n_data, N, replace = FALSE),MU)
  fitness = evaluateFitness(ctrl,population,x,d.cols,continuous_strata,cor_mat,N,theta) 
  
  
  ref.point = c(1500,500)
  logger = initLogger(ctrl,
                      log.stats = list(fitness = list("HV"=list(
                        fun = computeHV,
                        pars = list(ref.point = ref.point)))), 
                      log.pop = TRUE,
                      init.size = MAX.ITER+1)
  
  updateLogger(logger,population,fitness=fitness,n.evals = MU)
  
  for(i in seq_len(MAX.ITER)){
    offspring = ecr::mutate(ctrl,population,p.mut = 0.8,n_data=n_data)
    fitness.o = evaluateFitness(ctrl,population,x,d.cols,continuous_strata,cor_mat,N,theta) 
    sel = replaceMuPlusLambda(ctrl, population,offspring,fitness,fitness.o)
    population=sel$population
    fitness=sel$fitness
    updateLogger(logger,population,fitness,n.evals = MU)
  }
  
  stats = getStatistics(logger)
  pl.stats = plotStatistics(stats) + theme_bw()
  pl.front = plotFront(fitness,obj.names = fitness.names) + theme_bw()
  poplog = getPopulations(logger)
  
  #print(pl.stats)
  
  hv<-cbind(hv,stats)
  denom=10 # compact coefficient
  n.rows=ceiling(nrow(hv)/denom)
  hv<-hv[seq(1, n.rows*denom, denom),]
  
  popBest<-poplog[[MAX.ITER]]
  
  
  
  print(paste("End run ---",Sys.time()))
  run<-list(hv,popBest)
  return(run)
}



