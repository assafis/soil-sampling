
ranking<-function(obj){
  Nind = nrow(obj)
  RankMOV = rep(0, Nind)
  
  for (i in 1:Nind){
    for (j in 1:Nind){
      if (i!=j){
        RankMOV[i] = RankMOV[i] + rankplt(obj[j,], obj[i,]);
      }
    }
  }
  return(RankMOV)
}

rankplt <- function(ObjV1, ObjV2){
  Prefer = all(ObjV1 <= ObjV2, any(ObjV1 < ObjV2))
  return(Prefer)
}
