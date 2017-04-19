#' Simulation of patch-based movement trajectory
#'
#' Simulate a movement trajectory with user defined number of patches and interpatch movement
#' @param type whether movement within patches should be based on a 2states process (from package moveHMM) or a Bivariate Ornstein-Uhlenbeck process (OU) (from package adehabitatLT)
#' @param npatches Number of patches, default=5
#' @param ratio Ratio (in percent) of locations associated to interpatch movement, default=5
#' @param nswitch Number of switch/depart from patches, default=150
#' @param ncore Number of locations within a patch per visit, default=200
#' @param spacecore Minimum distance between center of patches, default=200
#' @param seq_visit Specify the sequence of visit among patches, default is random sequence
#' @param stepDist Distribution for step length if 2states specified in type, see simData of moveHMM package
#' @param angleDist Distribution for turn angle if 2states specified in type, see simData of moveHMM package
#' @param stepPar Parameters for step length distribution if 2states specified in type, see simData of moveHMM package
#' @param anglePar Parameters for turn angle distribution if 2states specified in type, see simData of moveHMM package
#' @param s Parameters for the OU process, see simm.mou of adehabitatLT package
#' @param grph Whether a graph of the trajectory should be produced, default=F
#' @keywords traj2adj adj2stack
#' @return A ltraj (adehabitatLT) object
#' @export
#' @examples
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' traj2<-sim_mov(type="2states", npatches=2, grph=T)


sim_mov<-function(type=c("2states", "OU"), npatches=5, ratio=5, nswitch=150, ncore=200,spacecore=200, seq_visit=sample(1:npatches, nswitch, replace=T), 
                  stepDist= "gamma", angleDist = "vm",  stepPar = c(0.5,3,1,5), anglePar = c(pi,0,0.5,2), s=diag(40,2), grph=F) {
  if(require("adehabitatLT") & require("moveHMM")){
  } else {
    print("trying to install R2jags")
    install.packages(c("adehabitatLT", "moveHMM"))
    if(require("adehabitatLT") & require("moveHMM")){
      print("Packages installed and loaded")
    } else {
      stop("Could not install required packages (moveHMM or adehabitatLT")
    }
  }
  coordx<-sample(seq(0,20,2), npatches, replace=F)*spacecore
  coordy<-sample(seq(0,20,2), npatches, replace=F)*spacecore
  nmig=ncore/ratio
  out<-data.frame()
  for (i in 1:(nswitch-1)){
    
    if(type=="2states") {
      core<-simData(nbAnimals=1,nbStates=2,stepDist=stepDist,angleDist=angleDist,stepPar=stepPar, anglePar=anglePar,zeroInflation=F,obsPerAnimal=ncore)
      corex<-core$x+coordx[seq_visit[i]]
      corey<-core$y+coordy[seq_visit[i]]
      Corri1<-rep(2, ncore)
    }
    
    if(type=="OU") {
      core<-simm.mou(date=1:ncore, b=c(coordx[seq_visit[i]],coordy[seq_visit[i]]), s=s)  
      corex<-ld(core)$x
      corey<-ld(core)$y
      Corri1<-rep(2, ncore)
    }
    
    if(seq_visit[i] != seq_visit[i+1]) {
      mig<-simm.bb(date=1:nmig, begin=c(tail(corex,1), tail(corey,1)), end=rnorm(2, c(coordx[seq_visit[i+1]],coordy[seq_visit[i+1]]), sd=25))
      Corri2<-rep(1, nmig)
      sub<-cbind(c(corex, ld(mig)$x), c(corey, ld(mig)$y), c(Corri1, Corri2))
      
    }
    if(seq_visit[i] == seq_visit[i+1]) {
      sub<-cbind(corex, corey, Corri1)
      colnames(sub)<-c("V1", "V2", "V3")
    }
    out<-rbind(out, sub)
  }
  names(out)<-c("x", "y", "Corri")
  out<-as.ltraj(out[,1:2], as.POSIXct(1:nrow(out), origin = "1960-01-01", tz="GMT"), id="id", infolocs=data.frame(out$Corri))
  if(grph==T) {plot(out)}
  return(out)
}



#' Generation of adjacency matrix from movement data
#'
#' Transform an ltraj object to an adjacency matrix using a user-specified grid size
#' @param mov Movement trajectory, need to be a ltraj object
#' @param res Grid size
#' @keywords adj2stack
#' @return A list of object containing the adjacency matrix, the grid use, and patch/corridor identification (only useful if sim_mov was used)
#' @export
#' @examples
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' adj<-traj2adj(traj1, res=100)

traj2adj<-function(mov, res=100) {
  if(require("adehabitatLT") & require("raster")){
    } else {
    print("trying to install packages")
    install.packages(c("adehabitatLT", "raster"))
    if(require("adehabitatLT") & require("raster")){
      print("Packages installed and loaded")
    } else {
      stop("Could not install required packages (raster or adehabitatLT")
    }
  }
  mov<-ld(mov)
  mov[,13]<-1:nrow(mov)
  tt<-SpatialPoints(mov[,1:2]) # HERE ADJUST FOR NA 
  tt1<-apply(coordinates(tt), 2, min)
  tt2<-apply(coordinates(tt), 2, max)
  ras<-raster(xmn=floor(tt1[1]), ymn=floor(tt1[2]),xmx=ceiling(tt2[1]), ymx=ceiling(tt2[2]), res=res)
  values(ras)<-1:ncell(ras)
  patch<-rasterize(tt, ras, field=mov[,13], fun = function(x, ...) round(mean(x)), na.rm=T)
  mov$pix_start<-extract(ras,tt)
  mov$pix_end<-c(mov$pix_start[-1], NA)
  mov$trans<-paste(mov$pix_start, mov$pix_end, sep="_")
  tab<-data.frame(table(mov$trans))
  mov<-merge(mov, tab, by.x="trans", by.y="Var1", all.x=T) # Weights
  mov2<-mov[!duplicated(mov$trans),]
  mat<-matrix(0, nrow=ncell(ras), ncol=ncell(ras))
  for (i in 1:nrow(mov2)) {
    mat[mov2$pix_start[i], mov2$pix_end[i]]<-mov2$Freq[i]
  }
 out<-list(mat, ras, patch)
 class(out)<-"adjmov"
 return(out)
}



#' Calculation of network metrics
#'
#' Transform an adjancency matrix to a series of network metrics at the node-level (weight, degree, betweenness, transitivity, eccenctricity) and graph level (diameter, transitivity, density, and modularity)
#' @param adjmov Adjacency matrix, need to be an object of class adjmov
#' @param grph Whether node level metrics are to be plotted 
#' @keywords traj2adj
#' @return A raster stack object
#' @export
#' @examples
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' stck<-adj2stack(traj2adj(traj1, res=100), grph=T)

adj2stack<-function(adjmov, grph=T) {
  if(require("igraph") & require("raster")){
    } else {
    print("trying to install packages")
    install.packages(c("igraph", "raster"))
    if(require("igraph") & require ("raster")){
      print("Packages installed and loaded")
    } else {
      stop("Could not install required packages (raster or igraph")
    }
  }
  g<-graph_from_adjacency_matrix(adjmov[[1]], mode="directed", weighted = T)
  grid<-stack(adjmov[[3]])
  grid[[2]]<-setValues(grid[[1]], rowSums(adjmov[[1]])/sum(adjmov[[1]]))
  grid[[3]]<-setValues(grid[[1]], degree(g))
  grid[[4]]<-setValues(grid[[1]], betweenness(g))
  grid[[5]]<-setValues(grid[[1]], transitivity(g, type="local")) # Weird so far
  grid[[6]]<-setValues(grid[[1]], eccentricity(g))
  grid[[7]]<-setValues(grid[[1]], diameter(g))
  grid[[8]]<-setValues(grid[[1]], transitivity(g, type="global"))
  grid[[9]]<-setValues(grid[[1]], edge_density(g))
  grid[[10]]<-setValues(grid[[1]], modularity(cluster_walktrap(g)))
  names(grid)<- c("Actual","Weight", "Degree", "Betweenness", "Transitivity", "Eccentricity",  "Diameter", "Global transitivity", "Density", "Modularity")
  
  if(grph==T) plot(grid[[2:6]])
  return(grid)
}



#' Normal mixture model for clustering of node level metrics
#'
#' Apply a normal mixture model to a node-level metric 
#' @param stack An object produce by the function adj2stack
#' @param id Metric to be used (2=Weight, 3=Degree, 4=Betweenness, 5=Transitivity, 6=Eccentricity)
#' @param grph Whether resulting classification should be plotted
#' @keywords adj2stack traj2adj Mclust
#' @return A list object containing a Mclust object and a raster object
#' @export
#' @examples
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' stck<-adj2stack(traj2adj(traj1, res=100), grph=T)
#' cl<-clustnet(stck, id=2, nclust=2, grph=T)
#' summary(cl[[1]])

clustnet<-function(stack, id=2, nclust=2, grph=T) {
  if(require("mclust") & require("raster")){
  } else {
    print("trying to install packages")
    install.packages(c("mclust", "raster"))
    if(require("mclust") & require ("raster")){
      print("Packages installed and loaded")
    } else {
      stop("Could not install required packages (raster or mclust")
    }
  }
  clip<-stack[[1]]*0+1
  clip1<-stack[[id]]*clip
  val<-values(clip1)
  valna<-val[!is.na(values(clip1))]
  clust<-Mclust(valna, nclust) 
  val[!is.na(val)]<-clust$classification
  values(clip)<-val
  if (grph==T) plot(clip)
  return(list(clust, clip))  
}


#' Sample quantile of distance for ltraj object
#'
#' Wrapper function that extract the sample quantile of distance 
#' @param x A ltraj object
#' @param p Probability, default=0.5 (median)
#' @keywords ltraj
#' @return A vector of length p
#' @export
#' @examples
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' stck<-adj2stack(traj2adj(traj1, res=quant(traj1)), grph=T)

quant<-function(x, p=0.5) {quantile(ld(x)$dist, probs=p, na.rm=T)}


#' Extract occupied cells in a raster object 
#'
#' Extract only occupied cells in a raster object, 
#' @param grid An object generated by the function adj2stack
#' @param id Metric to be used (2=Weight, 3=Degree, 4=Betweenness, 5=Transitivity, 6=Eccentricity)
#' @keywords adj2stack
#' @return A vector 
#' @export
#' @examples
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' stck<-adj2stack(traj2adj(traj1, res=quant(traj1)), grph=T)
#' mean(val(stck, 2))

val<-function(grid, id) {
  clip<-grid[[1]]*0+1
  clip1<-grid[[id]]*clip
  clip1<-clip1[!is.na(clip1)]
  return(clip1)
}


#' Summarize graph-level metrics 
#'
#' Summarize graph-level metrics from an object generated by adj2stack 
#' @param grid An object generated by the function adj2stack
#' @param id Metric to be used (2=Weight, 3=Degree, 4=Betweenness, 5=Transitivity, 6=Eccentricity)
#' @keywords adj2stack
#' @return A vector 
#' @export
#' @examples
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' stck<-adj2stack(traj2adj(traj1, res=quant(traj1)), grph=T)
#' graphmet(stck)


graphmet<-function(grid) {
 values(grid[[7:10]])[1,]
}


