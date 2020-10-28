# # #### Package creation
# # install.packages("devtools")
# library("devtools")
# devtools::install_github("klutometis/roxygen")
# library(roxygen2)
# create("moveNT")
# setwd("./moveNT")
# document()
# #
# # #Create pdf
# pack <- "moveNT"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),"CMD", "Rd2pdf", shQuote(path)))
# #
# #


#' Simulation of patch-based movement trajectory
#'
#' Simulate a movement trajectory with a user defined number of patches and interpatch movement
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

  coordx<-sample(seq(0,20,2), npatches, replace=F)*spacecore
  coordy<-sample(seq(0,20,2), npatches, replace=F)*spacecore
  nmig=ncore/ratio
  out<-data.frame()
  for (i in 1:(nswitch-1)){

    if(type=="2states") {
      core<-moveHMM::simData(nbAnimals=1,nbStates=2,stepDist=stepDist,angleDist=angleDist,stepPar=stepPar, anglePar=anglePar,zeroInflation=F,obsPerAnimal=ncore)
      corex<-core$x+coordx[seq_visit[i]]
      corey<-core$y+coordy[seq_visit[i]]
      Corri1<-rep(2, ncore)
    }

    if(type=="OU") {
      core<-adehabitatLT::simm.mou(date=1:ncore, b=c(coordx[seq_visit[i]],coordy[seq_visit[i]]), s=s)
      corex<-ld(core)$x
      corey<-ld(core)$y
      Corri1<-rep(2, ncore)
    }

    if(seq_visit[i] != seq_visit[i+1]) {
      mig<-adehabitatLT::simm.bb(date=1:nmig, begin=c(tail(corex,1), tail(corey,1)), end=rnorm(2, c(coordx[seq_visit[i+1]],coordy[seq_visit[i+1]]), sd=25))
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
  out<-adehabitatLT::as.ltraj(out[,1:2], as.POSIXct(1:nrow(out), origin = "1960-01-01", tz="GMT"), id="id", infolocs=data.frame(out$Corri))
  if(grph==T) {plot(out)}
  return(out)
}


#' Generation of adjacency matrix from movement data
#'
#' Transform an ltraj object to an adjacency matrix using a user-specified grid size
#' @param mov Movement trajectory, need to be a ltraj object
#' @param res Grid size (based on coordinate system of movement trajectory)
#' @param grid User specified grid (a raster), needs to have a larger extent than the movement trajectory
#' @keywords adj2stack
#' @return A list of objects containing the adjacency matrix, the grid use, and patch/corridor identification (only useful if sim_mov was used)
#' @export
#' @examples
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' adj<-traj2adj(traj1, res=100)



traj2adj<-function(mov, res=100, grid=NULL) {

  mov<-adehabitatLT::ld(mov)
  mov[,13]<-1:nrow(mov)
  tt<-sp::SpatialPoints(mov[,1:2])
  tt1<-apply(coordinates(tt), 2, min)
  tt2<-apply(coordinates(tt), 2, max)
  if(is.null(grid)){ras<-raster(xmn=floor(tt1[1]), ymn=floor(tt1[2]),xmx=ceiling(tt2[1]), ymx=ceiling(tt2[2]), res=res)}
  if(!is.null(grid)){ras<-crop(grid, tt)}
  values(ras)<-1:ncell(ras)
  patch<-rasterize(tt, ras, field=mov[,13], fun = function(x, ...) round(mean(x)), na.rm=T)
   mov$pix_start2<-extract(ras,tt)
 mov$pix_start<-as.numeric(as.factor(mov$pix_start2))
 tt<-values(patch)
 tt[!is.na(tt)]<-1:max(mov$pix_start, na.rm=T)
 values(patch)<-tt
  mov$pix_end<-c(mov$pix_start[-1], NA)
  mov$trans<-paste(mov$pix_start, mov$pix_end, sep="_")
  tab<-data.frame(table(mov$trans))
  mov<-merge(mov, tab, by.x="trans", by.y="Var1", all.x=T) # Weights
  mov2<-mov[!duplicated(mov$trans),]
  mat<-matrix(0, nrow=max(mov2$pix_start,na.rm=T), ncol=max(mov2$pix_end, na.rm=T))
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
#' @param adjmov Adjacency matrix, need to be an object produced by function traj2adj
#' @param grph Whether node level metrics are to be plotted
#' @param mode Whether the graph should be "directed" or "undirected. Default="directed". See "graph_from_adjacency_matrix" from package "igraph"
#' @param weighted Whether the graph should be weighted (=TRUE) or unweighted (= NULL). Default is weighted. See "graph_from_adjacency_matrix" from package "igraph"
#' @keywords traj2adj
#' @return A raster stack object
#' @export
#' @examples
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' stck<-adj2stack(traj2adj(traj1, res=100), grph=T)

adj2stack<-function(adjmov, grph=T, mode="directed", weighted=T, ...) {

  g<-igraph::graph_from_adjacency_matrix(adjmov[[1]], mode=mode, weighted = weighted)
  grid<-stack(adjmov[[3]])
  tt<-values(grid)
  tt[!is.na(tt)]<-rowSums(adjmov[[1]])/sum(adjmov[[1]])
  grid[[2]]<-setValues(grid[[1]], tt)
  tt<-values(grid[[1]])
  tt[!is.na(tt)]<-diag(adjmov[[1]])/sum(adjmov[[1]])
  grid[[3]]<-setValues(grid[[1]], tt)
    tt<-values(grid[[1]])
  tt[!is.na(tt)]<-igraph::degree(g)
  grid[[4]]<-setValues(grid[[1]], tt)
  tt<-values(grid[[1]])
  tt[!is.na(tt)]<-igraph::betweenness(g)
  grid[[5]]<-setValues(grid[[1]], tt)
  tt<-values(grid[[1]])
  tt[!is.na(tt)]<-igraph::transitivity(g, type="local")
  grid[[6]]<-setValues(grid[[1]], tt)
  tt<-values(grid[[1]])
  tt[!is.na(tt)]<-igraph::eccentricity(g)
  grid[[7]]<-setValues(grid[[1]], tt)
  grid[[8]]<-setValues(grid[[1]], igraph::diameter(g))
  grid[[9]]<-setValues(grid[[1]], igraph::transitivity(g, type="global"))
  grid[[10]]<-setValues(grid[[1]], igraph::edge_density(g))
  grid[[11]]<-setValues(grid[[1]], igraph::modularity(igraph::cluster_walktrap(g)))
  names(grid)<- c("Actual","Weight", "Self-loop", "Degree",  "Betweenness", "Transitivity", "Eccentricity",  "Diameter", "Global transitivity", "Density", "Modularity")
  if(grph==T) plot(grid[[2:7]])
  return(grid)
}

#' Looping over all individuals
#'
#' Extract the adjancency matrix and calculate network metrics for all individuals in a trajectory object. Also calculate mean speed, mean direction, and dot product of turning angles
#' @param traj An object produce by the function adehabitatLT with multiple individuals
#' @param res Grid size, will be apply to all individuals
#' @keywords adj2stack traj2adj
#' @return A list object containing a raster stack object for each individual
#' @export
#' @examples
#' data(puechabonsp)
#' locs <- puechabonsp$relocs
#' xy <- coordinates(locs)
#' df <- as.data.frame(locs)
#' da <- as.character(df$Date)
#' da <- as.POSIXct(strptime(as.character(df$Date),"%y%m%d", tz="Europe/Paris"))
#' litr <- as.ltraj(xy, da, id = df$Name)
#' out1<-loop(litr)
loop<-function(traj, res=100 ){
  tt<-SpatialPoints(ld(traj)[,1:2])
  tt1<-apply(coordinates(tt), 2, min)
  tt2<-apply(coordinates(tt), 2, max)
  ras<-raster(xmn=floor(tt1[1]), ymn=floor(tt1[2]),xmx=ceiling(tt2[1]), ymx=ceiling(tt2[2]), res=res)
  id<-unique(adehabitatLT::id(traj))
  id2<-adehabitatLT::id(traj)
  out<-list()
  for (i in 1:length(id)) {
    try(out[[i]]<-adj2stack(traj2adj(traj[which(id2==id[i])], res=res, grid=ras), grph=F))
    pt<-ld(traj[which(id2==id[i])])
    points<-SpatialPoints(pt[,1:2])
    try(out[[i]][[11]]<-rasterize(points, out[[i]][[1]], pt$dist, fun=mean))
    try(out[[i]][[12]]<-rasterize(points, out[[i]][[1]], pt$abs.angle, fun=mean))
    try(out[[i]][[13]]<-rasterize(points, out[[i]][[1]], pt$rel.angle, fun=dot))
    cat(id[i], '\n')
    names(out[[i]])[11:13]<-c("Speed", "Abs angle","DotP")
  }

  return(out)
  }

#' Interpolation based on movement steps for all individuals
#'
#' Use movement steps to linearly interpolate raster produced by loop. User can select if the mean or max is taken when multiple steps overlap in a single pixel. Function need to be applied following the loop function.  This process is very slow.
#' @param traj An object produce by the function adj2stack
#' @param ls An object produced by the loop
#' @param wei Whether mean or max should be used for weight (default = mean)
#' @param deg Whether mean or max should be used for degree (default = mean)
#' @param bet Whether mean or max should be used for betweeness (default = max)
#' @param spe Whether mean or max should be used for speed (default = mean)
#' @param dt Whether mean, max, or dot product should be used for turning angle (default = dot)
#' @keywords adj2stack traj2adj loop
#' @return A list of object containing a raster stack object for each individual
#' @export
#' @examples
#' data(puechabonsp)
#' locs <- puechabonsp$relocs
#' xy <- coordinates(locs)
#' df <- as.data.frame(locs)
#' da <- as.character(df$Date)
#' da <- as.POSIXct(strptime(as.character(df$Date),"%y%m%d", tz="Europe/Paris"))
#' litr <- as.ltraj(xy, da, id = id)
#' out1<-loop(litr)
#' out2<-interpolation(litr, out1)

interpolation<-function(traj, ls, wei=mean, deg=mean, bet=max, spe=mean, dt=dot) {
id<-unique(adehabitatLT::id(traj))
id2<-id(adehabitatLT::traj)
out_fill<-list()
pb <- txtProgressBar(min = 1, max = length(id), style = 3)
for (i in 1:length(id)) {
  setTxtProgressBar(pb, i)
  try({data<-ld(traj[which(id2==id[i])])
  data2<-na.omit(cbind(data[,1],data[,2], data[,1]+data[,4],data[,2]+data[,5], data[,6:7], data[,9]))
  steps<-SpatialLines(apply(data2, 1, function(r) {
    Lines(list(sp::Line(cbind(r[c(1,3)], r[c(2,4)]))), uuid::UUIDgenerate())
  }))
  pts<-extract(ls[[i]], data2[,1:2])
  weight<-rasterize(steps, ls[[i]][[1]], field=pts[,2], fun=wei)
  degree<-rasterize(steps, ls[[i]][[1]], field=pts[,4], fun=deg)
  between<-rasterize(steps, ls[[i]][[1]], field=pts[,5], fun=bet)
  speed<-rasterize(steps, ls[[i]][[1]], field=data2$dist/data2$dt, fun=spe)
  dotp<-rasterize(steps, ls[[i]][[1]], field=data2[,7], fun=dt)
  out_fill[[i]]<-stack(weight, degree, between, speed, dotp)
  })}
return(out_fill)}


#' Mosaic (combine) individual raster together for a given variable
#'
#' Use output of loop or interpolation and combine all individuals (mosaic) together using the mean or max values.
#' @param ls An object produced by the loop or interpolate functions
#' @param index Index indicating which layer to take in the stack
#' @param sc Whether to scale all individual rasters (default = TRUE)
#' @param fun Whether mean or max should be used as the mosaic function (default = mean)
#' @keywords adj2stack traj2adj loop interpolation
#' @return A raster layer object.
#' @export
#' @examples
#' data(puechabonsp)
#' locs <- puechabonsp$relocs
#' xy <- coordinates(locs)
#' df <- as.data.frame(locs)
#' da <- as.character(df$Date)
#' da <- as.POSIXct(strptime(as.character(df$Date),"%y%m%d", tz="Europe/Paris"))
#' litr <- as.ltraj(xy, da, id = id)
#' out1<-loop(litr)
#' mean_weight<-mosaic_network(out1, index=2, sc=T, fun=mean) #Perform mean weight (not-interpolated)
#' plot(mean_weight)
mosaic_network<-function(ls, index=2, sc=T, fun=mean){
  layers<-lapply(ls, function(x) x[[index]])
  if(sc==T) {layers<-lapply(layers, raster::scale)}
  names(layers)[1:2]<-c("x", "y")
  layers$fun<-fun
  layers$na.rm<-TRUE
  layers_mosaic<-do.call(mosaic, layers)
  return(layers_mosaic)
}


#' dot product
dot<-function(x, ...) {
  sum(abs(cos(x)), na.rm=T)/(sum(!is.na(x)))
}




#' Normal mixture model for clustering of single node level metric
#'
#' Apply a normal mixture model to a single node-level metric
#' @param stack An object produce by the function adj2stack (not compatible with loop or interpolation)
#' @param id Metric to be used (2=Weight, 3=Degree, 4=Betweenness, 5=Transitivity, 6=Eccentricity)
#' @param grph Whether resulting classification should be plotted
#' @keywords adj2stack traj2adj Mclust
#' @return A list of object containing a Mclust object and a raster object
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
  clust<-mclust::Mclust(valna, nclust)
  val[!is.na(val)]<-clust$classification
  values(clip)<-val
  if (grph==T) plot(clip)
  return(list(clust, clip))
}


#' Sample quantile of distance for ltraj object
#'
#' Wrapper function that extract the sample quantile of distance of a trajectory object
#' @param x A ltraj object
#' @param p Probability, default=0.5 (median)
#' @keywords ltraj
#' @return A vector of length p
#' @export
#' @examples
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' stck<-adj2stack(traj2adj(traj1, res=quant(traj1)), grph=T)

quant<-function(x, p=0.5) {quantile(adehabitatLT::ld(x)$dist, probs=p, na.rm=T)}


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
#' @keywords adj2stack
#' @return A vector
#' @export
#' @examples
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' stck<-adj2stack(traj2adj(traj1, res=quant(traj1)), grph=T)
#' graphmet(stck)

graphmet<-function(grid) {
 values(grid[[8:11]])[1,]
}




######################
### Movescape functions
######################


#' Convert a list of adj2stack object to a data.frame for clustering
#'
#' Convert output of loop function to a data.frame.
#' @param traj The trajectory used in loop (a traj object)
#' @param grid The output of the loop function
#' @keywords adj2stack traj2adj loop
#' @return A data.frame object.
#' @export
#' @examples
#' data(albatross)
#' grid<-loop(albatross, 35000)
#' table_grid<-table_cluster(albatross, grid)
#' head(table_grid)
table_cluster<-function(traj, grid) {
  id<-unique(adehabitatLT::id(traj))
  if(length(id)!=length(grid)) {stop("traj and grid don't have the same number of individuals")}
  out<-data.frame()
  for (i in 1:length(id)) {
    tt1<-data.frame(values(grid[[i]]))
    tt3<-data.frame(na.omit(tt1))
    tt3$ID<-id[i]
    out<-rbind(out, tt3)
  }
  names(out)[11:13]<-c("Speed", "Abs angle","DotP")
  return(out)
}


#' Individual-level clustering of movement metrics
#'
#' Perform individual-level clustering (first step) of movement metrics. This function uses the output of table_cluster and perform a mixture-model. Users can select which variables will be used and the maximum number of clusters. See also mclust
#' @param table An output from the table_cluster function
#' @param max.n.clust The maximum number of clusters to test, see the documentation for mclust for more information. Default = 8.
#' @param modelname The model structure of the clustering, see the documentation for mclust for more information. Default is equal mean and variance for each clusters (EEV).
#' @param vars The variable to be included. Default = c("Weight", "Degree", "Betweenness", "Speed", "DotP")
#' @keywords adj2stack traj2adj loop table_cluster pop_clust
#' @return A list object with each element representing an individual.
#' @export
#' @examples
#' data(albatross)
#' grid<-loop(albatross, 35000)
#' table_grid<-table_cluster(albatross, grid)
#' ls_ind<-ind_clust(table_grid, max.n.clust=8)
#' table(unlist(lapply(ls_ind, function(x) x$G)))
ind_clust<-function(table, max.n.clust=8, modelname="EEV", vars=c("Weight", "Degree", "Betweenness", "Speed", "DotP")) {
  if(require("mclust")){
    print("mclust is loaded correctly")
  } else {
    print("trying to install mclust")
    install.packages("mclust")
    if(require(mclust)){
      print("mclust installed and loaded")
    } else {
      stop("could not install mclust")
    }
  }
  id<-unique(table$ID)
  ls<-list()
  for (i in 1:length(id)) {
    tt<-scale(table[table$ID==id[i],vars])
    try(ls[[i]]<-Mclust(tt, G=2:max.n.clust, modelNames=modelname))
    print(id[i])
  }
  return(ls)
}

#' Population-level clustering of movement metrics
#'
#' Combine individual-level clustering of movement metrics into a population-level clustering (second step). Users can  the maximum number of clusters. See also mclust
#' @param traj The trajectory object
#' @param ls_ind Individual-level clustering object, the output of ind_clust.
#' @param max.n.clust The maximum number of clusters to test, see the documentation for mclust for more information. Default = 8.
#' @keywords adj2stack traj2adj loop table_cluster ind_clust
#' @return A list object with each element representing an individual.
#' @export
#' @examples
#' data(albatross)
#' grid<-loop(albatross, 35000)
#' table_grid<-table_cluster(albatross, grid)
#' ls_ind<-ind_clust(table_grid, max.n.clust=8)
#' pop<-pop_clust(albatross, ls_ind)
#' pop[[1]]$parameters$mean
#' pop[[1]]$parameters$pro
pop_clust<-function(traj, ls, max.n.clust=8) {
  id<-unique(adehabitatLT::id(traj))
  nvar<-nrow(ls[[1]]$parameters$mean)
  coef<-data.frame()
  for (i in 1:length(id)) {
    G<-ls[[i]]$G
    gg<-data.frame(cbind(t(ls[[i]]$parameters$mean), ls[[i]]$parameters$pro, rep(id[i], G))   )
    coef<-rbind(coef, gg)
  }
  coef[,-ncol(coef)] <- sapply(coef[-ncol(coef)],function(x) as.numeric(as.character(x)))
  names(coef)[nvar+1]<-"Prop"
  names(coef)[nvar+2]<-"ID"
  clust<-Mclust(coef[,1:nvar], G=1:max.n.clust)
  coef$clust<-clust$classification
  out<-list(clust, coef)
  return(out)
}



#' Back-association of population-level clustering to individual clusters
#'
#' Generate individual-level rasters of population-level clustering. For each individual, the fuction generate a raster stack containing a raster of the most likely cluster, and several rasters giving the probability of observing each cluster.
#' @param grid The output of the loop function
#' @param pop_clust The output of the pop_clust function
#' @param ind_clust The output of the ind_clust function
#' @param table The output of table_cluster
#' @keywords adj2stack traj2adj loop table_cluster ind_clust pop_clust
#' @return A list of raster stack object.
#' @export
#' @examples
#' data(albatross)
#' grid<-loop(albatross, 35000)
#' table_grid<-table_cluster(albatross, grid)
#' ls_ind<-ind_clust(table_grid, max.n.clust=8)
#' pop<-pop_clust(albatross, ls_ind)
#' clust_stack<-clust_stack(grid, pop, ls_ind, table_grid)
#' plot(clust_stack[[1]])
clust_stack<-function(grid, pop_clust, ind_clust, table) {
  id<-unique(table$ID)
  coef2<-cbind(pop_clust[[2]], pop_clust[[1]]$z)
  out_ls<-list()
  n.clust1<-length(unique(coef2$clust))

  for (i in 1:length(grid)) {
    coef3<-coef2[coef2$ID==id[i],]
    n.clust<-length(coef3$clust)
    cl<-((coef3$clust))
    class<-data.frame(cbind(1:length(cl), cl))
    prop<-coef3[,as.character(cl)]
    out2<-table[table$ID==id[i],]
    out2$clust<-ind_clust[[i]]$classification
    out2$id<-1:nrow(out2)
    out2<-merge(out2, class, by.x="clust", by.y="V1", sort=F)
    out2<-out2[order(out2$id),]
    out2<-cbind(out2, ind_clust[[i]]$z)

    for (j in 1:nrow(out2)) {
      out2[j,(ncol(out2)-n.clust+1):ncol(out2)]<- out2[j,(ncol(out2)-n.clust+1):ncol(out2)]*prop[out2$clust[j],] #Multiply the two probability
    }

    dd<-values(grid[[i]])
    dd2<-apply(dd, 1, function(x) max(is.na(x)))
    ind<-which(dd2==0)
    r0<-grid[[i]][[1]]
    values(r0)<-0
    tt <- values(r0)
    gr<-stack(r0)
    tt[ind]<-out2$cl
    gr[[1]]<-setValues(gr[[1]],tt)

    for (z in 2:(n.clust1+1)) { gr[[z]]<-r0 }
    for (k in 1:n.clust) {
      tt[ind]<-out2[,(ncol(out2)-n.clust)+k]
      gg<-setValues(gr[[1]],tt)
      gr[[cl[k]+1]]<-mosaic(gr[[cl[k]+1]], gg, fun=max)
    }
    names(gr)<-c("Clust", paste("Prop", unique(class$cl), sep=""))
    out_ls[[i]]<-gr
    print(id[i])
  }
  return(out_ls)
}

#' Population-level single-use maps of each cluster
#'
#' Produce maps (raster) indicating if at least one individual is using a pixel as a given cluster
#' @param clust_stack The output of clust_stack
#' @keywords adj2stack traj2adj loop table_cluster ind_clust
#' @return A stack object with each raster showing use of each cluster
#' @export
#' @examples
#' data(albatross)
#' grid<-loop(albatross, 35000)
#' table_grid<-table_cluster(albatross, grid)
#' ls_ind<-ind_clust(table_grid, max.n.clust=8)
#' pop<-pop_clust(albatross, ls_ind)
#' clust_stack<-clust_stack(grid, pop, ls_ind, table_grid)
#' pop_stack<-pop_stack(clust_stack)
#' plot(pop_stack)
pop_stack<-function(clust_stack) {
  n.clust<-length(names(clust_stack[[1]]))-1
  ls<-lapply(clust_stack, function(x) x[[1]])
  fct<-function(rast, val) {
    values(rast) <- ifelse(values(rast) %in% c(val), 1, 0) # just replace
    return(rast)
  }
  ls1<-list()
  for (i in 1:n.clust) {ls1[[i]]<-lapply(ls, function(x) fct(x, i)) }

  #Mosaic all together
  x <- ls1[[1]]
  names(x)[1:2] <- c('x', 'y')
  x$fun <- max
  x$na.rm <- TRUE
  y <- do.call(mosaic, x)
  out<-stack(y)
  for (i in 2:n.clust) {
    x <- ls1[[i]]
    names(x)[1:2] <- c('x', 'y')
    x$fun <- max
    x$na.rm <- TRUE
    y <- do.call(mosaic, x)
    out[[i]]<-y
  }
  names(out)<-paste("Clust", c(1:n.clust), sep="_")
  return(out)
}

#' Population-level multi-use map
#'
#' Produce a single map (raster) indicating all types of use found in a cluster.
#' @param clust_stack The output of clust_stack
#' @keywords adj2stack traj2adj loop table_cluster ind_clust
#' @return A stack object with each raster showing use of each cluster
#' @export
#' @examples
#' data(albatross)
#' grid<-loop(albatross, 35000)
#' table_grid<-table_cluster(albatross, grid)
#' ls_ind<-ind_clust(table_grid, max.n.clust=8)
#' pop<-pop_clust(albatross, ls_ind)
#' clust_stack<-clust_stack(grid, pop, ls_ind, table_grid)
#' pop_overl<-pop_overl(clust_stack)
#' table(values(pop_overl))
#' plot(pop_overl)
pop_overl<-function(clust_stack) {
  ls<-lapply(clust_stack, function(x) x[[1]])
  ff<-function(x, na.rm=T) {
    if (na.rm) {x<-na.omit(x)}
    x1<-x
    x1<-x1[!duplicated(x1)]
    x1<-sort(x1)
    as.numeric(paste(x1, collapse=""))
  }
  x <- ls
  names(x)[1:2] <- c('x', 'y')
  x$fun <- ff
  x$na.rm <- TRUE
  y <- do.call(mosaic, x)
}


