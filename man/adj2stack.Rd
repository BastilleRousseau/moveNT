% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcts packages.R
\name{adj2stack}
\alias{adj2stack}
\title{Calculation of network metrics}
\usage{
adj2stack(adjmov, grph = T, mode = "directed", weighted = T, ...)
}
\arguments{
\item{adjmov}{Adjacency matrix, need to be an object produced by function traj2adj}

\item{grph}{Whether node level metrics are to be plotted}

\item{mode}{Whether the graph should be "directed" or "undirected. Default="directed". See "graph_from_adjacency_matrix" from package "igraph"}

\item{weighted}{Whether the graph should be weighted (=TRUE) or unweighted (= NULL). Default is weighted. See "graph_from_adjacency_matrix" from package "igraph"}
}
\value{
A raster stack object
}
\description{
Transform an adjancency matrix to a series of network metrics at the node-level (weight, degree, betweenness, transitivity, eccenctricity) and graph level (diameter, transitivity, density, and modularity)
}
\examples{
traj1<-sim_mov(type="OU", npatches=3, grph=T)
stck<-adj2stack(traj2adj(traj1, res=100), grph=T)
}
\keyword{traj2adj}
