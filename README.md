## moveNT  ##

**moveNT** is an R package that analyzes animal movement using Network Theory.

This is the development area for the package `moveNT`, which provides a series of functions to analyze movement data using Network Theory. 

*References*: Bastille-Rousseau, G., Douglas-Hamilton, I., Blake, S., Northrup, J., Wittemyer, G. (2018) Applying network theory to animal movements to identify properties of landscape space use. Ecological Applications. 

For questions: gbr |at| colostate.edu

## Installation of the development version  ##

You need to use the package `devtools`from Hadley Wickham. 
    
    library(devtools)
    install_github("BastilleRousseau/moveNT")


## Getting started ##

The package main functions are `traj2adj` and `adj2stack`. The package also contains a function to simulate movement data `sim_mov`. For a list of documented functions see the Reference manual. 
For examples of how to use the main functions, you can also look at the vignette. 

Alternatively, here is a quick example to get you going: 

    traj1<-sim_mov(type="OU", npatches=3, grph=T)
    stck<-adj2stack(traj2adj(traj1, res=100), grph=T)



