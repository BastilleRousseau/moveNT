## moveNT  ##

**moveNT** is an R package that analyzes animal movement using Network Theory.

This is the development area for the package `moveNT`, which provides a series of functions to analyze movement data using Network Theory. 

*References*: Bastille-Rousseau, G., Douglas-Hamilton, I., Blake, S., Northrup, J., Wittemyer, G. (2018) Applying network theory to animal movements to identify properties of landscape space use. Ecological Applications. 

For questions: gbr |at| siu.edu

## Installation of the development version  ##

You need to use the package `devtools`from Hadley Wickham. 
    
    library(devtools)
    install_github("BastilleRousseau/moveNT")


## Getting started ##

The package main functions are `traj2adj` and `adj2stack`. The package also contains a function to simulate movement data `sim_mov`. For a list of documented functions see the Reference manual. 
Two vignettes have been created. `moveNT`shows the basic functions and also presents the use of new functions such as `loop`, `interpolation`, and `mosaic_network` that are helpful for mapping. 
The `movescape`vignette introduce a series of new functions that perform machine learning and combine different metrics into a single representation. The manuscript is currently in review. 





