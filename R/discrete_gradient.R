#######################################################################
### Some helper functions to make colour ramps look a little better ###
#######################################################################

# Have a look at help on stack overflow here:
# https://stackoverflow.com/questions/50148905/how-to-create-colorbars-in-ggplot-similar-to-those-created-by-lattice
# https://github.com/tidyverse/ggplot2/issues/2673#issuecomment-402878574

discrete_gradient_pal <- function(colours, bins) {
  ramp <- scales::colour_ramp(colours)
  
  function(x) {
    if (length(x) == 0) return(character())
    
    i <- floor(x * bins)
    i <- ifelse(i > bins-1, bins-1, i)
    ramp(i/(bins-1))
  }
}

scale_colour_discrete_gradient <- function(..., colours, bins, na.value = "grey50", guide = "colourbar", aesthetics = c("colour", "fill"), colors)  {
  colours <- if (missing(colours)) 
    colors
  else colours
  continuous_scale(
    aesthetics,
    "discrete_gradient",
    discrete_gradient_pal(colours, bins),
    na.value = na.value,
    guide = guide,
    ...
  )
}

# ends