library(rpanel)

## example 1
m <- pcrfit(reps, 1, 2, l5)
x11(width=4,height=4)

draw  <-  function(panel)  {
thresh <- panel$t
efficiency(m, thresh = thresh)
panel
}

panel  <-  rp.control(t  =  0)
rp.slider(panel, t,  0,  5,  draw, showvalue  =  TRUE)

