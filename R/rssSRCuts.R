rssSRCuts <-
function(time, response, cutlim=c(0,5))
{
    ssr <- sapply(cutlim[1]:cutlim[2], function(coupure) {
        te <- time*(time>(coupure+0.5))
        rssSR(lm(response~te))
    })
    return(ssr)
}
