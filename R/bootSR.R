bootSR <-
function(time, response, cutlim=c(1,4), nBoot=1000) {

    re <- lapply(1:nBoot, function(k) {

        ind <- sample(c(1:length(time)), length(time), replace = TRUE)
        tempscr <- time[ind]
        responser <- response[ind]
        ssr <- sapply(cutlim[1]:cutlim[2], function(coupure) {
            te <- tempscr*(tempscr>(coupure+0.5))
            rssSR(lm(responser~te))
        })
        if (!any(tempscr>4.5)) {
            ssr[6] <- 1000000000000000
        }
        if (!any(tempscr>3.5)) {
            ssr[5:6] <- 1000000000000000
        }
        return(ssr)
    })
    return(table(sapply(re,which.min))/1000)
}
