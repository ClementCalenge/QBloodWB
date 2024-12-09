\name{QBloodWB-package}
\alias{QBloodWB-package}
\alias{QBloodWB}
\docType{package}
\title{
\packageTitle{QBloodWB}
}
\description{
\packageDescription{QBloodWB}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{QBloodWB}
\packageIndices{QBloodWB}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
\packageAuthor{QBloodWB}

Maintainer: \packageMaintainer{QBloodWB}
}
\references{
  Larrat, S. et al. (in prep.). Effect of Post-Mortem Interval on the Biochemical Quality of Blood in Wild Boar (Sus scrofa).
}
\keyword{ package }
\examples{

data(bloodWB)

## Get the biochemical parameters
a <- bloodWB[,1:15]

## The time
time <- bloodWB$TimeSinceColl

## Remove WB 19 before PCA (avoid leverage effect;
## (see vignette for further details)
ab <- a[-19,]

## PCA
pc <- ade4::dudi.pca(ab, scan=FALSE)

## Inserts 19 as supplementary row
axis1 <- c(pc$li[-19,1],ade4::suprow(pc, a[19,])$lisup[1,1],pc$li[19,1])

##  Show the SSR
ssr <- rssSRCuts(time, axis1, cutlim=c(0,5))
plot(c(0:5), ssr, ty="l", ylim=range(ssr),
     xlab="Threshold",
     ylab="Residual sum of squares")


## Seems minimal after two hours

## Bootstrap
(boo <- bootSR(time, axis1, cutlim=c(1,4), nBoot=1000))

}
