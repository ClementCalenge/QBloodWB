\name{bootSR}
\alias{bootSR}
\alias{rssSR}
\alias{rssSRCuts}
\title{
  Residual Sum of Squares of the Segmented Regression Model Used to
  Describe Blood Biochemistry of Wild Boar After Death.
}
\description{
  The function \code{rssSR} calculates the residual sum of squares of a
  linear model (class \code{"lm"}).  The function \code{rssSRCuts}
  calculates the residual sum of square of a segmented regression model
  for every possible threshold value between two limits.  The function \code{bootSR}
  implements the bootstrap approach described in Larrat et
  al. (calculates the proportion of \code{nBoot} bootstrap samples where
  the threshold value in the segmented regression is equal to 1 hour, 2
  hours, etc.). 
}
\usage{
rssSR(mod)

bootSR(time, response, cutlim = c(1, 4), nBoot = 1000)

rssSRCuts(time, response, cutlim = c(0, 5))

}
\arguments{
  \item{mod}{
    An object of class \code{"lm"}
  }
  \item{time}{
    numeric vector containing the number of hours elapsed since death.
  }
  \item{response}{
    numeric vector containing the response variable to be modelled by
    segmented regression.
  }
  \item{cutlim}{
    The min and max limits of the threshold in the segmented regression.
  }
  \item{nBoot}{
    The number of bootstrap samples
  }
}
\value{
  The function \code{rssSR} returns a numeric value (the residual sum of squares of
  the model).  The function \code{rssSRCuts} returns a vector of values
  (the residual sum of squares of a segmented regression model for all
  possible threshold values between \code{cutlim[1]} and
  \code{cutlim[2]}).  The function \code{bootSR} returns a table indicating the
  proportion of the bootstrap samples for which the residual sum of
  squares is minimal for every possible threshold value from
  \code{cutlim[1]} to \code{cutlim[2]}.
}
\references{
  Larrat, S. et al. (in prep.). Effect of Post-Mortem Interval on the Biochemical Quality of Blood in Wild Boar (Sus scrofa).
}
\author{
  Clement Calenge <clement.calenge@ofb.gouv.fr>
}
\seealso{
  \code{\link{bloodWB}}
}
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
\keyword{models}
