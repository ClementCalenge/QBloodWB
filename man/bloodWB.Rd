\name{bloodWB}
\alias{bloodWB}
\docType{data}
\title{
  Main Dataset Used to Explore How Wild Boar Blood Quality Changes After
  Death
}
\description{
  This dataset contains the description of blood biochemistry of 20 wild
  boar after their death.
}
\usage{data("bloodWB")}
\format{
  A data frame with 20 observations on the following 22 variables.
  \describe{
    \item{\code{Hemolysis}}{Degree of hemolysis (grade 0 to 4)}
    \item{\code{Albumin}}{Blood concentration of albumin}
    \item{\code{AlanineAT}}{Blood concentration of alanine aminotransferase}
    \item{\code{AspartateAT}}{Blood concentration of aspartate aminotransferase}
    \item{\code{CreatineKinase}}{Blood concentration of creatine kinase}
    \item{\code{Chloride}}{Blood concentration of chloride}
    \item{\code{Creatinine}}{Blood concentration of creatinine}
    \item{\code{Fructosamine}}{Blood concentration of fructosamine}
    \item{\code{Fe}}{Blood concentration of Iron}
    \item{\code{Globulins}}{Blood concentration of globulins}
    \item{\code{K}}{Blood concentration of potassium}
    \item{\code{Na}}{Blood concentration of sodium}
    \item{\code{Alk.Phosphatase}}{Blood concentration of alkaline phosphatase}
    \item{\code{TotalProtein}}{Blood concentration of total protein}
    \item{\code{UreaNitrogen}}{Blood concentration of urea nitrogen}
    \item{\code{age}}{a factor with levels \code{< 1 yo} \code{1-3 yo} \code{>3 yo}}
    \item{\code{sex}}{a character vector (either "M", or "F")}
    \item{\code{TimeColl}}{numeric vector indicating the hour at which
      the sample was collected (stored as a decimal number, e.g. 11.5
      stands for 11h30 AM).}
    \item{\code{Temperature}}{a numeric vector indicating the body
      temperature (NA indicates that the temperature was below detection
      limits).}
    \item{\code{Aspect}}{a character vector describing the aspect of the
      blood sample}
    \item{\code{Hour}}{a character vector storing the hour (same as
      TimeColl, but rounded to the nearest hour).}
    \item{\code{TimeSinceColl}}{a numeric vector storing the number of
      hours elapsed since death.}
  }
}
\source{
  Larrat, S. et al. (in prep.). Effect of Post-Mortem Interval on the Biochemical Quality of Blood in Wild Boar (Sus scrofa).
}
\references{
  Larrat, S. et al. (in prep.). Effect of Post-Mortem Interval on the Biochemical Quality of Blood in Wild Boar (Sus scrofa).
}
\examples{
data(bloodWB)
head(bloodWB)
}
\keyword{datasets}
