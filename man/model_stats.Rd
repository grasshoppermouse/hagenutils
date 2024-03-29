% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/other_utils.R
\name{model_stats}
\alias{model_stats}
\title{model_stats}
\usage{
model_stats(...)
}
\arguments{
\item{...}{Either a single regression model object, or multiple named models, e.g., m1 = m1, m2 = m2, ...}
}
\value{
A named list of model parameters and statistics
}
\description{
Return a named list of estimates, standard errors, and p-values for one or more regression models.
}
\details{
Uses the tidy functions from broom and broom.mixed to return a named list of model stats that can be easily used inline in rmarkdown documents.
}
\examples{
\dontrun{
if(interactive()){
 m1 <- lm(mpg ~ hp, mtcars)
 m2 <- lm(mpg ~ wt, mtcars)
 x <- model_stats(m1 = m1, m2 = m2)
 
 x$m1$hp$estimate
 x$m1$hp$str # formatted beta coefficient with 95\% CI.
 }
}
}
\seealso{
\code{\link[broom.mixed]{reexports}}
 \code{\link[glue]{glue}}
}
