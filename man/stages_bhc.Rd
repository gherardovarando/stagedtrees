% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stages_bhc.R
\name{stages_bhc}
\alias{stages_bhc}
\title{Backward hill-climbing}
\usage{
stages_bhc(
  object,
  score = function(x) {
     return(-BIC(x))
 },
  max_iter = Inf,
  scope = NULL,
  ignore = object$name_unobserved,
  trace = 0
)
}
\arguments{
\item{object}{an object of class \code{sevt} with fitted probabilities and
data, as returned by \code{full} or \code{sevt_fit}.}

\item{score}{the score function to be maximized.}

\item{max_iter}{the maximum number of iterations per variable.}

\item{scope}{names of variables that should be considered for the optimization.}

\item{ignore}{vector of stages which will be ignored and left untouched,
by default the name of the unobserved stages stored in
\code{object$name_unobserved}.}

\item{trace}{if >0 increasingly amount of info
is printed (via \code{cli::cli_text}).}
}
\value{
The final staged event tree obtained.
}
\description{
Greedy search of staged event trees with
iterative joining of stages.
}
\details{
For each variable the algorithm tries to join stages
and moves to the best model that increases the score. When no
increase is possible it moves to the next variable.
}
\examples{
DD <- generate_xor_dataset(p = 4, n = 100)
model <- stages_bhc(full(DD), trace = 2)
summary(model)
}
