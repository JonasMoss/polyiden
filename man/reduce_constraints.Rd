% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symmetry.R
\name{reduce_constraints}
\alias{reduce_constraints}
\title{Reduce a symmetric matrix problem.}
\usage{
reduce_constraints(constraints)
}
\arguments{
\item{constraints}{List of constraints.}
}
\value{
List of reduced constraints.
}
\description{
The inequality constraints in the symmetric matrix problem sometimes has
an empty interior. In this case \code{stats::constrOptim} will not work, and
we have to remove redunant variables.
}
