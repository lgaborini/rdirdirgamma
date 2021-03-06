% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{get_summary_statistics_cpp}
\alias{get_summary_statistics_cpp}
\title{Compute summary statistics.}
\usage{
get_summary_statistics_cpp(mtx, use_optimized_summary)
}
\arguments{
\item{mtx}{a data matrix (nxp)}

\item{use_optimized_summary}{if TRUE, return the optimized summary statistics (mean, sd, kurtosis, skewness), else standard (mean, sd)}
}
\value{
a kxp matrix of summary statistics
}
\description{
Compute summary statistics.
}
\seealso{
Other ABC summary functions: 
\code{\link{get_number_summary_statistics}()},
\code{\link{get_optimized_summary_statistics_cpp}()},
\code{\link{get_standard_summary_statistics_cpp}()}
}
\concept{ABC summary functions}
