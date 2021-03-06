% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_stress.R
\name{mapDistances}
\alias{mapDistances}
\title{Return calculated map distances for an acmap}
\usage{
mapDistances(map, optimization_number = NULL)
}
\arguments{
\item{map}{The acmap data object}

\item{optimization_number}{The optimization data to access (defaults to the currently selected optimization)}
}
\value{
Returns a matrix of map distances with antigens as rows and sera as columns.
}
\description{
Takes the acmap object and calculate distances between antigens and sera for the
currently selected or specified optimization.
}
\seealso{
Other {map diagnostic functions}{functions relating to map stress calculation}: 
\code{\link{logtiterTable}()},
\code{\link{mapResiduals}()},
\code{\link{pointStress}},
\code{\link{recalculateStress}()},
\code{\link{stressTable}()},
\code{\link{tableDistances}()}
}
\concept{{map diagnostic functions}{functions relating to map stress calculation}}
