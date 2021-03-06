% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_subset.R
\name{subsetPoints}
\alias{subsetPoints}
\alias{subsetAntigens}
\alias{subsetSera}
\title{Subset antigens and sera}
\usage{
subsetAntigens(map, antigens)

subsetSera(map, sera)
}
\arguments{
\item{map}{The map data object}

\item{antigens}{Antigens to keep (specified by name or index)}

\item{sera}{Sera to keep (specified by name or index)}
}
\description{
Functions to keep only a subset of antigens and sera from a map
}
\seealso{
Other {functions for working with map data}: 
\code{\link{acmap}()},
\code{\link{as.json}()},
\code{\link{cloneMap}()},
\code{\link{convertingMaps}},
\code{\link{orderPoints}},
\code{\link{read.acmap}()},
\code{\link{read.titerTable}()},
\code{\link{removePoints}},
\code{\link{save.acmap}()},
\code{\link{save.coords}()},
\code{\link{save.titerTable}()},
\code{\link{selectedOptimization}()},
\code{\link{subsetMap}()}
}
\concept{{functions for working with map data}}
