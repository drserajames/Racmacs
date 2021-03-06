% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_props_chart.R
\name{name}
\alias{name}
\alias{name<-}
\title{Getting and setting the map name}
\usage{
name(map, .name = TRUE)
name(map, .check = TRUE) <- value
}
\arguments{
\item{map}{The acmap data object}

\item{.name}{Should the output be labelled with corresponding antigen / sera names}

\item{.check}{Should normal checks be applied for validity be applied when setting this value (used internally)}
}
\value{
Returns either the requested attribute when using a getter function or the updated acmap object when using the setter function.
}
\description{
You can use the standard \code{name()} function to get and set the map name.
}
\seealso{
Other {map attribute functions}: 
\code{\link{acmapAttributes}},
\code{\link{titerTableLayers}()},
\code{\link{titerTable}()}
}
\concept{{map attribute functions}}
