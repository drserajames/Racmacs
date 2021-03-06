% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_making.R
\name{make.acmap}
\alias{make.acmap}
\title{Make an antigenic map from scratch}
\usage{
make.acmap(
  number_of_dimensions = 2,
  number_of_optimizations = 100,
  minimum_column_basis = "none",
  fixed_column_bases = NULL,
  move_trapped_points = NULL,
  parallel_optimization = TRUE,
  verbose = TRUE,
  vverbose = FALSE,
  ...
)
}
\arguments{
\item{number_of_dimensions}{The number of dimensions in the map}

\item{number_of_optimizations}{The number of optimization runs to perform}

\item{minimum_column_basis}{The minimum column basis for the map}

\item{move_trapped_points}{How should removal of trapped points be performed (see details)}

\item{...}{}
}
\value{
Returns an antigenic map object of the corresponding class.
}
\description{
Make an antigenic map from scratch
}
\details{
Hunting for and removal of trapped points can be performed for either none of
the optimization runs ("none"), all of the optimization runs ("all") or only
the best one found ("best").
}
\examples{
# Read in an acmap object from a file
map_file_path <- system.file("extdata/h3map2004.ace", package = "Racmacs")
map <- read.acmap(map_file_path)

# Read in a table of titer data and make a new acmap file from it
titer_file_path <- system.file("extdata/h3map2004_hitable.csv", package = "Racmacs")
titer_data <- read.titerTable(titer_file_path)
map <- acmap.cpp(table = titer_data)
}
\seealso{
Other {map optimization functions}: 
\code{\link{moveTrappedPoints}()},
\code{\link{optimizeMap}()},
\code{\link{randomizeCoords}()},
\code{\link{relaxMapOneStep}()},
\code{\link{relaxMap}()}
}
\concept{{map optimization functions}}
