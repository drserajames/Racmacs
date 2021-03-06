% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_methods.R
\name{view.rac}
\alias{view.rac}
\title{Viewing racmap objects}
\usage{
\method{view}{rac}(
  map,
  ...,
  .jsCode = NULL,
  .jsData = NULL,
  select_ags = NULL,
  select_sr = NULL,
  show_procrustes = NULL,
  show_stressblobs = NULL,
  keep_all_optimization_runs = FALSE,
  options = list()
)
}
\arguments{
\item{map}{The acmap data object}

\item{...}{Additional arguments to be passed to \code{\link{RacViewer}}}

\item{.jsCode}{Additional javascript code to be run after map has been loaded and rendered}

\item{.jsData}{Any data to supply to the .jsCode function}

\item{select_ags}{A vector of antigen indices to select in the plot}

\item{select_sr}{A vector of serum indices to select in the plot}

\item{show_procrustes}{If the map contains procrustes information, should procrustes lines be shown by default?}

\item{show_stressblobs}{If the map contains stress blob information, should stress blobs be shown by default?}

\item{keep_all_optimization_runs}{Should information on all the optimization runs be kept in the viewer, or just view the currently selected optimisation run.}

\item{options}{A named list of viewer options to pass to \code{RacViewer.options()}}
}
\description{
View a racmap object in the interactive viewer.
}
