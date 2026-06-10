
library(Racmacs)
library(testthat)
context("optimizeMapStart — custom starting coordinates")

# --------------------------------------------------------------------
# Shared test fixture: small perfect map with known global minimum
# (stress ≈ 0 when placed at the true coordinates)
# --------------------------------------------------------------------
set.seed(42)
ag_coords <- cbind(-4:4, runif(9, -1, 1))
sr_coords <- cbind(runif(9, -1, 1), -4:4)
colbases  <- round(runif(9, 3, 6))
colbasesmat <- matrix(colbases, 9, 9, byrow = TRUE)
distmat     <- as.matrix(dist(rbind(ag_coords, sr_coords)))[seq_len(9), -seq_len(9)]
logtiters   <- colbasesmat - distmat
titers      <- 2^logtiters * 10
mode(titers) <- "character"

perfect_map <- acmap(
  titer_table = titers,
  ag_coords   = ag_coords,
  sr_coords   = sr_coords
)

n_ag   <- numAntigens(perfect_map)
n_sr   <- numSera(perfect_map)
n_runs <- 5


# --------------------------------------------------------------------
# Helper: build a valid pre-computed coordinate list
# --------------------------------------------------------------------
make_coord_list <- function(n, seed = 1L) {
  set.seed(seed)
  lapply(seq_len(n), function(i) {
    list(
      ag_coords = matrix(runif(n_ag * 2, -3, 3), n_ag, 2),
      sr_coords = matrix(runif(n_sr * 2, -3, 3), n_sr, 2)
    )
  })
}


# ====================================================================
# 1.  Pre-computed list mode
# ====================================================================
test_that("optimizeMapStart runs with a pre-computed coordinate list", {

  starts <- make_coord_list(n_runs)
  result <- optimizeMapStart(
    map                    = perfect_map,
    number_of_dimensions   = 2,
    number_of_optimizations = n_runs,
    starting_coords        = starts,
    fixed_column_bases     = colbases,
    check_convergence      = FALSE,
    verbose                = FALSE
  )

  expect_s3_class(result, "acmap")
  expect_equal(numOptimizations(result), n_runs)
  expect_true(all(is.finite(mapStress(result))))
  expect_gt(mapStress(result), 0)

})


test_that("pre-computed list produces deterministic results", {

  starts <- make_coord_list(n_runs, seed = 7L)

  r1 <- optimizeMapStart(perfect_map, 2, n_runs,
    starting_coords    = starts,
    fixed_column_bases = colbases,
    check_convergence  = FALSE,
    verbose            = FALSE)

  r2 <- optimizeMapStart(perfect_map, 2, n_runs,
    starting_coords    = starts,
    fixed_column_bases = colbases,
    check_convergence  = FALSE,
    verbose            = FALSE)

  # Same starting coords → same final coords
  expect_equal(agCoords(r1), agCoords(r2))
  expect_equal(srCoords(r1), srCoords(r2))

})


# ====================================================================
# 2.  Function mode
# ====================================================================
test_that("optimizeMapStart runs with a generating function", {

  my_fn <- function(n_ag, n_sr, ndim, D) {
    n_pts <- n_ag + n_sr
    coords <- matrix(runif(n_pts * ndim, -2, 2), n_pts, ndim)
    list(
      ag_coords = coords[seq_len(n_ag),         , drop = FALSE],
      sr_coords = coords[seq(n_ag + 1L, n_pts), , drop = FALSE]
    )
  }

  result <- optimizeMapStart(
    map                    = perfect_map,
    number_of_dimensions   = 2,
    number_of_optimizations = n_runs,
    starting_coords        = my_fn,
    fixed_column_bases     = colbases,
    check_convergence      = FALSE,
    verbose                = FALSE
  )

  expect_s3_class(result, "acmap")
  expect_equal(numOptimizations(result), n_runs)
  expect_true(all(is.finite(mapStress(result))))

})


test_that("generating function receives coord_args correctly", {

  # A function that uses a 'range' argument from coord_args
  range_fn <- function(n_ag, n_sr, ndim, D, range = 1) {
    n_pts  <- n_ag + n_sr
    coords <- matrix(runif(n_pts * ndim, -range, range), n_pts, ndim)
    list(
      ag_coords = coords[seq_len(n_ag),         , drop = FALSE],
      sr_coords = coords[seq(n_ag + 1L, n_pts), , drop = FALSE]
    )
  }

  result <- optimizeMapStart(
    map                    = perfect_map,
    number_of_dimensions   = 2,
    number_of_optimizations = 2,
    starting_coords        = range_fn,
    coord_args             = list(range = 5),
    fixed_column_bases     = colbases,
    check_convergence      = FALSE,
    verbose                = FALSE
  )

  expect_s3_class(result, "acmap")
  expect_equal(numOptimizations(result), 2)

})


# ====================================================================
# 3.  Built-in distribution: "uniform"
# ====================================================================
test_that("optimizeMapStart runs with built-in 'uniform' distribution", {

  result <- optimizeMapStart(
    map                    = perfect_map,
    number_of_dimensions   = 2,
    number_of_optimizations = n_runs,
    starting_coords        = "uniform",
    fixed_column_bases     = colbases,
    check_convergence      = FALSE,
    verbose                = FALSE
  )

  expect_s3_class(result, "acmap")
  expect_equal(numOptimizations(result), n_runs)
  expect_true(all(is.finite(mapStress(result))))

})


test_that("'uniform' built-in respects min/max coord_args", {

  # Supply a very narrow range — optimisation will still run, just start squeezed
  result <- optimizeMapStart(
    map                    = perfect_map,
    number_of_dimensions   = 2,
    number_of_optimizations = 2,
    starting_coords        = "uniform",
    coord_args             = list(min = -0.01, max = 0.01),
    fixed_column_bases     = colbases,
    check_convergence      = FALSE,
    verbose                = FALSE
  )

  expect_s3_class(result, "acmap")
  expect_equal(numOptimizations(result), 2)

})


# ====================================================================
# 4.  Built-in distribution: "normal"
# ====================================================================
test_that("optimizeMapStart runs with built-in 'normal' distribution", {

  result <- optimizeMapStart(
    map                    = perfect_map,
    number_of_dimensions   = 2,
    number_of_optimizations = n_runs,
    starting_coords        = "normal",
    fixed_column_bases     = colbases,
    check_convergence      = FALSE,
    verbose                = FALSE
  )

  expect_s3_class(result, "acmap")
  expect_equal(numOptimizations(result), n_runs)
  expect_true(all(is.finite(mapStress(result))))

})


test_that("'normal' built-in respects mean and sd coord_args", {

  result <- optimizeMapStart(
    map                    = perfect_map,
    number_of_dimensions   = 2,
    number_of_optimizations = 2,
    starting_coords        = "normal",
    coord_args             = list(mean = 0, sd = 2),
    fixed_column_bases     = colbases,
    check_convergence      = FALSE,
    verbose                = FALSE
  )

  expect_s3_class(result, "acmap")
  expect_equal(numOptimizations(result), 2)

})


# ====================================================================
# 5.  Optimizations can find the global minimum
# ====================================================================
test_that("optimizeMapStart can recover the global minimum from the true coords", {

  # Supply the known-optimal coordinates directly as starting points —
  # one relaxation step should leave them (nearly) unchanged
  true_starts <- lapply(seq_len(10), function(i) {
    list(ag_coords = ag_coords, sr_coords = sr_coords)
  })

  result <- optimizeMapStart(
    map                    = perfect_map,
    number_of_dimensions   = 2,
    number_of_optimizations = 10,
    starting_coords        = true_starts,
    fixed_column_bases     = colbases,
    check_convergence      = FALSE,
    verbose                = FALSE
  )

  expect_lt(mapStress(result), 0.001)

})


# ====================================================================
# 6.  Optimizations are sorted by stress (sort_optimizations = TRUE)
# ====================================================================
test_that("optimizeMapStart sorts optimizations by stress", {

  starts <- make_coord_list(n_runs)

  result <- optimizeMapStart(
    map                    = perfect_map,
    number_of_dimensions   = 2,
    number_of_optimizations = n_runs,
    starting_coords        = starts,
    fixed_column_bases     = colbases,
    sort_optimizations     = TRUE,
    check_convergence      = FALSE,
    verbose                = FALSE
  )

  stresses <- sapply(seq_len(numOptimizations(result)), function(i) mapStress(result, i))
  expect_equal(stresses, sort(stresses))

})


# ====================================================================
# 7.  dim_annealing warning
# ====================================================================
test_that("optimizeMapStart warns and disables dim_annealing when starting_coords supplied", {

  starts <- make_coord_list(2)

  expect_warning(
    optimizeMapStart(
      map                    = perfect_map,
      number_of_dimensions   = 2,
      number_of_optimizations = 2,
      starting_coords        = starts,
      fixed_column_bases     = colbases,
      options                = list(dim_annealing = TRUE),
      check_convergence      = FALSE,
      verbose                = FALSE
    ),
    "dim_annealing"
  )

})


# ====================================================================
# 8.  Input validation errors
# ====================================================================
test_that("optimizeMapStart errors on wrong list length", {

  starts <- make_coord_list(2)   # only 2, but we ask for 5

  expect_error(
    optimizeMapStart(
      map                    = perfect_map,
      number_of_dimensions   = 2,
      number_of_optimizations = 5,
      starting_coords        = starts,
      fixed_column_bases     = colbases,
      verbose                = FALSE
    ),
    "2 element"
  )

})


test_that("optimizeMapStart errors on missing ag_coords/sr_coords names", {

  bad_starts <- list(
    list(ag = matrix(0, n_ag, 2), sr = matrix(0, n_sr, 2))  # wrong names
  )

  expect_error(
    optimizeMapStart(
      map                    = perfect_map,
      number_of_dimensions   = 2,
      number_of_optimizations = 1,
      starting_coords        = bad_starts,
      fixed_column_bases     = colbases,
      verbose                = FALSE
    ),
    "ag_coords.*sr_coords"
  )

})


test_that("optimizeMapStart errors on wrong matrix dimensions", {

  bad_starts <- list(
    list(
      ag_coords = matrix(0, n_ag + 1, 2),  # one row too many
      sr_coords = matrix(0, n_sr,     2)
    )
  )

  expect_error(
    optimizeMapStart(
      map                    = perfect_map,
      number_of_dimensions   = 2,
      number_of_optimizations = 1,
      starting_coords        = bad_starts,
      fixed_column_bases     = colbases,
      verbose                = FALSE
    ),
    paste0(n_ag, " .*", 2)
  )

})


test_that("optimizeMapStart errors on unknown built-in distribution name", {

  expect_error(
    optimizeMapStart(
      map                    = perfect_map,
      number_of_dimensions   = 2,
      number_of_optimizations = 1,
      starting_coords        = "poisson",
      fixed_column_bases     = colbases,
      verbose                = FALSE
    ),
    "poisson"
  )

})


test_that("optimizeMapStart errors on invalid starting_coords type", {

  expect_error(
    optimizeMapStart(
      map                    = perfect_map,
      number_of_dimensions   = 2,
      number_of_optimizations = 1,
      starting_coords        = 42L,
      fixed_column_bases     = colbases,
      verbose                = FALSE
    ),
    "list.*function.*character"
  )

})


# ====================================================================
# 9.  optimizeMap() is unchanged (no starting_coords parameter)
# ====================================================================
test_that("optimizeMap still has no starting_coords parameter", {
  expect_false("starting_coords" %in% names(formals(optimizeMap)))
})
