
library(Racmacs)
library(testthat)
context("Test optimization methods")

run.maptests(
  bothclasses = TRUE,
  loadlocally = FALSE,
  {

    # Create a toy HI table
    titer_table <- rbind(
      c("<10", "40",  "*",  "1280"),
      c("<10", "<10", "80", "640")
    )

    num_antigens <- nrow(titer_table)
    num_sera     <- ncol(titer_table)

    # Create a new map
    map <- make.map(
      table = titer_table
    )

    ag_names   <- agNames(map)
    sera_names <- srNames(map)

    # Add some optimizations
    minimum_column_bases <- c("640", "none")
    dimensions           <- c(2, 3)

    for(x in seq_along(minimum_column_bases)){

      minimum_column_basis <- minimum_column_bases[x]

      ag_coords <- matrix(1, num_antigens, dimensions[x])
      sr_coords <- matrix(1, num_sera,     dimensions[x])

      map <- addOptimization(
        map                  = map,
        ag_coords            = ag_coords,
        sr_coords            = sr_coords,
        minimum_column_basis = minimum_column_basis
      )

    }

    # Set the main optimization number
    selectedOptimization(map) <- 2

    # Wrong length mincolbasis
    test_that("Min column bases wrong length throws error", {
      expect_error(
        addOptimization(map       = map,
                        ag_coords = ag_coords,
                        sr_coords = sr_coords,
                        minimum_column_basis = c("1280", "none")),
        regexp = "minumum_column_basis must be provided as a vector of length 1"
      )
    })

    # Colbases alongside min col basis
    test_that("Column bases match when specified alongside min col basis", {
      expect_error(
        addOptimization(map       = map,
                        ag_coords = ag_coords,
                        sr_coords = sr_coords,
                        column_bases  = c(1,2,3,7),
                        minimum_column_basis = "none"),
        regexp = "Column bases provided do not match up with minimum_column_basis specification of 'none'"
      )
    })

    # Colbases fixed without specifying what they are
    test_that("Error when fixed col bases specified without specifying column bases", {
      expect_error(
        addOptimization(map       = map,
                        ag_coords = ag_coords,
                        sr_coords = sr_coords,
                        minimum_column_basis = "fixed"),
        regexp = "Column bases must be provided .*"
      )
    })

    # Column bases wrong length
    test_that("Error when column bases the wrong length", {
      expect_error(
        addOptimization(map       = testmap,
                        ag_coords = ag_coords,
                        sr_coords = sr_coords,
                        minimum_column_basis = "fixed",
                        column_bases = c(4,3,4,5,9))
      )
    })

    # Add an optimization correctly
    testmap <- cloneMap(map)
    testmap <- addOptimization(map       = testmap,
                               ag_coords = ag_coords,
                               sr_coords = sr_coords,
                               minimum_column_basis = "fixed",
                               column_bases = c(4,3,4,5))

    optimization_num <- numOptimizations(testmap)

    # Check column bases worked
    test_that("Fixed column bases specified correctly upon adding optimization", {
      expect_equal(
        unname(colBases(testmap, optimization_num)),
        c(4,3,4,5)
      )
    })

    # Check min colbasis changes to fixed
    test_that("Min column basis changes to fixed when column bases specified", {
      expect_equal(
        minColBasis(testmap, optimization_num),
        "fixed"
      )
    })

  }
)

