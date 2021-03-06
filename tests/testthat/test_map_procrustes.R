
library(Racmacs)
library(testthat)
context("Test procrustes methods")

# Test for raccharts and racmaps
run.maptests(
  bothclasses = TRUE,
  loadlocally = FALSE,
  {

    # Setup rotation and translation matrices
    rot_mat <- matrix(data = c(cos(0.24), sin(0.24), -sin(0.24), cos(0.24)),
                      nrow = 2,
                      ncol = 2)

    inv_rot_mat <- t(rot_mat)

    coords1 <- matrix(runif(10), ncol = 2)
    coords2 <- coords1%*%rot_mat
    coords3 <- coords2%*%inv_rot_mat

    trans_mat <- matrix(c(2.4,3.8), nrow = 1)


    # Create new maps
    num_ags <- c(10, 14)
    num_sr  <- c(8,  6)

    ag_names1 <- paste("Map1 antigen", seq_len(num_ags[1]))
    sr_names1 <- paste("Map1 sera",    seq_len(num_sr[1]))

    ag_names2 <- paste("Map2 antigen", seq_len(num_ags[2]))
    sr_names2 <- paste("Map2 sera",    seq_len(num_sr[2]))


    # Define matching antigens
    matching_ags <- rbind(
      c(10, 2),
      c(1,  14),
      c(3,  4),
      c(6,  7),
      c(4,  11)
    )

    matching_sr <- rbind(
      c(3, 1),
      c(1, 3),
      c(8, 6),
      c(2, 5)
    )

    for(x in seq_len(nrow(matching_ags))){
      ag_names1[matching_ags[x,1]] <- paste("Matched antigen", x)
      ag_names2[matching_ags[x,2]] <- paste("Matched antigen", x)
    }

    for(x in seq_len(nrow(matching_sr))){
      sr_names1[matching_sr[x,1]] <- paste("Matched sera", x)
      sr_names2[matching_sr[x,2]] <- paste("Matched sera", x)
    }


    # Generate coordinates
    ag_coords1 <- matrix(runif(num_ags[1]*2)*10, ncol = 2)
    ag_coords2 <- matrix(runif(num_ags[2]*2)*10, ncol = 2)
    rownames(ag_coords1) <- ag_names1
    rownames(ag_coords2) <- ag_names2

    sr_coords1 <- matrix(runif(num_sr[1]*2)*10, ncol = 2)
    sr_coords2 <- matrix(runif(num_sr[2]*2)*10, ncol = 2)
    rownames(sr_coords1) <- sr_names1
    rownames(sr_coords2) <- sr_names2


    # Create the test maps
    map1 <- make.map(ag_coords = ag_coords1,
                     sr_coords = sr_coords1,
                     ag_names  = ag_names1,
                     sr_names  = sr_names1,
                     minimum_column_basis = "none")

    map2 <- make.map(ag_coords = ag_coords2,
                     sr_coords = sr_coords2,
                     ag_names  = ag_names2,
                     sr_names  = sr_names2,
                     minimum_column_basis = "none")

    # Create a rotated and shuffled version
    ## Shuffle antigens and sera
    ag_order1rot <- sample(seq_along(ag_names1))
    sr_order1rot <- sample(seq_along(sr_names1))

    ## Add some name mismatches
    ag_mismatches1rot <- c(1,6)
    sr_mismatches1rot <- c(2,3)
    ag_names1rot <- ag_names1[ag_order1rot]
    sr_names1rot <- sr_names1[sr_order1rot]
    ag_names1rot[ag_mismatches1rot] <- paste("Mismatched antigen", ag_mismatches1rot)
    sr_names1rot[sr_mismatches1rot] <- paste("Mismatched sera",    sr_mismatches1rot)

    ## Rotate and translate the coordinates
    ag_coords1rot <- ag_coords1%*%rot_mat + matrix(trans_mat, num_ags[1], 2, byrow = TRUE)
    sr_coords1rot <- sr_coords1%*%rot_mat + matrix(trans_mat, num_sr[1],  2, byrow = TRUE)
    ag_coords1rot <- ag_coords1rot[ag_order1rot,]
    sr_coords1rot <- sr_coords1rot[sr_order1rot,]

    map1rot <- make.map(ag_coords = ag_coords1rot,
                        sr_coords = sr_coords1rot,
                        ag_names  = ag_names1rot,
                        sr_names  = sr_names1rot,
                        minimum_column_basis = "none")


    # Test procrustes of map to itself
    test_that("Realign a map to itself", {

      omap1 <- cloneMap(map1)
      omap1 <- realignMap(omap1, map1)
      expect_equal(agCoords(omap1), agCoords(map1))
      expect_equal(srCoords(omap1), srCoords(map1))

    })


    test_that("Procrustes a map to itself", {

      pc1 <- procrustesMap(map1, map1)
      expected_ag_dists <- rep(0, num_ags[1])
      expected_sr_dists <- rep(0, num_sr[1])

      expect_equal(pc1$procrustes$ag_dists, expected_ag_dists)
      expect_equal(pc1$procrustes$sr_dists, expected_sr_dists)
      expect_equal(pc1$procrustes$ag_rmsd, 0)
      expect_equal(pc1$procrustes$sr_rmsd, 0)
      expect_equal(pc1$procrustes$pc_coords$ag, unname(ag_coords1))
      expect_equal(pc1$procrustes$pc_coords$sr, unname(sr_coords1))

    })

    test_that("Realign to a transformed version", {

      omap1 <- cloneMap(map1)
      expect_warning(omap1 <- realignMap(omap1, map1rot))
      expect_equal(unname(srCoords(omap1)[sr_order1rot,]), unname(srCoords(map1rot)))
      expect_equal(unname(agCoords(omap1)[ag_order1rot,]), unname(agCoords(map1rot)))

    })

    test_that("Procrustes to a transformed version", {

      expect_warning(pc1 <- procrustesMap(map1rot, map1))

      expected_ag_dists <- rep(0, num_ags[1])
      expected_sr_dists <- rep(0, num_sr[1])
      expected_ag_dists[ag_mismatches1rot] <- NA
      expected_sr_dists[sr_mismatches1rot] <- NA

      expected_pc_coords_ag <- agCoords(map1rot)
      expected_pc_coords_sr <- srCoords(map1rot)
      expected_pc_coords_ag[ag_mismatches1rot,] <- NA
      expected_pc_coords_sr[sr_mismatches1rot,] <- NA

      expect_equal(round(pc1$procrustes$ag_dists, 5), expected_ag_dists)
      expect_equal(round(pc1$procrustes$sr_dists, 5), expected_sr_dists)
      expect_equal(round(pc1$procrustes$ag_rmsd, 5), 0)
      expect_equal(round(pc1$procrustes$sr_rmsd, 5), 0)
      expect_equal(unname(pc1$procrustes$pc_coords$ag), unname(expected_pc_coords_ag))
      expect_equal(unname(pc1$procrustes$pc_coords$sr), unname(expected_pc_coords_sr))

    })


    # Convert to chart so we can compare against acmacs
    chart1 <- as.cpp(map1)
    chart2 <- as.cpp(map2)

    test_that("Procrustes gets same as acmacs", {

      ## Perform Racmacs procrustes
      expect_warning({ racmacs_pc              <- procrustesMap(map1, map2) })
      expect_warning({ racmacs_pc_ags          <- procrustesMap(map1, map2, sera     = FALSE) })
      expect_warning({ racmacs_pc_sr           <- procrustesMap(map1, map2, antigens = FALSE) })
      expect_warning({ racmacs_pc_scaling      <- procrustesMap(map1, map2, scaling = TRUE) })
      expect_warning({ racmacs_pc_ags_scaling  <- procrustesMap(map1, map2, scaling = TRUE, sera     = FALSE) })
      expect_warning({ racmacs_pc_sr_scaling   <- procrustesMap(map1, map2, scaling = TRUE, antigens = FALSE) })
      expect_warning({ racmacs_pc_rotation     <- procrustesMap(map1, map2, translation = FALSE) })
      expect_warning({ racmacs_pc_ags_rotation <- procrustesMap(map1, map2, translation = FALSE, sera     = FALSE) })
      expect_warning({ racmacs_pc_sr_rotation  <- procrustesMap(map1, map2, translation = FALSE, antigens = FALSE) })

      ## Perform acmacs procrustes
      optimization1 <- chart1$chart$projections[[1]]
      optimization2 <- chart2$chart$projections[[1]]
      acmacs_pc             <- acmacs.r::acmacs.procrustes(optimization1, optimization2, scaling = FALSE, match = "ignore")
      acmacs_pc_ags         <- acmacs.r::acmacs.procrustes(optimization1, optimization2, scaling = FALSE, match = "ignore", subset = "antigens")
      acmacs_pc_sr          <- acmacs.r::acmacs.procrustes(optimization1, optimization2, scaling = FALSE, match = "ignore", subset = "sera")
      acmacs_pc_scaling     <- acmacs.r::acmacs.procrustes(optimization1, optimization2, scaling = TRUE, match = "ignore")
      acmacs_pc_ags_scaling <- acmacs.r::acmacs.procrustes(optimization1, optimization2, scaling = TRUE, match = "ignore", subset = "antigens")
      acmacs_pc_sr_scaling  <- acmacs.r::acmacs.procrustes(optimization1, optimization2, scaling = TRUE, match = "ignore", subset = "sera")

      ## Check equivalence
      ### rmsd
      expect_equal(racmacs_pc$procrustes$total_rmsd,             acmacs_pc$rms             )
      expect_equal(racmacs_pc_ags$procrustes$total_rmsd,         acmacs_pc_ags$rms         )
      expect_equal(racmacs_pc_sr$procrustes$total_rmsd,          acmacs_pc_sr$rms          )
      expect_equal(racmacs_pc_scaling$procrustes$total_rmsd,     acmacs_pc_scaling$rms     )
      expect_equal(racmacs_pc_ags_scaling$procrustes$total_rmsd, acmacs_pc_ags_scaling$rms )
      expect_equal(racmacs_pc_sr_scaling$procrustes$total_rmsd,  acmacs_pc_sr_scaling$rms  )

      ### Transformed coordinate positions
      acmacs_matched_optimization <- function(agnames1,
                                              agnames2,
                                              srnames1,
                                              srnames2,
                                              optimization1,
                                              optimization2,
                                              point_type,
                                              acmacs_pc){
        # Rotation
        coords <- optimization2$layout %*% acmacs_pc$transformation[1:2, 1:2]
        # Translation
        coords <- coords + matrix(acmacs_pc$transformation[3,1:2], num_ags[2] + num_sr[2], ncol = 2, byrow = TRUE)
        # Antigen and sera coords
        ag_coords <- coords[seq_along(agnames2),,drop=FALSE]
        sr_coords <- coords[-seq_along(agnames2),,drop=FALSE]
        # Name matching to 1
        ag_coords <- ag_coords[match(agnames1, agnames2),,drop=FALSE]
        sr_coords <- sr_coords[match(srnames1, srnames2),,drop=FALSE]
        # Drop unmatched
        if(point_type == "antigens") sr_coords[] <- NA
        if(point_type == "sera")    ag_coords[] <- NA
        # Return list
        list(ag = ag_coords,
             sr = sr_coords)
      }

      expect_equal(racmacs_pc$procrustes$pc_coords[c("ag", "sr")],     acmacs_matched_optimization(ag_names1, ag_names2, sr_names1, sr_names2, optimization1, optimization2, "both", acmacs_pc))
      expect_equal(racmacs_pc_ags$procrustes$pc_coords[c("ag", "sr")], acmacs_matched_optimization(ag_names1, ag_names2, sr_names1, sr_names2, optimization1, optimization2, "antigens", acmacs_pc_ags))
      expect_equal(racmacs_pc_sr$procrustes$pc_coords[c("ag", "sr")],  acmacs_matched_optimization(ag_names1, ag_names2, sr_names1, sr_names2, optimization1, optimization2, "sera", acmacs_pc_sr))
      expect_equal(racmacs_pc_scaling$procrustes$pc_coords[c("ag", "sr")],     acmacs_matched_optimization(ag_names1, ag_names2, sr_names1, sr_names2, optimization1, optimization2, "both", acmacs_pc_scaling))
      expect_equal(racmacs_pc_ags_scaling$procrustes$pc_coords[c("ag", "sr")], acmacs_matched_optimization(ag_names1, ag_names2, sr_names1, sr_names2, optimization1, optimization2, "antigens", acmacs_pc_ags_scaling))
      expect_equal(racmacs_pc_sr_scaling$procrustes$pc_coords[c("ag", "sr")],  acmacs_matched_optimization(ag_names1, ag_names2, sr_names1, sr_names2, optimization1, optimization2, "sera", acmacs_pc_sr_scaling))

    })


    # Realign a map that's been rotated into 3D
    test_that("Realigning 2D to 3D and back", {

      coords2d <- matrix(c(2,3,1,8,3,3,2,9,1,0), 5, 2)
      coords3d <- coords2d %*% rotation_matrix_3D(1.2, "y")[1:2,]

      map2d <- make.map(
        ag_coords = coords2d[1:3,],
        sr_coords = coords2d[4:5,],
        minimum_column_basis = "none"
      )

      map3d <- make.map(
        ag_coords = coords3d[1:3,],
        sr_coords = coords3d[4:5,],
        minimum_column_basis = "none"
      )

      pc2d3d <- procrustesMap(
        map2d,
        map3d
      )

      pc3d2d <- procrustesMap(
        map3d,
        map2d
      )

      expect_equal(round(pc2d3d$procrustes$total_rmsd, 5), 0)
      expect_equal(round(pc3d2d$procrustes$total_rmsd, 5), 0)

    })


    # Realign a map that's been rotated into 3D
    test_that("Realigning 2D to 3D and back in a rotated map", {

      coords2d <- matrix(c(2,3,1,8,3,3,2,9,1,0), 5, 2)

      map2d <- make.map(
        ag_coords = coords2d[1:3,],
        sr_coords = coords2d[4:5,],
        minimum_column_basis = "none"
      )

      map3d <- make.map(
        ag_coords = coords2d[1:3,],
        sr_coords = coords2d[4:5,],
        minimum_column_basis = "none"
      )

      mapTransformation(map3d) <- rotation_matrix_3D(1, "x")
      mapTranslation(map3d)    <- c(1,2)

      pc2d3d <- procrustesMap(
        map2d,
        map3d
      )

      pc3d2d <- procrustesMap(
        map3d,
        map2d
      )

      map3d <- realignMap(
        map3d,
        map2d
      )

      map2d <- realignMap(
        map2d,
        map3d
      )

      expect_equal(pc2d3d$procrustes$total_rmsd, 0)
      expect_equal(pc3d2d$procrustes$total_rmsd, 0)
      expect_equal(mapTransformation(map3d), diag(nrow = 3))
      expect_equal(mapTranslation(map3d), matrix(0, nrow = 1, ncol = 3))
      expect_equal(ncol(agBaseCoords(map3d)), 2)

    })


    # Testing realigning optimizations
    mapA <- read.map(test_path("../testdata/testmap.ace"))
    mapB <- cloneMap(mapA)

    test_that("Realigning map optimizations 3D to 2D", {

      mapB <- realignOptimizations(mapB)

      expect_lt(
        sum((agCoords(mapB, 1) - MCMCpack::procrustes(agCoords(mapB, 2), agCoords(mapB, 1))$X.new)^2),
        sum((agCoords(mapA, 1) - MCMCpack::procrustes(agCoords(mapA, 2), agCoords(mapA, 1))$X.new)^2)
      )

      expect_lt(
        sum((cbind(agCoords(mapB, 1), 0) - MCMCpack::procrustes(agCoords(mapB, 3), cbind(agCoords(mapB, 1), 0))$X.new)^2),
        sum((cbind(agCoords(mapA, 1), 0) - MCMCpack::procrustes(agCoords(mapA, 3), cbind(agCoords(mapA, 1), 0))$X.new)^2)
      )

    })


    mapA <- read.map(test_path("../testdata/testmap.ace"))
    mapB <- cloneMap(mapA)

    test_that("Realigning map optimizations 2D to 3D", {

      selectedOptimization(mapB) <- 3
      mapB <- realignOptimizations(mapB)

      pcA <- procrustesMap(
        map = mapA,
        comparison_map = mapA,
        optimization_number = 1,
        comparison_optimization_number = 3
      )

      pcB <- procrustesMap(
        map = mapB,
        comparison_map = mapB,
        optimization_number = 3,
        comparison_optimization_number = 1
      )

      export.viewer.test(
        view(
          pcB
        ),
        "procrustes_3d_to_2d.html"
      )

      expect_equal(pcA$total_rmsd, pcB$total_rmsd)

      expect_lt(
        sum(calc_coord_dist(agCoords(mapB, 1), agCoords(mapB, 2))^2),
        sum(calc_coord_dist(agCoords(mapA, 1), agCoords(mapA, 2))^2)
      )

    })

    test_that("Procrustes maps with na coords", {

      map1na <- cloneMap(map1)
      agCoords(map1na)[1:2,] <- NA
      srCoords(map1na)[1,] <- NA
      export.viewer.test(
        view(map1na),
        "na_map.html"
      )

      expect_warning({
        pcmap <- procrustesMap(
          map1na,
          map2
        )
      })
      export.viewer.test(
        view(
          pcmap
        ),
        "na_map_procrustes.html"
      )

    })

})


