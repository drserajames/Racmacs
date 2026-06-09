
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]

#include "acmap_optimization.h"
#include "ac_stress_blobs.h"
#include "ac_stress.h"
#include "ac_relax_coords.h"
#include "ac_optimizer_options.h"

// Stress contribution of one antigen across all sera
static double calc_ag_stress(
    int ag,
    const arma::mat &ag_coords,
    const arma::mat &sr_coords,
    const arma::mat &tabledists,
    const arma::imat &titertypes,
    double dilution_stepsize
) {
  double stress = 0.0;
  for (arma::uword sr = 0; sr < sr_coords.n_rows; sr++) {
    arma::sword ttype = titertypes(ag, sr);
    if (ttype <= 0) continue;
    double map_dist = arma::norm(arma::vectorise(ag_coords.row(ag) - sr_coords.row(sr)));
    double td = tabledists(ag, sr);
    stress += ac_ptStress(map_dist, td, ttype, dilution_stepsize);
  }
  return stress;
}

// LispMDS method: for each antigen, try num_randomizations random perturbations
// and keep the position with the lowest stress (if lower than original).
static AcOptimization move_trapped_points_lispmds(
    AcOptimization optimization,
    const arma::mat &tabledists,
    const arma::imat &titertypes,
    const AcOptimizerOptions &options,
    int num_randomizations,
    double randomize_distance,
    double dilution_stepsize
) {

  int num_ags = optimization.num_ags();
  int num_sr  = optimization.num_sr();
  arma::mat ag_coords = optimization.get_ag_base_coords();
  arma::mat sr_coords = optimization.get_sr_base_coords();
  arma::uword ndims   = ag_coords.n_cols;

  // All sera are fixed throughout — only the test antigen is free
  arma::uvec fixed_sera = arma::regspace<arma::uvec>(0, num_sr - 1);

  for (int ag = 0; ag < num_ags; ag++) {

    double orig_stress = calc_ag_stress(ag, ag_coords, sr_coords, tabledists, titertypes, dilution_stepsize);

    arma::rowvec best_coords = ag_coords.row(ag);
    double best_stress = orig_stress;

    // All antigens fixed except this one
    arma::uvec fixed_ags = arma::regspace<arma::uvec>(0, num_ags - 1);
    fixed_ags.shed_row(ag);

    for (int i = 0; i < num_randomizations; i++) {

      // Perturb uniformly by up to randomize_distance in each dimension
      arma::mat ag_coords_trial = ag_coords;
      ag_coords_trial.row(ag) += (arma::randu<arma::rowvec>(ndims) * 2.0 - 1.0) * randomize_distance;

      // Relax only this antigen
      ac_relax_coords(
        tabledists, titertypes, ag_coords_trial, sr_coords,
        options, fixed_ags, fixed_sera, arma::mat(), dilution_stepsize
      );

      double trial_stress = calc_ag_stress(ag, ag_coords_trial, sr_coords, tabledists, titertypes, dilution_stepsize);

      if (trial_stress < best_stress) {
        best_stress   = trial_stress;
        best_coords   = ag_coords_trial.row(ag);
      }
    }

    if (best_stress < orig_stress) {
      ag_coords.row(ag) = best_coords;
    }
  }

  optimization.set_ag_base_coords(ag_coords);
  return optimization;
}

// Check for trapped antigens
arma::mat check_ag_trapped_points(
    const AcOptimization &optimization,
    const arma::mat &tabledists,
    const arma::imat &titertypes,
    const double &grid_spacing,
    AcOptimizerOptions options
){

  // Variables
  double stress_lim = 0;
  int num_ags = optimization.num_ags();
  arma::mat ag_coords = optimization.get_ag_base_coords();
  arma::mat sr_coords = optimization.get_sr_base_coords();

  arma::mat trapped_ag_improved_coords(arma::size(ag_coords));
  trapped_ag_improved_coords.fill(arma::datum::nan);

  // Check trapped antigens
  #pragma omp parallel for schedule(dynamic) num_threads(options.num_cores)
  for(int ag=0; ag<num_ags; ag++){

    // Do a grid search
    StressBlobGrid grid_results = ac_stress_blob_grid(
      ag_coords.row(ag).as_col(),
      sr_coords,
      tabledists.row(ag).as_col(),
      titertypes.row(ag).as_col(),
      stress_lim,
      grid_spacing
    );

    // Check if any grid points have lower stress than the minimum
    if(grid_results.grid.min() < 0.0){
      arma::uword index = grid_results.grid.index_min();
      arma::uvec sub = arma::ind2sub( arma::size(grid_results.grid), index );
      trapped_ag_improved_coords(ag,0) = grid_results.xcoords( sub(0) );
      trapped_ag_improved_coords(ag,1) = grid_results.ycoords( sub(1) );
    }

  }

  // Return trapped point information
  return trapped_ag_improved_coords;

}


// Check for trapped sera
arma::mat check_sr_trapped_points(
    const AcOptimization &optimization,
    const arma::mat &tabledists,
    const arma::imat &titertypes,
    const double &grid_spacing,
    AcOptimizerOptions options
){

  // Variables
  double stress_lim = 0;
  int num_sr = optimization.num_sr();
  arma::mat ag_coords = optimization.get_ag_base_coords();
  arma::mat sr_coords = optimization.get_sr_base_coords();

  arma::mat trapped_sr_improved_coords(arma::size(sr_coords));
  trapped_sr_improved_coords.fill(arma::datum::nan);

  // Check trapped sera
  #pragma omp parallel for schedule(dynamic) num_threads(options.num_cores)
  for(int sr=0; sr<num_sr; sr++){

    // Do a grid search
    StressBlobGrid grid_results = ac_stress_blob_grid(
      sr_coords.row(sr).as_col(),
      ag_coords,
      tabledists.col(sr),
      titertypes.col(sr),
      stress_lim,
      grid_spacing
    );

    // Check if any grid points have lower stress than the minimum
    if(grid_results.grid.min() < 0.0){
      arma::uword index = grid_results.grid.index_min();
      arma::uvec sub = arma::ind2sub( arma::size(grid_results.grid), index );
      trapped_sr_improved_coords(sr,0) = grid_results.xcoords( sub(0) );
      trapped_sr_improved_coords(sr,1) = grid_results.ycoords( sub(1) );
    }

  }

  // Return trapped point information
  return trapped_sr_improved_coords;

}


// Function to find and move trapped coordinates
// [[Rcpp::export]]
AcOptimization ac_move_trapped_points(
  AcOptimization optimization,
  AcTiterTable titertable,
  double grid_spacing,
  AcOptimizerOptions options,
  int max_iterations = 10,
  double dilution_stepsize = 1.0,
  std::string method = "racmacs",
  int num_randomizations = 10,
  double randomize_distance = 20.0
){


  // Dispatch to LispMDS method if requested
  if (method == "lispmds") {
    arma::imat titertypes_l = titertable.get_titer_types();
    arma::mat tabledists_l = titertable.numeric_table_distances(
      optimization.get_min_column_basis(),
      optimization.get_fixed_column_bases(),
      optimization.get_ag_reactivity_adjustments()
    );
    return move_trapped_points_lispmds(
      optimization, tabledists_l, titertypes_l, options,
      num_randomizations, randomize_distance, dilution_stepsize
    );
  }

  // Check antigen and sera trapped points recursively
  if(options.report_progress) REprintf("Checking for trapped points recursively:");

  arma::imat titertypes = titertable.get_titer_types();
  arma::mat tabledists = titertable.numeric_table_distances(
    optimization.get_min_column_basis(),
    optimization.get_fixed_column_bases(),
    optimization.get_ag_reactivity_adjustments()
  );

  int num_iterations = 0;
  while(num_iterations < max_iterations){

    // Variables
    arma::mat ag_coords = optimization.get_ag_base_coords();
    arma::mat sr_coords = optimization.get_sr_base_coords();

    // Check for any improved coordinates
    arma::mat ag_trapped_improved_coords = check_ag_trapped_points(optimization, tabledists, titertypes, grid_spacing, options);
    arma::mat sr_trapped_improved_coords = check_sr_trapped_points(optimization, tabledists, titertypes, grid_spacing, options);

    // Get any improved indices
    arma::uvec ag_trapped_coord_indices = arma::find_finite(ag_trapped_improved_coords);
    arma::uvec sr_trapped_coord_indices = arma::find_finite(sr_trapped_improved_coords);

    // Break if no improvements found
    if((ag_trapped_coord_indices.n_elem == 0) && (sr_trapped_coord_indices.n_elem == 0)){
      break;
    }

    // Move antigen and serum coordinates to improved positions
    ag_coords.elem(ag_trapped_coord_indices) = ag_trapped_improved_coords.elem(ag_trapped_coord_indices);
    sr_coords.elem(sr_trapped_coord_indices) = sr_trapped_improved_coords.elem(sr_trapped_coord_indices);

    optimization.set_ag_base_coords(ag_coords);
    optimization.set_sr_base_coords(sr_coords);

    // Relax the optimization
    optimization.relax_from_raw_matrices(
      tabledists,
      titertypes,
      options,
      arma::uvec(),
      arma::uvec(),
      arma::mat(),
      dilution_stepsize
    );

    // Increment loop num
    if(options.report_progress) REprintf(".");
    num_iterations++;

  }

  // Output message indicating if some were found
  if(options.report_progress){
    if(num_iterations == 0){
      REprintf(" no trapped points found.\n");
    } else if(num_iterations == max_iterations){
      REprintf(" maximum iteration number reached.\n");
    } else {
      REprintf(" all trapped points moved.\n");
    }
  }

  // Return the improved optimization
  return optimization;

}

