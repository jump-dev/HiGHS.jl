# HiGHS.jl

[![Build Status](https://github.com/jump-dev/HiGHS.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/HiGHS.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/HiGHS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/HiGHS.jl)

HiGHS.jl is a wrapper for the [HiGHS](https://highs.dev) linear solver.

It has two components:
 - a thin wrapper around the complete C API
 - an interface to [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)

The C API can be accessed via `HiGHS.Highs_xxx` functions, where the names and
arguments are identical to the C API.

## Installation

Install HiGHS as follows:
```julia
import Pkg
Pkg.add("HiGHS")
```

In addition to installing the HiGHS.jl package, this will also download and
install the HiGHS binaries. (You do not need to install HiGHS separately.)

To use a custom binary, read the [Custom solver binaries](https://jump.dev/JuMP.jl/stable/developers/custom_solver_binaries/)
section of the JuMP documentation.

## Issues and feedback

To report a problem (e.g., incorrect results, or a crash of the solver),
or make a suggestion for how to improve HiGHS, please [file a GitHub issue](https://github.com/jump-dev/HiGHS.jl).

If you use HiGHS from JuMP, use `JuMP.write_to_file(model, "filename.mps")`
to write your model an MPS file, then upload the MPS file to [https://gist.github.com](https://gist.github.com)
and provide a link to the gist in the GitHub issue.

## Use with JuMP

Pass `HiGHS.Optimizer` to `JuMP.Model` to create a JuMP model with HiGHS as the
optimizer. Set options using `set_optimizer_attribute`.

```julia
using JuMP
import HiGHS
model = Model(HiGHS.Optimizer)
set_optimizer_attribute(model, "presolve", "on")
set_optimizer_attribute(model, "time_limit", 60.0)
```

## Options

To print the list of options supported by HiGHS, do:
```julia
using HiGHS
model = HiGHS.Optimizer()
Highs_writeOptions(model, "options.txt")
println(read("options.txt", String))
```

The current option list is:
```
# Presolve option: "off", "choose" or "on"
# [type: string, advanced: false, default: "choose"]
presolve = choose

# Solver option: "simplex", "choose" or "ipm"
# [type: string, advanced: false, default: "choose"]
solver = choose

# Parallel option: "off", "choose" or "on"
# [type: string, advanced: false, default: "choose"]
parallel = choose

# Time limit
# [type: double, advanced: false, range: [0, inf], default: inf]
time_limit = inf

# Compute cost, bound, RHS and basic solution ranging: "off" or "on"
# [type: string, advanced: false, default: "off"]
ranging = off

# Limit on cost coefficient: values larger than this will be treated as infinite
# [type: double, advanced: false, range: [1e+15, inf], default: 1e+20]
infinite_cost = 1e+20

# Limit on |constraint bound|: values larger than this will be treated as infinite
# [type: double, advanced: false, range: [1e+15, inf], default: 1e+20]
infinite_bound = 1e+20

# Lower limit on |matrix entries|: values smaller than this will be treated as zero
# [type: double, advanced: false, range: [1e-12, inf], default: 1e-09]
small_matrix_value = 1e-09

# Upper limit on |matrix entries|: values larger than this will be treated as infinite
# [type: double, advanced: false, range: [1, inf], default: 1e+15]
large_matrix_value = 1e+15

# Primal feasibility tolerance
# [type: double, advanced: false, range: [1e-10, inf], default: 1e-07]
primal_feasibility_tolerance = 1e-07

# Dual feasibility tolerance
# [type: double, advanced: false, range: [1e-10, inf], default: 1e-07]
dual_feasibility_tolerance = 1e-07

# IPM optimality tolerance
# [type: double, advanced: false, range: [1e-12, inf], default: 1e-08]
ipm_optimality_tolerance = 1e-08

# Objective bound for termination
# [type: double, advanced: false, range: [-inf, inf], default: inf]
objective_bound = inf

# Objective target for termination
# [type: double, advanced: false, range: [-inf, inf], default: -inf]
objective_target = -inf

# random seed used in HiGHS
# [type: HighsInt, advanced: false, range: {0, 2147483647}, default: 0]
random_seed = 0

# number of threads used by HiGHS (0: automatic)
# [type: HighsInt, advanced: false, range: {0, 2147483647}, default: 0]
threads = 0

# Debugging level in HiGHS
# [type: HighsInt, advanced: false, range: {0, 3}, default: 0]
highs_debug_level = 0

# Analysis level in HiGHS
# [type: HighsInt, advanced: false, range: {0, 63}, default: 0]
highs_analysis_level = 0

# Strategy for simplex solver 0 => Choose; 1 => Dual (serial); 2 => Dual (PAMI); 3 => Dual (SIP); 4 => Primal
# [type: HighsInt, advanced: false, range: {0, 4}, default: 1]
simplex_strategy = 1

# Simplex scaling strategy: off / choose / equilibration / forced equilibration / max value 0 / max value 1 (0/1/2/3/4/5)
# [type: HighsInt, advanced: false, range: {0, 5}, default: 1]
simplex_scale_strategy = 1

# Strategy for simplex crash: off / LTSSF / Bixby (0/1/2)
# [type: HighsInt, advanced: false, range: {0, 9}, default: 0]
simplex_crash_strategy = 0

# Strategy for simplex dual edge weights: Choose / Dantzig / Devex / Steepest Edge (-1/0/1/2)
# [type: HighsInt, advanced: false, range: {-1, 2}, default: -1]
simplex_dual_edge_weight_strategy = -1

# Strategy for simplex primal edge weights: Choose / Dantzig / Devex / Steepest Edge (-1/0/1/2)
# [type: HighsInt, advanced: false, range: {-1, 2}, default: -1]
simplex_primal_edge_weight_strategy = -1

# Iteration limit for simplex solver
# [type: HighsInt, advanced: false, range: {0, 2147483647}, default: 2147483647]
simplex_iteration_limit = 2147483647

# Limit on the number of simplex UPDATE operations
# [type: HighsInt, advanced: false, range: {0, 2147483647}, default: 5000]
simplex_update_limit = 5000

# Iteration limit for IPM solver
# [type: HighsInt, advanced: false, range: {0, 2147483647}, default: 2147483647]
ipm_iteration_limit = 2147483647

# Minimum level of concurrency in parallel simplex
# [type: HighsInt, advanced: false, range: {1, 8}, default: 1]
simplex_min_concurrency = 1

# Maximum level of concurrency in parallel simplex
# [type: HighsInt, advanced: false, range: {1, 8}, default: 8]
simplex_max_concurrency = 8

# Enables or disables solver output
# [type: bool, advanced: false, range: {false, true}, default: true]
output_flag = true

# Enables or disables console logging
# [type: bool, advanced: false, range: {false, true}, default: true]
log_to_console = true

# Solution file
# [type: string, advanced: false, default: ""]
solution_file = ""

# Log file
# [type: string, advanced: false, default: ""]
log_file = ""

# Write the primal and dual solution to a file
# [type: bool, advanced: false, range: {false, true}, default: false]
write_solution_to_file = false

# Write the solution in style: 0=>Raw (computer-readable); 1=>Pretty (human-readable) 
# [type: HighsInt, advanced: false, range: {0, 2}, default: 0]
write_solution_style = 0

# Write model file
# [type: string, advanced: false, default: ""]
write_model_file = "" 

# Write the model to a file
# [type: bool, advanced: false, range: {false, true}, default: false]
write_model_to_file = false

# Whether symmetry should be detected
# [type: bool, advanced: false, range: {false, true}, default: true]
mip_detect_symmetry = true

# MIP solver max number of nodes
# [type: HighsInt, advanced: false, range: {0, 2147483647}, default: 2147483647]
mip_max_nodes = 2147483647

# MIP solver max number of nodes where estimate is above cutoff bound
# [type: HighsInt, advanced: false, range: {0, 2147483647}, default: 2147483647]
mip_max_stall_nodes = 2147483647

# MIP solver max number of leave nodes
# [type: HighsInt, advanced: false, range: {0, 2147483647}, default: 2147483647]
mip_max_leaves = 2147483647

# limit on the number of improving solutions found to stop the MIP solver prematurely
# [type: HighsInt, advanced: false, range: {1, 2147483647}, default: 2147483647]
mip_max_improving_sols = 2147483647

# maximal age of dynamic LP rows before they are removed from the LP relaxation
# [type: HighsInt, advanced: false, range: {0, 32767}, default: 10]
mip_lp_age_limit = 10

# maximal age of rows in the cutpool before they are deleted
# [type: HighsInt, advanced: false, range: {0, 1000}, default: 30]
mip_pool_age_limit = 30

# soft limit on the number of rows in the cutpool for dynamic age adjustment
# [type: HighsInt, advanced: false, range: {1, 2147483647}, default: 10000]
mip_pool_soft_limit = 10000

# minimal number of observations before pseudo costs are considered reliable
# [type: HighsInt, advanced: false, range: {0, 2147483647}, default: 8]
mip_pscost_minreliable = 8

# minimal number of entries in the cliquetable before neighborhood queries of the conflict graph use parallel processing
# [type: HighsInt, advanced: false, range: {0, 2147483647}, default: 100000]
mip_min_cliquetable_entries_for_parallelism = 100000

# MIP solver reporting level
# [type: HighsInt, advanced: false, range: {0, 2}, default: 1]
mip_report_level = 1

# MIP feasibility tolerance
# [type: double, advanced: false, range: [1e-10, inf], default: 1e-06]
mip_feasibility_tolerance = 1e-06

# effort spent for MIP heuristics
# [type: double, advanced: false, range: [0, 1], default: 0.05]
mip_heuristic_effort = 0.05

# tolerance on relative gap, |ub-lb|/|ub|, to determine whether optimality has been reached for a MIP instance
# [type: double, advanced: false, range: [0, inf], default: 0.0001]
mip_rel_gap = 0.0001

# tolerance on absolute gap of MIP, |ub-lb|, to determine whether optimality has been reached for a MIP instance
# [type: double, advanced: false, range: [0, inf], default: 1e-06]
mip_abs_gap = 1e-06

# Output development messages: 0 => none; 1 => info; 2 => verbose
# [type: HighsInt, advanced: true, range: {0, 3}, default: 0]
log_dev_level = 0

# Run the crossover routine for IPX
# [type: bool, advanced: true, range: {false, true}, default: true]
run_crossover = true

# Allow ModelStatus::kUnboundedOrInfeasible
# [type: bool, advanced: true, range: {false, true}, default: false]
allow_unbounded_or_infeasible = false

# Use relaxed implied bounds from presolve
# [type: bool, advanced: true, range: {false, true}, default: false]
use_implied_bounds_from_presolve = false

# Prevents LP presolve steps for which postsolve cannot maintain a basis
# [type: bool, advanced: true, range: {false, true}, default: true]
lp_presolve_requires_basis_postsolve = true

# Use the free format MPS file reader
# [type: bool, advanced: true, range: {false, true}, default: true]
mps_parser_type_free = true

# For multiple N-rows in MPS files: delete rows / delete entries / keep rows (-1/0/1)
# [type: HighsInt, advanced: true, range: {-1, 1}, default: -1]
keep_n_rows = -1

# Scaling factor for costs
# [type: HighsInt, advanced: true, range: {-20, 20}, default: 0]
cost_scale_factor = 0

# Largest power-of-two factor permitted when scaling the constraint matrix
# [type: HighsInt, advanced: true, range: {0, 30}, default: 20]
allowed_matrix_scale_factor = 20

# Largest power-of-two factor permitted when scaling the costs
# [type: HighsInt, advanced: true, range: {0, 20}, default: 0]
allowed_cost_scale_factor = 0

# Strategy for dualising before simplex
# [type: HighsInt, advanced: true, range: {-1, 1}, default: -1]
simplex_dualise_strategy = -1

# Strategy for permuting before simplex
# [type: HighsInt, advanced: true, range: {-1, 1}, default: -1]
simplex_permute_strategy = -1

# Max level of dual simplex cleanup
# [type: HighsInt, advanced: true, range: {0, 2147483647}, default: 1]
max_dual_simplex_cleanup_level = 1

# Max level of dual simplex phase 1 cleanup
# [type: HighsInt, advanced: true, range: {0, 2147483647}, default: 2]
max_dual_simplex_phase1_cleanup_level = 2

# Strategy for PRICE in simplex
# [type: HighsInt, advanced: true, range: {0, 3}, default: 3]
simplex_price_strategy = 3

# Strategy for solving unscaled LP in simplex
# [type: HighsInt, advanced: true, range: {0, 2}, default: 1]
simplex_unscaled_solution_strategy = 1

# Perform initial basis condition check in simplex
# [type: bool, advanced: true, range: {false, true}, default: true]
simplex_initial_condition_check = true

# No unnecessary refactorization on simplex rebuild
# [type: bool, advanced: true, range: {false, true}, default: true]
no_unnecessary_rebuild_refactor = true

# Tolerance on initial basis condition in simplex
# [type: double, advanced: true, range: [1, inf], default: 1e+14]
simplex_initial_condition_tolerance = 1e+14

# Tolerance on solution error when considering refactorization on simplex rebuild
# [type: double, advanced: true, range: [-inf, inf], default: 1e-08]
rebuild_refactor_solution_error_tolerance = 1e-08

# Tolerance on dual steepest edge weight errors
# [type: double, advanced: true, range: [0, inf], default: inf]
dual_steepest_edge_weight_error_tolerance = inf

# Threshold on dual steepest edge weight errors for Devex switch
# [type: double, advanced: true, range: [1, inf], default: 10]
dual_steepest_edge_weight_log_error_threshold = 10.0

# Dual simplex cost perturbation multiplier: 0 => no perturbation
# [type: double, advanced: true, range: [0, inf], default: 1]
dual_simplex_cost_perturbation_multiplier = 1.0

# Primal simplex bound perturbation multiplier: 0 => no perturbation
# [type: double, advanced: true, range: [0, inf], default: 1]
primal_simplex_bound_perturbation_multiplier = 1.0

# Dual simplex pivot growth tolerance
# [type: double, advanced: true, range: [1e-12, inf], default: 1e-09]
dual_simplex_pivot_growth_tolerance = 1e-09

# Matrix factorization pivot threshold for substitutions in presolve
# [type: double, advanced: true, range: [0.0008, 0.5], default: 0.01]
presolve_pivot_threshold = 0.01

# Maximal fillin allowed for substitutions in presolve
# [type: HighsInt, advanced: true, range: {0, 2147483647}, default: 10]
presolve_substitution_maxfillin = 10

# Matrix factorization pivot threshold
# [type: double, advanced: true, range: [0.0008, 0.5], default: 0.1]
factor_pivot_threshold = 0.1

# Matrix factorization pivot tolerance
# [type: double, advanced: true, range: [0, 1], default: 1e-10]
factor_pivot_tolerance = 1e-10

# Tolerance to be satisfied before IPM crossover will start
# [type: double, advanced: true, range: [1e-12, inf], default: 1e-08]
start_crossover_tolerance = 1e-08

# Use original HFactor logic for sparse vs hyper-sparse TRANs
# [type: bool, advanced: true, range: {false, true}, default: true]
use_original_HFactor_logic = true

# Check whether LP is candidate for LiDSE
# [type: bool, advanced: true, range: {false, true}, default: true]
less_infeasible_DSE_check = true

# Use LiDSE if LP has right properties
# [type: bool, advanced: true, range: {false, true}, default: true]
less_infeasible_DSE_choose_row = true
```
