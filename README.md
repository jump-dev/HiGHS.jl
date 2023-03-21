# HiGHS.jl

[![Build Status](https://github.com/jump-dev/HiGHS.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/HiGHS.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/HiGHS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/HiGHS.jl)

HiGHS.jl is a wrapper for the [HiGHS](https://highs.dev) solver.

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

See the [HiGHS documentation](https://ergo-code.github.io/HiGHS/dev/options/definitions/)
for a full list of the available options.
