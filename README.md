_While HiGHS is ready to use today, you should be aware that it is still under
active development. That means it might be slow or unstable on some problems.
However, by reporting these problems you can help make HiGHS better!_

_To report a problem (e.g., incorrect results, or a crash of the solver),
or make a suggestion for how to improve HiGHS, please [file a GitHub issue](https://github.com/jump-dev/HiGHS.jl)._

_To make the problem reproducible, use `JuMP.write_to_file(model, "filename.mps")`
to write your model an MPS file, then upload the MPS file to [https://gist.github.com](https://gist.github.com)
and provide a link to the gist in the GitHub issue._

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

**Minimum version requirement:** HiGHS.jl requres at least Julia v1.3.

Install HiGHS as follows:
```julia
import Pkg
Pkg.add("HiGHS")
```

In addition to installing the HiGHS.jl package, this will also download and
install the HiGHS binaries. (You do not need to install HiGHS separately.)

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
