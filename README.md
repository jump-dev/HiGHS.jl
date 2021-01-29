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

The package is not registered in the [General registry](https://github.com/JuliaRegistries/General/)
and so must be installed as follows:

```julia
import Pkg
Pkg.add(Pkg.PackageSpec(url="https://github.com/jump-dev/HiGHS.jl"))
```

In addition to installing the HiGHS.jl package, this will also download and
install the HiGHS binaries. (You do not need to install HiGHS separately.)
