# HiGHS.jl

[![Build Status](https://github.com/jump-dev/HiGHS.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/HiGHS.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/HiGHS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/HiGHS.jl)

[HiGHS.jl](https://github.com/jump-dev/HiGHS.jl) is a wrapper for the
[HiGHS](https://highs.dev) solver.

It has two components:

 - a thin wrapper around the complete C API
 - an interface to [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)

## Affiliation

This wrapper is maintained by the JuMP community and is not an official project
of the HiGHS developers.

## Getting help

If you need help, please ask a question on the [JuMP community forum](https://jump.dev/forum).

If you have a reproducible example of a bug, please [open a GitHub issue](https://github.com/jump-dev/HiGHS.jl/issues/new).

## License

`HiGHS.jl` is licensed under the [MIT License](https://github.com/jump-dev/HiGHS.jl/blob/master/LICENSE.md).

The underlying solver, [ERGO-Code/HiGHS](https://github.com/ERGO-Code/HiGHS), is
licensed under the [MIT license](https://github.com/ERGO-Code/HiGHS/blob/master/LICENSE).

## Installation

Install HiGHS as follows:
```julia
import Pkg
Pkg.add("HiGHS")
```

In addition to installing the HiGHS.jl package, this will also download and
install the HiGHS binaries. You do not need to install HiGHS separately.

To use a custom binary, read the [Custom solver binaries](https://jump.dev/JuMP.jl/stable/developers/custom_solver_binaries/)
section of the JuMP documentation.

## Use with JuMP

To use HiGHS with JuMP, use `HiGHS.Optimizer`:

```julia
using JuMP, HiGHS
model = Model(HiGHS.Optimizer)
# Set options as needed, for example:
set_attribute(model, "presolve", "on")
set_attribute(model, "time_limit", 60.0)
```

## MathOptInterface API

The HiGHS optimizer supports the following constraints and attributes.

List of supported objective functions:

 * [`MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}`](@ref)
 * [`MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}`](@ref)

List of supported variable types:

 * [`MOI.Reals`](@ref)

List of supported constraint types:

 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.EqualTo{Float64}`](@ref)
 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.GreaterThan{Float64}`](@ref)
 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.Interval{Float64}`](@ref)
 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.LessThan{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.EqualTo{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.GreaterThan{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.Integer`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.Interval{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.LessThan{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.Semicontinuous{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.Semiinteger{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.ZeroOne`](@ref)

List of supported model attributes:

 * [`MOI.Name()`](@ref)
 * [`MOI.ObjectiveSense()`](@ref)

## Options

See the [HiGHS documentation](https://ergo-code.github.io/HiGHS/dev/options/definitions/)
for a full list of the available options.

## C API

The C API can be accessed via `HiGHS.Highs_xxx` functions, where the names and
arguments are identical to the C API.

## Threads

HiGHS uses a global scheduler that is shared between threads.

Before changing the number of threads using `MOI.Threads()`, you must call
`Highs_resetGlobalScheduler(1)`:
```julia
using JuMP, HiGHS
model = Model(HiGHS.Optimizer)
Highs_resetGlobalScheduler(1)
set_attribute(model, MOI.NumberOfThreads(), 1)
```

If modifying the number of HiGHS threads across different Julia threads, be sure
to read the docstring of `Highs_resetGlobalScheduler`. In particular, resetting
the scheduler is not thread-safe.
