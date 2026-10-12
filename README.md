[![codecov](https://codecov.io/github/ProjectTorreyPines/TGLF.jl/graph/badge.svg?token=1NY1YYLBWH)](https://codecov.io/github/ProjectTorreyPines/TGLF.jl)
![Docs](https://github.com/ProjectTorreyPines/TGLF.jl/actions/workflows/make_docs.yml/badge.svg)

# TGLF.jl

**TGLF.jl** is the Julia implementation of
[TGLF](https://gafusion.github.io/doc/tglf.html), the quasi-linear model of
gyrokinetic turbulent transport. Given a local plasma equilibrium and its
gradients, TGLF.jl predicts the turbulent particle, energy, and momentum fluxes
that drive cross-field transport in tokamaks.

It is a faithful translation of the Fortran TGLF code — **verified against it**
through a golden-output regression suite — and adds two things the Fortran does
not have: an optional **GPU eigensolver** and **full automatic differentiation**,
so you can take exact gradients of the predicted fluxes with respect to any input.

For the full API reference, see the
[online documentation](https://projecttorreypines.github.io/TGLF.jl/dev).

## Capabilities

1. **Quasi-linear turbulent fluxes** for electrons and an arbitrary number of ion
   species — particle flux, energy flux, toroidal and parallel momentum stress,
   and exchange, returned per `[field, species, type]`.
2. **All TGLF saturation rules** — `SAT_RULE` 0, 1, 2, and 3.
3. **Electrostatic and electromagnetic** runs (`USE_BPER`, `USE_BPAR`).
4. **Verified against Fortran TGLF** — a regression suite checks TGLF.jl's fluxes
   against archived TGLF `out.tglf.gbflux` golden outputs across several cases.
5. **Standard TGLF inputs** — read a plain `input.tglf` file with `readInput`, or
   build/convert from an `InputTGLF`/`InputTJLF` struct.
6. **Multithreaded** over the `ky` spectrum and over a `Vector{InputTJLF}`, with
   automatic, scoped BLAS-thread management so it scales without manual setup.
7. **Optional GPU eigensolver** — pass `use_gpu=true` to offload the dense
   eigenproblem to CUDA (CUSOLVER).
8. **Differentiable** — exact flux gradients via `ForwardDiff` (works on CPU and
   GPU), validated against finite differences.
9. **Fast re-runs** — reuse a previously found width / eigenvalue spectrum via the
   `FIND_WIDTH` and `FIND_EIGEN` flags to skip the expensive search.
10. **Linear wavefunctions** via `get_wavefunction` for mode-structure analysis.

## Installation

TGLF.jl requires Julia **1.11+** and is registered in the
[FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry).

```julia
using Pkg
Pkg.add("TJLF")
```

The GPU and autodiff features are optional
[package extensions](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)):
they activate automatically when you also load `CUDA` (requires a CUDA
toolkit/runtime >= 12.6 for `cusolverDnXgeev`) or `ForwardDiff`. Neither is needed
for ordinary CPU runs.

## Quick start

> Note: `run` and `run_tjlf` are not exported (`run` would clash with
> `Base.run`), so call them as `TJLF.run` / `TJLF.run_tjlf`. The flux helpers
> (`Qe`, `Qi`, …) are likewise accessed as `TJLF.Qe`, etc.

```julia
using TJLF

# 1. Build the input from a standard input.tglf file (NOT input.tglf.gen)
input = readInput("input.tglf")

# 2. Run and get the summed fluxes you usually care about: [field, species, type]
fluxes = TJLF.run_tjlf(input)

# 3. Convenience accessors for the common channels (each takes that array)
Qe = TJLF.Qe(fluxes)   # electron energy flux
Qi = TJLF.Qi(fluxes)   # total ion energy flux
Γe = TJLF.Γe(fluxes)   # electron particle flux
Πi = TJLF.Πi(fluxes)   # ion toroidal momentum stress
```

Need everything (weights, spectra, eigenvalues, mode structure)? Call `TJLF.run`,
which returns a NamedTuple:

```julia
out = TJLF.run(input)
# out.QL_weights       quasi-linear weights      [field, species, mode, ky, type]
# out.eigenvalue       growth rate & frequency   [mode, ky, type]
# out.QL_flux_out      summed fluxes             [field, species, type]   (== run_tjlf output)
# out.flux_spectrum    flux per ky               [field, species, mode, ky, type]
# out.field_weight_out mode (eigen) structure    [field, basis, mode, ky]
```

Run many cases in parallel — just pass a vector (threaded automatically):

```julia
# Vector of [field, species, type] arrays, one per input
fluxes = TJLF.run_tjlf([input1, input2, input3])
```

### GPU

```julia
using TJLF, CUDA
fluxes = TJLF.run_tjlf(input; use_gpu=true)   # dense eigensolves run on the GPU
```

### Gradients (autodiff)

Loading `ForwardDiff` makes the flux calculation differentiable. Because the AD
type flows through the whole `InputTJLF`, rebuild the input at the differentiated
element type `T` rather than mutating a `Float64` struct:

```julia
using TJLF, ForwardDiff

# rebuild a Float64 InputTJLF as an InputTJLF{T}
function as_eltype(::Type{T}, base::TJLF.InputTJLF{Float64}) where {T<:Real}
    inp = TJLF.InputTJLF{T}(base.NS, length(base.KY_SPECTRUM))
    for fn in fieldnames(TJLF.InputTJLF)
        v = getfield(base, fn)
        v isa Float64         ? setfield!(inp, fn, T(v))          :
        v isa Vector{Float64} ? setfield!(inp, fn, T.(v))         :
        v isa Vector{ComplexF64} ? setfield!(inp, fn, Complex{T}.(v)) :
                                  setfield!(inp, fn, v)
    end
    return inp
end

# d(electron energy flux) / d(ion temperature gradient)
function Qe_of_RLTS2(x::T) where {T<:Real}
    inp = as_eltype(T, input)
    inp.RLTS[2] = x
    return TJLF.Qe(TJLF.run_tjlf(inp))
end
dQe = ForwardDiff.derivative(Qe_of_RLTS2, input.RLTS[2])
```

`ForwardDiff.gradient` works the same way for a vector of inputs. See
[`test/runtests_ad.jl`](test/runtests_ad.jl) for the full helper and more
examples (gradients, and AD across all saturation rules).

## Outputs and array indices

Several outputs are multi-dimensional; the index conventions are:

| Output | Shape | Indices |
|--------|-------|---------|
| `QL_weights` | 5D | `[field, species, mode, ky, type]` |
| `flux_spectrum` | 5D | `[field, species, mode, ky, type]` |
| `eigenvalue` | 3D | `[mode, ky, type]` |
| `QL_flux_out` (= `run_tjlf` output) | 3D | `[field, species, type]` |
| `field_weight_out` | 4D | `[field, basis, mode, ky]` |

Where:

- **type** (fluxes): `1`=particle, `2`=energy, `3`=toroidal stress, `4`=parallel stress, `5`=exchange
- **type** (eigenvalue): `1`=gamma (growth rate), `2`=frequency
- **species**: `1`=electron, `2+`=ions
- **mode**: `1`=most unstable
- **field**: `1`=electrostatic, `2`=perpendicular magnetic (`USE_BPER`), `3`=parallel magnetic (`USE_BPAR`)

To reduce `run_tjlf`'s `[field, species, type]` output to per-species totals,
sum over the field dimension, e.g. `sum(fluxes; dims=1)[1, :, :]` gives
`[species, type]` (this is exactly what the flux helpers do internally).

Index order is chosen to take advantage of Julia's column-major memory layout.

---

# Advanced topics

The sections below cover internals and tuning. New users can skip them.

## The InputTJLF struct

`InputTJLF` mostly mirrors the `input.tglf` variables. The translation from an
`InputTGLF` struct is handled in `tjlf_modules.jl`. By design, `InputTJLF` fields
have **no meaningful defaults** — `checkInput` will error on an incompletely
populated struct, so you are always aware of exactly which parameters you are
using.

### TGLF.jl-specific parameters

- `FIND_WIDTH::Bool` — repurposed from the TGLF flag: tells TGLF.jl whether to solve
  for the Gaussian mode width (`true`) or use an externally provided width
  spectrum (`false`), avoiding extra width-search calls.
- `WIDTH_SPECTRUM::Vector` — the per-ky width spectrum. Populated automatically
  when `FIND_WIDTH=true`; consumed as input when `FIND_WIDTH=false`. A common
  workflow: run once with `FIND_WIDTH=true`, tweak gradients on the same struct,
  then set `FIND_WIDTH=false` to reuse the saved widths.
- `FIND_EIGEN::Bool` — selects the eigenvalue solver (see below). `true` (default)
  is the robust choice. Setting it to `false` is only sensible when `FIND_WIDTH`
  is also `false`.
- `EIGEN_SPECTRUM::Vector` — external eigenvalue spectrum used as the initial
  guess (the `sigma` shift) for the iterative solver when `FIND_EIGEN=false`.
- `SMALL::Float` — small value used to rewrite the eigenvalue problem as a linear
  system, mirroring the Fortran (default `1e-13`).

### Deleted parameters

Unused legacy parameters were removed from `InputTJLF`:
`USE_TRANSPORT_MODEL`, `GEOMETRY_FLAG`, `B_MODEL_SA`, `FT_MODEL_SA`,
`VPAR_SHEAR_MODEL`, `WRITE_WAVEFUNCTION_FLAG`, `VTS_SHEAR`, `VNS_SHEAR`, `VEXB`,
`RMIN_SA`, `RMAJ_SA`, `Q_SA`, `SHAT_SA`, `ALPHA_SA`, `XWELL_SA`, `THETA0_SA`,
`NN_MAX_ERROR`. See commits
[b9a9c99](https://github.com/ProjectTorreyPines/TGLF.jl/commit/b9a9c99b40a3ba5fcbc683eb94df4ec8b6fd5671)
and
[18f810f](https://github.com/ProjectTorreyPines/TGLF.jl/commit/18f810fdf277a691b584f987781b976cccdcbb5b).

## Eigenvalue solvers

TGLF.jl supports multiple eigenvalue-solver strategies with different
robustness/performance trade-offs:

- **Direct LAPACK (default, `FIND_EIGEN=true`)** — `gesv!` + `geev!` from
  `LinearAlgebra.LAPACK`. Converts the generalized eigenproblem to a linear system
  and finds all eigenvalues. Robust, fast, and scales well with threads.
- **Iterative KrylovKit (`FIND_EIGEN=false`)** — `KrylovKit.eigsolve` with
  shift-and-invert, using `EIGEN_SPECTRUM` as the initial guess. Fast when a good
  guess is available, but less robust; falls back to LAPACK with a warning if it
  fails.

Solver selection logic:

1. `FIND_EIGEN=true` → robust LAPACK solver.
2. `FIND_EIGEN=false` with a valid `EIGEN_SPECTRUM` → KrylovKit iterative solver.
3. KrylovKit failure → automatic fallback to LAPACK.

A legacy Arpack `eigs()` path has been discontinued in favor of KrylovKit.

## GPU

Pass `use_gpu=true` to `run`/`run_tjlf` to move the dense eigensolves onto the
GPU via CUSOLVER (`getrf`/`getrs` for the solve, `Xgeev` for the diagonalization).
`pick_device(:cpu | :gpu | :auto)` resolves the device, erroring if `:gpu` is
requested without a functional CUDA runtime. Requires a CUDA toolkit/runtime
>= 12.6.

## Automatic differentiation

Loading `ForwardDiff` activates `TJLFForwardDiffExt`, which makes the whole flux
calculation differentiable. LAPACK has no `Dual`-number support, so for
`ForwardDiff.Dual` element types TGLF.jl computes the `Float64` eigensystem with
LAPACK and propagates the perturbation analytically via the implicit function
theorem. These specializations are defined on **TGLF.jl-owned** eigen wrappers, not
on `LinearAlgebra.eigen`, so there is no type piracy and no effect on other
packages. With `use_gpu=true` the same `Dual` eigensolve + gradient runs on the
GPU. Gradients match finite differences (and the GPU path matches the CPU path to
~1e-10).

## Multithreading

See `test/test_multithreads` for benchmarks.

TGLF.jl parallelizes over `ky` with Julia threads, each doing small dense LAPACK
solves. If BLAS is also multithreaded, those kernels oversubscribe cores (e.g.
8 Julia threads × 128 BLAS threads) and per-solve time can balloon ~100×.

**TGLF.jl handles this automatically.** It wraps its threaded entry points (`run`,
`run_tjlf`, and the transport-model loop) with `with_blas_threads(1)`, which sets
the BLAS thread count to 1 for the duration of the call and restores it
afterwards. This is scoped to TGLF.jl's own drivers, so loading TGLF.jl does not change
BLAS threading for the rest of a host application (e.g. FUSE) — **no manual setup
is required**.

If you thread over your *own* TGLF.jl calls, set BLAS=1 once outside the threaded
region (nested scopes short-circuit, so there are no BLAS races):

```julia
TJLF.with_blas_threads(1) do
    Threads.@threads for i in eachindex(cases)
        TJLF.run(cases[i])
    end
end
```

Setting `export OMP_NUM_THREADS=1` in the shell is a reasonable
belt-and-suspenders but is no longer required for TGLF.jl's own scaling.

## Verification

TGLF.jl's trust anchor is its test suite. `test/runtests_regressions.jl` runs each
case in `test/tglf_regression/` through TGLF.jl and compares the fluxes against the
archived Fortran TGLF `out.tglf.gbflux` golden output (to `atol=1e-2`). Companion
suites cover saturation rules (`runtests_sat.jl`), electromagnetic runs
(`runtests_EM.jl`), the transport model (`runtests_tm.jl`), the eigensolver
(`runtests_eigen.jl`), the ky grid (`runtests_kygrid.jl`), and autodiff
(`runtests_ad.jl`).

## Citation

If this software contributes to an academic publication, please refer to it as **TGLF.jl** and cite it as follows:

> T.F. Neiser, D. Sun, B. Agnew, T. Slendebroek, O. Meneghini, B.C. Lyons, A. Ghiozzi, J. McClenaghan, G. Staebler and J. Candy, _TJLF: The quasi-linear model of gyrokinetic transport TGLF translated to Julia_, APS Meeting Abstracts (2024)
