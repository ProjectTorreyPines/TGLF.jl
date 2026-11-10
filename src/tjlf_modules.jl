abstract type AbstractAve{T<:Number} end

Base.@kwdef mutable struct InputTGLF{T<:Real}
    # Concretely typed (no Union{...,Missing}) for type stability — same recipe and
    # rationale as InputTJLF below. Float defaults are NaN "unset" sentinels
    # (detected with `is_unset`), written as plain literals rather than T(NaN) so the
    # @kwdef constructor can infer T from them: `InputTGLF()` and `InputTGLF(NS=3)`
    # give an InputTGLF{Float64}, `InputTGLF{Float32}()` converts. Ints/bools/strings
    # carry the Fortran TGLF defaults,
    # so an omitted switch means "TGLF default", never a silent zero (an unset
    # SIGN_IT=0 used to silently kill VEXB_SHEAR).
    SIGN_BT::Int = 1
    SIGN_IT::Int = 1
    NS::Int = 2
    ZMAJ_LOC::T = NaN
    DRMINDX_LOC::T = NaN
    DZMAJDX_LOC::T = NaN
    S_DELTA_LOC::T = NaN
    ZETA_LOC::T = NaN
    S_ZETA_LOC::T = NaN

    MASS_1::T = NaN
    ZS_1::T = NaN
    AS_1::T = NaN
    TAUS_1::T = NaN

    MASS_2::T = NaN
    ZS_2::T = NaN
    VPAR_2::T = NaN
    VPAR_SHEAR_2::T = NaN

    MASS_3::T = NaN
    ZS_3::T = NaN
    RLTS_3::T = NaN
    TAUS_3::T = NaN
    VPAR_3::T = NaN
    VPAR_SHEAR_3::T = NaN

    # TGLF-NN uses 3 species
    # This is why parameters for species 1:3 are sorted differently than 4:10
    MASS_4::T = NaN
    AS_4::T = NaN
    ZS_4::T = NaN
    RLNS_4::T = NaN
    RLTS_4::T = NaN
    TAUS_4::T = NaN
    VPAR_4::T = NaN
    VPAR_SHEAR_4::T = NaN

    MASS_5::T = NaN
    AS_5::T = NaN
    ZS_5::T = NaN
    RLNS_5::T = NaN
    RLTS_5::T = NaN
    TAUS_5::T = NaN
    VPAR_5::T = NaN
    VPAR_SHEAR_5::T = NaN

    MASS_6::T = NaN
    AS_6::T = NaN
    ZS_6::T = NaN
    RLNS_6::T = NaN
    RLTS_6::T = NaN
    TAUS_6::T = NaN
    VPAR_6::T = NaN
    VPAR_SHEAR_6::T = NaN

    MASS_7::T = NaN
    AS_7::T = NaN
    ZS_7::T = NaN
    RLNS_7::T = NaN
    RLTS_7::T = NaN
    TAUS_7::T = NaN
    VPAR_7::T = NaN
    VPAR_SHEAR_7::T = NaN

    MASS_8::T = NaN
    AS_8::T = NaN
    ZS_8::T = NaN
    RLNS_8::T = NaN
    RLTS_8::T = NaN
    TAUS_8::T = NaN
    VPAR_8::T = NaN
    VPAR_SHEAR_8::T = NaN

    MASS_9::T = NaN
    AS_9::T = NaN
    ZS_9::T = NaN
    RLNS_9::T = NaN
    RLTS_9::T = NaN
    TAUS_9::T = NaN
    VPAR_9::T = NaN
    VPAR_SHEAR_9::T = NaN

    MASS_10::T = NaN
    AS_10::T = NaN
    ZS_10::T = NaN
    RLNS_10::T = NaN
    RLTS_10::T = NaN
    TAUS_10::T = NaN
    VPAR_10::T = NaN
    VPAR_SHEAR_10::T = NaN

    AS_2::T = NaN
    AS_3::T = NaN
    BETAE::T = NaN
    DEBYE::T = NaN
    DELTA_LOC::T = NaN
    DRMAJDX_LOC::T = NaN
    KAPPA_LOC::T = NaN
    P_PRIME_LOC::T = NaN
    Q_LOC::T = NaN
    Q_PRIME_LOC::T = NaN
    RLNS_1::T = NaN
    RLNS_2::T = NaN
    RLNS_3::T = NaN
    RLTS_1::T = NaN
    RLTS_2::T = NaN
    RMAJ_LOC::T = NaN
    RMIN_LOC::T = NaN
    S_KAPPA_LOC::T = NaN
    TAUS_2::T = NaN
    VEXB_SHEAR::T = NaN
    VPAR_1::T = NaN
    VPAR_SHEAR_1::T = NaN
    XNUE::T = NaN
    ZEFF::T = NaN

    # switches (bool/int/string defaults = Fortran tglf_interface.f90 defaults, see NOTE above)
    UNITS::String = "GYRO"
    ALPHA_ZF::T = NaN
    USE_MHD_RULE::Bool = true
    NKY::Int = 12
    SAT_RULE::Int = 0
    KYGRID_MODEL::Int = 1
    NMODES::Int = 2
    NBASIS_MIN::Int = 2
    NBASIS_MAX::Int = 4
    XNU_MODEL::Int = 2
    USE_AVE_ION_GRID::Bool = false
    ALPHA_QUENCH::Int = 0
    ALPHA_MACH::T = NaN
    WDIA_TRAPPED::T = NaN
    USE_BPAR::Bool = false
    USE_BPER::Bool = false
    USE_PRESETS::Bool = true

    _Qgb::T = NaN

    # missing
    USE_BISECTION::Bool = true
    USE_INBOARD_DETRAPPED::Bool = false
    NEW_EIKONAL::Bool = true
    FIND_WIDTH::Bool = true
    IFLUX::Bool = true
    ADIABATIC_ELEC::Bool = false

    NWIDTH::Int = 21
    NXGRID::Int = 16
    VPAR_MODEL::Int = 0
    VPAR_SHEAR_MODEL::Int = 1
    IBRANCH::Int = -1

    KY::T = 0.3
    ALPHA_E::T = 1.0
    ALPHA_P::T = 1.0
    XNU_FACTOR::T = 1.0
    DEBYE_FACTOR::T = 1.0
    RLNP_CUTOFF::T = 18.0
    WIDTH::T = 1.65
    WIDTH_MIN::T = 0.3
    BETA_LOC::T = 0.0
    KX0_LOC::T = 0.0
    PARK::T = 1.0
    GHAT::T = 1.0
    GCHAT::T = 1.0
    WD_ZERO::T = 0.1
    LINSKER_FACTOR::T = 0.0
    GRADB_FACTOR::T = 0.0
    FILTER::T = 2.0
    THETA_TRAPPED::T = 0.7
    ETG_FACTOR::T = 1.25
    DAMP_PSI::T = 0.0
    DAMP_SIG::T = 0.0

    # MXH params. NaN sentinel (not zero) so "unset by geometry" stays detectable:
    # conversions (update_input_tjlf!, tglf_to_cgyro) skip unset fields and let the
    # target's own 0.0 defaults stand.
    SHAPE_COS0::T = NaN
    SHAPE_COS1::T = NaN
    SHAPE_COS2::T = NaN
    SHAPE_COS3::T = NaN
    SHAPE_COS4::T = NaN
    SHAPE_COS5::T = NaN
    SHAPE_COS6::T = NaN

    SHAPE_SIN3::T = NaN
    SHAPE_SIN4::T = NaN
    SHAPE_SIN5::T = NaN
    SHAPE_SIN6::T = NaN

    SHAPE_S_COS0::T = NaN
    SHAPE_S_COS1::T = NaN
    SHAPE_S_COS2::T = NaN
    SHAPE_S_COS3::T = NaN
    SHAPE_S_COS4::T = NaN
    SHAPE_S_COS5::T = NaN
    SHAPE_S_COS6::T = NaN

    SHAPE_S_SIN3::T = NaN
    SHAPE_S_SIN4::T = NaN
    SHAPE_S_SIN5::T = NaN
    SHAPE_S_SIN6::T = NaN

    # Bitmask (by field index) of fields assigned via setproperty! since construction:
    # ints/bools have no unset sentinel, so "did the user set this?" (FUSE's
    # custom_input_files mask) cannot be inferred from the value. Inline UInt128s
    # rather than a Set{Symbol} so tracking adds zero heap allocations. The underscore
    # prefix keeps them out of save/serialize loops; `setfield!` writes (internal
    # copies/transforms, not user intent) are deliberately not recorded.
    _setmask_lo::UInt128 = UInt128(0)
    _setmask_hi::UInt128 = UInt128(0)
end

@assert fieldcount(InputTGLF) <= 256 "InputTGLF outgrew the 256-bit _setmask; widen the mask fields"

# Fortran namelist parsers often hand string values back with their quotes intact
# (UNITS='CGYRO'). Normalize on assignment so programmatic sets behave like readInput,
# which strips quotes while parsing input.tglf.
@inline _strip_quotes(v::AbstractString) = String(strip(v, ('\'', '"', ' ')))

@inline function Base.setproperty!(input_tglf::InputTGLF, field::Symbol, value)
    if field !== :_setmask_lo && field !== :_setmask_hi
        # fieldindex const-folds when `field` is a literal (input_tglf.NKY = 12);
        # for runtime symbols it is the same C lookup setfield! performs anyway.
        idx = Base.fieldindex(typeof(input_tglf), field)
        if idx <= 128
            setfield!(input_tglf, :_setmask_lo,
                getfield(input_tglf, :_setmask_lo) | (UInt128(1) << (idx - 1)))
        else
            setfield!(input_tglf, :_setmask_hi,
                getfield(input_tglf, :_setmask_hi) | (UInt128(1) << (idx - 129)))
        end
    end
    if value isa AbstractString
        value = _strip_quotes(value)
    end
    return setfield!(input_tglf, field, convert(fieldtype(typeof(input_tglf), field), value))
end

"""
    was_set(input_tglf::InputTGLF, field::Symbol) -> Bool

`true` when `field` was assigned via `input_tglf.field = value` after construction.
Use this instead of comparing against a default when you need to distinguish "the user
asked for this value" from "this is just the struct default" — the two are
indistinguishable by value for ints, bools and strings. Zero-allocation.
"""
@inline function was_set(input_tglf::InputTGLF, field::Symbol)
    idx = Base.fieldindex(typeof(input_tglf), field, false)
    idx == 0 && return false
    if idx <= 128
        return getfield(input_tglf, :_setmask_lo) & (UInt128(1) << (idx - 1)) != 0
    else
        return getfield(input_tglf, :_setmask_hi) & (UInt128(1) << (idx - 129)) != 0
    end
end

"""
    set_fields(input_tglf::InputTGLF) -> Set{Symbol}

The set of fields assigned since construction (see [`was_set`](@ref)). Empty for a
struct populated exclusively through `setfield!`. Materializes a fresh Set on each
call — use [`was_set`](@ref) in loops that must not allocate.
"""
function set_fields(input_tglf::InputTGLF)
    T = typeof(input_tglf)
    return Set{Symbol}(f for f in fieldnames(T) if was_set(input_tglf, f))
end

"""
    is_unset(v) -> Bool

`true` when a sentinel-typed input field holds its "never populated" sentinel:
`missing` (legacy), a NaN float (including `ForwardDiff.Dual`), or an empty
string. Ints and bools have no detectable unset state — they carry the Fortran
TGLF interface defaults instead of sentinels, so an "unset" int/bool is simply
the TGLF default value. This is the replacement for `ismissing` checks on
`InputTGLF`/`InputTJLF` fields now that both structs are concretely typed.
"""
is_unset(::Missing) = true
is_unset(v::Real) = v isa Integer ? false : isnan(v)
is_unset(v::AbstractString) = isempty(v)
is_unset(::Any) = false

"""
    InputTJLF{T<:Real}

The complete set of inputs for a TJLF run — the Julia counterpart of a TGLF
`input.tglf` file. Fields largely mirror the TGLF input variables (physics
switches like `USE_BPER`/`USE_BPAR`/`USE_MHD_RULE`, the saturation rule
`SAT_RULE`, species arrays `ZS`/`MASS`/`AS`/`TAUS`/`RLNS`/`RLTS`/`VPAR`, Miller
geometry, and the `ky` grid controls). Build one with `readInput("input.tglf")`,
convert from a `TurbulentTransport.InputTGLF`, or construct directly via
`InputTJLF{T}(ns, nky)` and populate the fields.

Fields are **concretely typed** (no `Union{Missing,...}`) for type stability in
the hot eigensolver path. Float defaults are `NaN` sentinels that construction
is expected to overwrite — `checkInput` errors on any leftover `NaN`, so floats
deliberately have no silent physics defaults. Bools and ints have no detectable
unset state, so they carry the Fortran TGLF defaults
(gacode `tglf/src/tglf_interface.f90`) instead of `0`/`false` sentinels: an
omitted switch means "TGLF default", never a silent zero.

A few fields are TJLF-specific (repurposed or added relative to TGLF):
- `FIND_WIDTH::Bool` — solve for the Gaussian mode width (`true`) or reuse an
  externally supplied `WIDTH_SPECTRUM` (`false`).
- `WIDTH_SPECTRUM::Vector` — per-`ky` widths; filled when `FIND_WIDTH=true`,
  consumed when `false`.
- `FIND_EIGEN::Bool` — select the eigenvalue solver (robust LAPACK when `true`).
- `EIGEN_SPECTRUM::Vector` — initial eigenvalue guess used when `FIND_EIGEN=false`.
- `SMALL::Float` — small value used to rewrite the eigenproblem as a linear
  system (default `1e-13`).

See also [`readInput`](@ref), [`run`](@ref), and `run_tjlf`.
"""
Base.@kwdef mutable struct InputTJLF{T<:Real}
    # NOTE: fields are concretely typed (no Union{...,Missing}): every inputs.FIELD
    # read in the hot eigensolver path must infer concretely, or inference poisons
    # all of TJLF and heap-boxes millions of scalars. Float defaults are NaN "unset"
    # sentinels enforced by checkInput; bools/ints/strings carry the Fortran TGLF
    # defaults (tglf_interface.f90), so an omitted switch means "TGLF default",
    # never a silent zero. UNITS: GYRO default; apply_presets! forces CGYRO for SAT2/3.
    UNITS::String = "GYRO"

    USE_BPER::Bool = false
    USE_BPAR::Bool = false
    USE_MHD_RULE::Bool = true
    USE_BISECTION::Bool = true
    USE_INBOARD_DETRAPPED::Bool = false
    USE_AVE_ION_GRID::Bool = false
    USE_PRESETS::Bool = true  # gates apply_presets! (hard-coded .TRUE. in Fortran tglf_startup.f90)
    NEW_EIKONAL::Bool = true
    FIND_WIDTH::Bool = true
    IFLUX::Bool = true
    ADIABATIC_ELEC::Bool = false

    SAT_RULE::Int = 0
    NS::Int = 2
    NMODES::Int = 2
    NWIDTH::Int = 21
    NBASIS_MAX::Int = 4
    NBASIS_MIN::Int = 2
    NXGRID::Int = 16
    NKY::Int = 12
    KYGRID_MODEL::Int = 1
    XNU_MODEL::Int = 2
    VPAR_MODEL::Int = 0
    IBRANCH::Int = -1

    ZS::Vector{T} = T[]
    MASS::Vector{T} = T[]
    RLNS::Vector{T} = T[]
    RLTS::Vector{T} = T[]
    TAUS::Vector{T} = T[]
    AS::Vector{T} = T[]
    VPAR::Vector{T} = T[]
    VPAR_SHEAR::Vector{T} = T[]

    # NOT IN TGLF (or missing from InputTGLF structure)
    WIDTH_SPECTRUM::Vector{T} = T[]    # TJLF-specific
    KY_SPECTRUM::Vector{T} = T[]       # TJLF-specific  
    EIGEN_SPECTRUM::Vector{Complex{T}} = Complex{T}[]  # TJLF-specific
    FIND_EIGEN::Bool = true             # TGLF parameter but missing from InputTGLF
    # NOT IN TGLF (or missing from InputTGLF structure)

    SIGN_BT::Int = 1
    SIGN_IT::Int = 1
    KY::T = T(NaN)

    VEXB_SHEAR::T = T(NaN)
    BETAE::T = T(NaN)
    XNUE::T = T(NaN)
    ZEFF::T = T(NaN)
    DEBYE::T = T(NaN)

    ALPHA_MACH::T = T(NaN)
    ALPHA_E::T = T(NaN)
    ALPHA_P::T = T(NaN)
    ALPHA_QUENCH::Int = 0
    ALPHA_ZF::T = T(NaN)
    XNU_FACTOR::T = T(NaN)
    DEBYE_FACTOR::T = T(NaN)
    ETG_FACTOR::T = T(NaN)
    RLNP_CUTOFF::T = T(NaN)

    WIDTH::T = T(NaN)
    WIDTH_MIN::T = T(NaN)

    RMIN_LOC::T = T(NaN)
    RMAJ_LOC::T = T(NaN)
    ZMAJ_LOC::T = T(NaN)
    DRMINDX_LOC::T = T(NaN)
    DRMAJDX_LOC::T = T(NaN)
    DZMAJDX_LOC::T = T(NaN)
    Q_LOC::T = T(NaN)
    KAPPA_LOC::T = T(NaN)
    S_KAPPA_LOC::T = T(NaN)
    DELTA_LOC::T = T(NaN)
    S_DELTA_LOC::T = T(NaN)
    ZETA_LOC::T = T(NaN)
    S_ZETA_LOC::T = T(NaN)
    P_PRIME_LOC::T = T(NaN)
    Q_PRIME_LOC::T = T(NaN)
    BETA_LOC::T = T(NaN)
    KX0_LOC::T = T(NaN)

    DAMP_PSI::T = T(NaN)
    DAMP_SIG::T = T(NaN)
    WDIA_TRAPPED::T = T(NaN)
    PARK::T = T(NaN)
    GHAT::T = T(NaN)
    GCHAT::T = T(NaN)
    WD_ZERO::T = T(NaN)
    LINSKER_FACTOR::T = T(NaN)
    GRADB_FACTOR::T = T(NaN)
    FILTER::T = T(NaN)
    THETA_TRAPPED::T = T(NaN)
    SMALL::T = T(NaN)

    #MXH params (optional; default to 0 so checkInput passes when unset by geometry)
    SHAPE_COS0::T = zero(T)
    SHAPE_COS1::T = zero(T)
    SHAPE_COS2::T = zero(T)
    SHAPE_COS3::T = zero(T)
    SHAPE_COS4::T = zero(T)
    SHAPE_COS5::T = zero(T)
    SHAPE_COS6::T = zero(T)

    SHAPE_SIN3::T = zero(T)
    SHAPE_SIN4::T = zero(T)
    SHAPE_SIN5::T = zero(T)
    SHAPE_SIN6::T = zero(T)

    SHAPE_S_COS0::T = zero(T)
    SHAPE_S_COS1::T = zero(T)
    SHAPE_S_COS2::T = zero(T)
    SHAPE_S_COS3::T = zero(T)
    SHAPE_S_COS4::T = zero(T)
    SHAPE_S_COS5::T = zero(T)
    SHAPE_S_COS6::T = zero(T)

    SHAPE_S_SIN3::T = zero(T)
    SHAPE_S_SIN4::T = zero(T)
    SHAPE_S_SIN5::T = zero(T)
    SHAPE_S_SIN6::T = zero(T)

    USE_TRANSPORT_MODEL::Bool = true
end

# String-valued assignments only (UNITS); numeric hot-path assignments keep the
# default setproperty!. See _strip_quotes above for the namelist-quote rationale.
@inline Base.setproperty!(input_tjlf::InputTJLF, field::Symbol, value::AbstractString) =
    setfield!(input_tjlf, field, convert(fieldtype(typeof(input_tjlf), field), _strip_quotes(value)))

function InputTJLF{T}(ns::Int, nky::Int) where {T<:Real}
    InputTJLF{T}(;
        UNITS = "GYRO",
        NS = ns,
        ZS = fill(T(NaN), ns),
        MASS = fill(T(NaN), ns),
        RLNS = fill(T(NaN), ns),
        RLTS = fill(T(NaN), ns),
        TAUS = fill(T(NaN), ns),
        AS = fill(T(NaN), ns),
        VPAR = fill(T(NaN), ns),
        VPAR_SHEAR = fill(T(NaN), ns),
        WIDTH_SPECTRUM = fill(T(NaN), nky),
        KY_SPECTRUM = fill(T(NaN), nky),
        EIGEN_SPECTRUM = fill(Complex{T}(NaN*im), nky),
        SIGN_BT = 1,
        SIGN_IT = 1,
        KY = T(NaN),
        VEXB_SHEAR = T(NaN),
        BETAE = T(NaN),
        XNUE = T(NaN),
        ZEFF = T(NaN),
        DEBYE = T(NaN),
        ALPHA_MACH = T(NaN),
        ALPHA_QUENCH = 0,
        ALPHA_E = T(NaN),
        ALPHA_P = T(NaN),
        ALPHA_ZF = T(NaN),
        XNU_FACTOR = T(NaN),
        DEBYE_FACTOR = T(NaN),
        ETG_FACTOR = T(NaN),
        RLNP_CUTOFF = T(NaN),
        WIDTH = T(NaN),
        WIDTH_MIN = T(NaN),
        RMIN_LOC = T(NaN),
        RMAJ_LOC = T(NaN),
        ZMAJ_LOC = T(NaN),
        DRMINDX_LOC = T(NaN),
        DRMAJDX_LOC = T(NaN),
        DZMAJDX_LOC = T(NaN),
        Q_LOC = T(NaN),
        KAPPA_LOC = T(NaN),
        S_KAPPA_LOC = T(NaN),
        DELTA_LOC = T(NaN),
        S_DELTA_LOC = T(NaN),
        ZETA_LOC = T(NaN),
        S_ZETA_LOC = T(NaN),
        P_PRIME_LOC = T(NaN),
        Q_PRIME_LOC = T(NaN),
        BETA_LOC = T(NaN),
        KX0_LOC = T(NaN),
        DAMP_PSI = T(NaN),
        DAMP_SIG = T(NaN),
        WDIA_TRAPPED = T(NaN),
        PARK = T(NaN),
        GHAT = T(NaN),
        GCHAT = T(NaN),
        WD_ZERO = T(NaN),
        LINSKER_FACTOR = T(NaN),
        GRADB_FACTOR = T(NaN),
        FILTER = T(NaN),
        THETA_TRAPPED = T(NaN),
        SMALL = T(1.0e-13),
        SHAPE_COS0 = T(0.0),
        SHAPE_COS1 = T(0.0),
        SHAPE_COS2 = T(0.0),
        SHAPE_COS3 = T(0.0),
        SHAPE_COS4 = T(0.0),
        SHAPE_COS5 = T(0.0),
        SHAPE_COS6 = T(0.0),
        SHAPE_SIN3 = T(0.0),
        SHAPE_SIN4 = T(0.0),
        SHAPE_SIN5 = T(0.0),
        SHAPE_SIN6 = T(0.0),
        SHAPE_S_COS0 = T(0.0),
        SHAPE_S_COS1 = T(0.0),
        SHAPE_S_COS2 = T(0.0),
        SHAPE_S_COS3 = T(0.0),
        SHAPE_S_COS4 = T(0.0),
        SHAPE_S_COS5 = T(0.0),
        SHAPE_S_COS6 = T(0.0),
        SHAPE_S_SIN3 = T(0.0),
        SHAPE_S_SIN4 = T(0.0),
        SHAPE_S_SIN5 = T(0.0),
        SHAPE_S_SIN6 = T(0.0),
        USE_TRANSPORT_MODEL = true
    )
end

"""
    InputTJLF(input_tglf::InputTGLF)

Generates an InputTJLF from a InputTGLF
"""
function InputTJLF{T}(input_tglf::InputTGLF{T}) where {T<:Real}
    nky = get_ky_spectrum_size(input_tglf.NKY, input_tglf.KYGRID_MODEL)
    input_tjlf = InputTJLF{T}(input_tglf.NS, nky)
    return update_input_tjlf!(input_tjlf, input_tglf)
end

# Copy every field shared by InputTGLF and InputTJLF, skipping unset (NaN/empty)
# source fields so the target's default stands. Unrolled at compile time to typed
# getfield/setfield! — a runtime-Symbol loop infers `Any` and boxes each scalar.
@generated function _copy_shared_fields!(input_tjlf::InputTJLF{T}, input_tglf::InputTGLF{T}) where {T<:Real}
    shared = intersect(fieldnames(InputTGLF), fieldnames(InputTJLF))
    assigns = map(collect(shared)) do f
        qf = QuoteNode(f)
        :(let v = getfield(input_tglf, $qf)
              is_unset(v) || setfield!(input_tjlf, $qf, v)
          end)
    end
    quote
        $(assigns...)
        return input_tjlf
    end
end

# Copy the per-species scalars (ZS_i, AS_i, ...) into the InputTJLF vectors, unrolled
# so the field symbols are compile-time constants (the old `Symbol("ZS_", i)` form
# allocated a Symbol per species per field and inferred to `Any`).
@generated function _copy_species_fields!(input_tjlf::InputTJLF{T}, input_tglf::InputTGLF{T}) where {T<:Real}
    stmts = Expr[]
    for i in 1:10
        block = map((:ZS, :AS, :MASS, :RLNS, :RLTS, :TAUS, :VPAR, :VPAR_SHEAR)) do name
            :(getfield(input_tjlf, $(QuoteNode(name)))[$i] = getfield(input_tglf, $(QuoteNode(Symbol(name, :_, i)))))
        end
        push!(stmts, :(if ns >= $i
            $(block...)
        end))
    end
    quote
        ns = input_tglf.NS
        $(stmts...)
        return input_tjlf
    end
end

"""
    update_input_tjlf!(input_tjlf::InputTJLF, input_tglf::InputTGLF)

Modifies an InputTJLF from a InputTGLF
"""
function update_input_tjlf!(input_tjlf::InputTJLF{T}, input_tglf::InputTGLF{T}) where {T<:Real}
    input_tjlf.NWIDTH = 21

    _copy_shared_fields!(input_tjlf, input_tglf)
    _copy_species_fields!(input_tjlf, input_tglf)

    input_tjlf.WIDTH_SPECTRUM .= input_tjlf.WIDTH

    # Set defaults ONLY for parameters that exist in InputTJLF but NOT in InputTGLF
    # These are either TJLF-specific parameters OR TGLF parameters missing from InputTGLF structure
    input_tjlf.FIND_EIGEN = true    # TGLF parameter but missing from InputTGLF structure
    input_tjlf.NXGRID = 16         # TJLF-specific parameter

    # check converison
    checkInput(input_tjlf)

    return input_tjlf
end

# Per-ky working copy taken inside the threaded ky loop on every solve (and every AD
# evaluation). Unrolled to typed field copies for the same reason as above. NOTE: this
# deliberately aliases the parent's Vector fields into the copy (`setfield!` shares the
# reference) — the ky loop relies on writing per-ky results into the shared spectra.
@generated function minimal_scalar_copy(inputs::InputTJLF{T}) where {T<:Real}
    assigns = [:(setfield!(local_inputs, $(QuoteNode(f)), getfield(inputs, $(QuoteNode(f)))))
               for f in fieldnames(InputTJLF)]
    quote
        local_inputs = InputTJLF{T}(inputs.NS, length(inputs.KY_SPECTRUM))
        $(assigns...)
        return local_inputs
    end
end

##########################################################
#       Get from tjlf_hermite
##########################################################

struct OutputHermite{T<:Real}
    x::Vector{T}
    wx::Vector{T}
    h::Matrix{T}
    _dvec::Matrix{T}
end

function OutputHermite(x, wx, h, nky::Int)
    _dvec = Matrix{typeof(wx[1])}(undef, size(wx, 1), nky)
    return OutputHermite(x, wx, h, _dvec)
end

##########################################################
#       Get from xgrid_functions_geo
##########################################################

struct OutputGeometry{T<:Real}

    kx0_e::T

    fts::Vector{T}

    kxx::Vector{T}

    wdx::Vector{T}
    wdpx::Vector{T}

    b0x::Vector{T}
    b2x::Vector{T}

    cx_tor_par::Vector{T}
    cx_tor_per::Vector{T}
    cx_par_par::Vector{T}

end

struct SaturationParameters{T<:Real}

    SAT_geo0::T
    SAT_geo1::T
    SAT_geo2::T

    y::Vector{T}
    R_unit::T
    B_unit::T
    q_unit::T

    R::Vector{T}
    Bp::Vector{T}
    Bt::Vector{T}

    Bt0::T
    grad_r0::T

    S_prime::Vector{T}
    kx_factor::Vector{T}

    B_geo::Vector{T}
    qrat_geo::Vector{T}

    sintheta_geo::Vector{T}
    costheta_geo::Vector{T}
    costheta_p_geo::Vector{T}

    theta::Vector{T}

end



##########################################################
#       Used for get_matrx()
##########################################################

mutable struct Ave{T<:Real} <: AbstractAve{T}
    kx::Matrix{T}

    wdh::Matrix{T}
    modwdh::Matrix{T}
    wdg::Matrix{T}
    modwdg::Matrix{T}

    gradB::Matrix{T}
    b0::Matrix{T}
    b0inv::Matrix{T}
    lnB::Matrix{T}
    p0::Matrix{T}
    p0inv::Matrix{T}
    bp::Matrix{T}
    bpinv::Matrix{T}

    c_par_par::Matrix{T}
    c_tor_par::Matrix{T}
    c_tor_per::Matrix{T}

    kpar::Matrix{T}
    modkpar::Matrix{T}
    kpar_eff::Array{Complex{T},3}
    modkpar_eff::Array{Complex{T},3}
end

function Ave{T}(ns::Int, nbasis::Int) where {T<:Real}
    kx = zeros(T, nbasis, nbasis)
    wdh = zeros(T, nbasis, nbasis)
    modwdh = zeros(T, nbasis, nbasis)
    wdg = zeros(T, nbasis, nbasis)
    modwdg = zeros(T, nbasis, nbasis)
    gradB = zeros(T, nbasis, nbasis)
    b0 = zeros(T, nbasis, nbasis)
    b0inv = zeros(T, nbasis, nbasis)
    lnB = zeros(T, nbasis, nbasis)
    p0 = zeros(T, nbasis, nbasis)
    p0inv = zeros(T, nbasis, nbasis)
    bp = zeros(T, nbasis, nbasis)
    bpinv = zeros(T, nbasis, nbasis)
    c_par_par = zeros(T, nbasis, nbasis)
    c_tor_par = zeros(T, nbasis, nbasis)
    c_tor_per = zeros(T, nbasis, nbasis)
    kpar = zeros(T, nbasis, nbasis)
    modkpar = zeros(T, nbasis, nbasis)
    kpar_eff = zeros(Complex{T}, ns, nbasis, nbasis)
    modkpar_eff = zeros(Complex{T}, ns, nbasis, nbasis)

    Ave{T}(kx,
        wdh, modwdh, wdg, modwdg,
        gradB, b0, b0inv, lnB, p0, p0inv, bp, bpinv,
        c_par_par, c_tor_par, c_tor_per,
        kpar, modkpar, kpar_eff, modkpar_eff)
end


mutable struct AveH{T<:Real} <: AbstractAve{T}

    # average h-bessel functions
    hn::Array{T,3}
    hp1::Array{T,3}
    hp3::Array{T,3}
    hr11::Array{T,3}
    hr13::Array{T,3}
    hr33::Array{T,3}
    hw113::Array{T,3}
    hw133::Array{T,3}
    hw333::Array{T,3}

    ht1::Array{T,3}
    ht3::Array{T,3}
    hu1::Array{T,3}
    hu3::Array{T,3}
    hu33::Array{T,3}
    hu3ht1::Array{T,3}
    hu33ht1::Array{T,3}
    hu3ht3::Array{T,3}
    hu33ht3::Array{T,3}

    hnp0::Array{T,3}
    hp1p0::Array{T,3}
    hp3p0::Array{T,3}
    hr11p0::Array{T,3}
    hr13p0::Array{T,3}
    hr33p0::Array{T,3}
    c_tor_par_hp1p0::Array{T,3}
    c_tor_par_hr11p0::Array{T,3}
    c_tor_par_hr13p0::Array{T,3}

    hnb0::Array{T,3}
    hp1b0::Array{T,3}
    hp3b0::Array{T,3}
    hr11b0::Array{T,3}
    hr13b0::Array{T,3}
    hr33b0::Array{T,3}
    hw113b0::Array{T,3}
    hw133b0::Array{T,3}
    hw333b0::Array{T,3}

    hnbp::Array{T,3}
    hp1bp::Array{T,3}
    hp3bp::Array{T,3}
    hr11bp::Array{T,3}
    hr13bp::Array{T,3}
    hr33bp::Array{T,3}
    hw113bp::Array{T,3}
    hw133bp::Array{T,3}
    hw333bp::Array{T,3}
end

function AveH{T}(ns::Int, nbasis::Int) where {T<:Real}
    hn = zeros(T, ns, nbasis, nbasis)
    hp1 = zeros(T, ns, nbasis, nbasis)
    hp3 = zeros(T, ns, nbasis, nbasis)
    hr11 = zeros(T, ns, nbasis, nbasis)
    hr13 = zeros(T, ns, nbasis, nbasis)
    hr33 = zeros(T, ns, nbasis, nbasis)
    hw113 = zeros(T, ns, nbasis, nbasis)
    hw133 = zeros(T, ns, nbasis, nbasis)
    hw333 = zeros(T, ns, nbasis, nbasis)
    ht1 = zeros(T, ns, nbasis, nbasis)
    ht3 = zeros(T, ns, nbasis, nbasis)
    hu1 = zeros(T, ns, nbasis, nbasis)
    hu3 = zeros(T, ns, nbasis, nbasis)
    hu33 = zeros(T, ns, nbasis, nbasis)
    hu3ht1 = zeros(T, ns, nbasis, nbasis)
    hu33ht1 = zeros(T, ns, nbasis, nbasis)
    hu3ht3 = zeros(T, ns, nbasis, nbasis)
    hu33ht3 = zeros(T, ns, nbasis, nbasis)
    hnp0 = zeros(T, ns, nbasis, nbasis)
    hp1p0 = zeros(T, ns, nbasis, nbasis)
    hp3p0 = zeros(T, ns, nbasis, nbasis)
    hr11p0 = zeros(T, ns, nbasis, nbasis)
    hr13p0 = zeros(T, ns, nbasis, nbasis)
    hr33p0 = zeros(T, ns, nbasis, nbasis)
    c_tor_par_hp1p0 = zeros(T, ns, nbasis, nbasis)
    c_tor_par_hr11p0 = zeros(T, ns, nbasis, nbasis)
    c_tor_par_hr13p0 = zeros(T, ns, nbasis, nbasis)
    hnb0 = zeros(T, ns, nbasis, nbasis)
    hp1b0 = zeros(T, ns, nbasis, nbasis)
    hp3b0 = zeros(T, ns, nbasis, nbasis)
    hr11b0 = zeros(T, ns, nbasis, nbasis)
    hr13b0 = zeros(T, ns, nbasis, nbasis)
    hr33b0 = zeros(T, ns, nbasis, nbasis)
    hw113b0 = zeros(T, ns, nbasis, nbasis)
    hw133b0 = zeros(T, ns, nbasis, nbasis)
    hw333b0 = zeros(T, ns, nbasis, nbasis)
    hnbp = zeros(T, ns, nbasis, nbasis)
    hp1bp = zeros(T, ns, nbasis, nbasis)
    hp3bp = zeros(T, ns, nbasis, nbasis)
    hr11bp = zeros(T, ns, nbasis, nbasis)
    hr13bp = zeros(T, ns, nbasis, nbasis)
    hr33bp = zeros(T, ns, nbasis, nbasis)
    hw113bp = zeros(T, ns, nbasis, nbasis)
    hw133bp = zeros(T, ns, nbasis, nbasis)
    hw333bp = zeros(T, ns, nbasis, nbasis)

    AveH{T}(
        hn, hp1, hp3, hr11, hr13, hr33, hw113, hw133, hw333,
        ht1, ht3, hu1, hu3, hu33, hu3ht1, hu33ht1, hu3ht3, hu33ht3,
        hnp0, hp1p0, hp3p0, hr11p0, hr13p0, hr33p0, c_tor_par_hp1p0, c_tor_par_hr11p0, c_tor_par_hr13p0,
        hnb0, hp1b0, hp3b0, hr11b0, hr13b0, hr33b0, hw113b0, hw133b0, hw333b0,
        hnbp, hp1bp, hp3bp, hr11bp, hr13bp, hr33bp, hw113bp, hw133bp, hw333bp
    )
end

mutable struct AveWH{T<:Real} <: AbstractAve{T}

    wdhp1p0::Array{T,3}
    wdhr11p0::Array{T,3}
    wdhr13p0::Array{T,3}
    wdht1::Array{T,3}
    wdht3::Array{T,3}
    wdhu1::Array{T,3}
    wdhu3::Array{T,3}
    wdhu3ht1::Array{T,3}
    wdhu3ht3::Array{T,3}
    wdhu33::Array{T,3}
    wdhu33ht1::Array{T,3}
    wdhu33ht3::Array{T,3}

    modwdht1::Array{T,3}
    modwdht3::Array{T,3}
    modwdhu1::Array{T,3}
    modwdhu3::Array{T,3}
    modwdhu3ht1::Array{T,3}
    modwdhu3ht3::Array{T,3}
    modwdhu33::Array{T,3}
    modwdhu33ht1::Array{T,3}
    modwdhu33ht3::Array{T,3}

    wdhp1b0::Array{T,3}
    wdhr11b0::Array{T,3}
    wdhr13b0::Array{T,3}

    wdhp1bp::Array{T,3}
    wdhr11bp::Array{T,3}
    wdhr13bp::Array{T,3}
end

function AveWH{T}(ns::Int, nbasis::Int) where {T<:Real}
    wdhp1p0 = zeros(T, ns, nbasis, nbasis)
    wdhr11p0 = zeros(T, ns, nbasis, nbasis)
    wdhr13p0 = zeros(T, ns, nbasis, nbasis)
    wdht1 = zeros(T, ns, nbasis, nbasis)
    wdht3 = zeros(T, ns, nbasis, nbasis)
    wdhu1 = zeros(T, ns, nbasis, nbasis)
    wdhu3 = zeros(T, ns, nbasis, nbasis)
    wdhu3ht1 = zeros(T, ns, nbasis, nbasis)
    wdhu3ht3 = zeros(T, ns, nbasis, nbasis)
    wdhu33 = zeros(T, ns, nbasis, nbasis)
    wdhu33ht1 = zeros(T, ns, nbasis, nbasis)
    wdhu33ht3 = zeros(T, ns, nbasis, nbasis)

    modwdht1 = zeros(T, ns, nbasis, nbasis)
    modwdht3 = zeros(T, ns, nbasis, nbasis)
    modwdhu1 = zeros(T, ns, nbasis, nbasis)
    modwdhu3 = zeros(T, ns, nbasis, nbasis)
    modwdhu3ht1 = zeros(T, ns, nbasis, nbasis)
    modwdhu3ht3 = zeros(T, ns, nbasis, nbasis)
    modwdhu33 = zeros(T, ns, nbasis, nbasis)
    modwdhu33ht1 = zeros(T, ns, nbasis, nbasis)
    modwdhu33ht3 = zeros(T, ns, nbasis, nbasis)

    wdhp1b0 = zeros(T, ns, nbasis, nbasis)
    wdhr11b0 = zeros(T, ns, nbasis, nbasis)
    wdhr13b0 = zeros(T, ns, nbasis, nbasis)

    wdhp1bp = zeros(T, ns, nbasis, nbasis)
    wdhr11bp = zeros(T, ns, nbasis, nbasis)
    wdhr13bp = zeros(T, ns, nbasis, nbasis)

    AveWH{T}(
        wdhp1p0, wdhr11p0, wdhr13p0, wdht1, wdht3, wdhu1, wdhu3, wdhu3ht1, wdhu3ht3, wdhu33, wdhu33ht1, wdhu33ht3,
        modwdht1, modwdht3, modwdhu1, modwdhu3, modwdhu3ht1, modwdhu3ht3, modwdhu33, modwdhu33ht1, modwdhu33ht3,
        wdhp1b0, wdhr11b0, wdhr13b0,
        wdhp1bp, wdhr11bp, wdhr13bp
    )
end

mutable struct AveKH{K<:Complex} <: AbstractAve{K}

    kparhnp0::Array{K,3}
    kparhp1p0::Array{K,3}
    kparhp3p0::Array{K,3}
    kparhu1::Array{K,3}
    kparhu3::Array{K,3}
    kparht1::Array{K,3}
    kparht3::Array{K,3}
    modkparhu1::Array{K,3}
    modkparhu3::Array{K,3}

    kparhp1b0::Array{K,3}
    kparhr11b0::Array{K,3}
    kparhr13b0::Array{K,3}

    kparhnbp::Array{K,3}
    kparhp3bp::Array{K,3}
    kparhp1bp::Array{K,3}
    kparhr11bp::Array{K,3}
    kparhr13bp::Array{K,3}
end

function AveKH{K}(ns::Int, nbasis::Int) where {K<:Complex}
    kparhnp0 = zeros(K, ns, nbasis, nbasis)
    kparhp1p0 = zeros(K, ns, nbasis, nbasis)
    kparhp3p0 = zeros(K, ns, nbasis, nbasis)
    kparhu1 = zeros(K, ns, nbasis, nbasis)
    kparhu3 = zeros(K, ns, nbasis, nbasis)
    kparht1 = zeros(K, ns, nbasis, nbasis)
    kparht3 = zeros(K, ns, nbasis, nbasis)
    modkparhu1 = zeros(K, ns, nbasis, nbasis)
    modkparhu3 = zeros(K, ns, nbasis, nbasis)

    kparhp1b0 = zeros(K, ns, nbasis, nbasis)
    kparhr11b0 = zeros(K, ns, nbasis, nbasis)
    kparhr13b0 = zeros(K, ns, nbasis, nbasis)

    kparhnbp = zeros(K, ns, nbasis, nbasis)
    kparhp3bp = zeros(K, ns, nbasis, nbasis)
    kparhp1bp = zeros(K, ns, nbasis, nbasis)
    kparhr11bp = zeros(K, ns, nbasis, nbasis)
    kparhr13bp = zeros(K, ns, nbasis, nbasis)

    AveKH{K}(
        kparhnp0, kparhp1p0, kparhp3p0, kparhu1, kparhu3, kparht1, kparht3, modkparhu1, modkparhu3,
        kparhp1b0, kparhr11b0, kparhr13b0,
        kparhnbp, kparhp3bp, kparhp1bp, kparhr11bp, kparhr13bp
    )
end

# Convenience constructor defaulting to ComplexF64
AveKH(ns::Int, nbasis::Int) = AveKH{ComplexF64}(ns, nbasis)


mutable struct AveG{T<:Real} <: AbstractAve{T}

    # average g-bessel functions
    gn::Array{T,3}
    gp1::Array{T,3}
    gp3::Array{T,3}
    gr11::Array{T,3}
    gr13::Array{T,3}
    gr33::Array{T,3}
    gw113::Array{T,3}
    gw133::Array{T,3}
    gw333::Array{T,3}

    gt1::Array{T,3}
    gt3::Array{T,3}
    gu1::Array{T,3}
    gu3::Array{T,3}
    gu33::Array{T,3}
    gu3gt1::Array{T,3}
    gu3gt3::Array{T,3}
    gu33gt1::Array{T,3}
    gu33gt3::Array{T,3}

    gnp0::Array{T,3}
    gp1p0::Array{T,3}
    gp3p0::Array{T,3}
    gr11p0::Array{T,3}
    gr13p0::Array{T,3}
    gr33p0::Array{T,3}
    c_tor_par_gp1p0::Array{T,3}
    c_tor_par_gr11p0::Array{T,3}
    c_tor_par_gr13p0::Array{T,3}

    gnb0::Array{T,3}
    gp1b0::Array{T,3}
    gp3b0::Array{T,3}
    gr11b0::Array{T,3}
    gr13b0::Array{T,3}
    gr33b0::Array{T,3}
    gw113b0::Array{T,3}
    gw133b0::Array{T,3}
    gw333b0::Array{T,3}

    gnbp::Array{T,3}
    gp1bp::Array{T,3}
    gp3bp::Array{T,3}
    gr11bp::Array{T,3}
    gr13bp::Array{T,3}
    gr33bp::Array{T,3}
    gw113bp::Array{T,3}
    gw133bp::Array{T,3}
    gw333bp::Array{T,3}
end

function AveG{T}(ns::Int, nbasis::Int) where {T<:Real}
    gn = zeros(T, ns, nbasis, nbasis)
    gp1 = zeros(T, ns, nbasis, nbasis)
    gp3 = zeros(T, ns, nbasis, nbasis)
    gr11 = zeros(T, ns, nbasis, nbasis)
    gr13 = zeros(T, ns, nbasis, nbasis)
    gr33 = zeros(T, ns, nbasis, nbasis)
    gw113 = zeros(T, ns, nbasis, nbasis)
    gw133 = zeros(T, ns, nbasis, nbasis)
    gw333 = zeros(T, ns, nbasis, nbasis)

    gt1 = zeros(T, ns, nbasis, nbasis)
    gt3 = zeros(T, ns, nbasis, nbasis)
    gu1 = zeros(T, ns, nbasis, nbasis)
    gu3 = zeros(T, ns, nbasis, nbasis)
    gu33 = zeros(T, ns, nbasis, nbasis)
    gu3gt1 = zeros(T, ns, nbasis, nbasis)
    gu3gt3 = zeros(T, ns, nbasis, nbasis)
    gu33gt1 = zeros(T, ns, nbasis, nbasis)
    gu33gt3 = zeros(T, ns, nbasis, nbasis)

    gnp0 = zeros(T, ns, nbasis, nbasis)
    gp1p0 = zeros(T, ns, nbasis, nbasis)
    gp3p0 = zeros(T, ns, nbasis, nbasis)
    gr11p0 = zeros(T, ns, nbasis, nbasis)
    gr13p0 = zeros(T, ns, nbasis, nbasis)
    gr33p0 = zeros(T, ns, nbasis, nbasis)
    c_tor_par_gp1p0 = zeros(T, ns, nbasis, nbasis)
    c_tor_par_gr11p0 = zeros(T, ns, nbasis, nbasis)
    c_tor_par_gr13p0 = zeros(T, ns, nbasis, nbasis)

    gnb0 = zeros(T, ns, nbasis, nbasis)
    gp1b0 = zeros(T, ns, nbasis, nbasis)
    gp3b0 = zeros(T, ns, nbasis, nbasis)
    gr11b0 = zeros(T, ns, nbasis, nbasis)
    gr13b0 = zeros(T, ns, nbasis, nbasis)
    gr33b0 = zeros(T, ns, nbasis, nbasis)
    gw113b0 = zeros(T, ns, nbasis, nbasis)
    gw133b0 = zeros(T, ns, nbasis, nbasis)
    gw333b0 = zeros(T, ns, nbasis, nbasis)

    gnbp = zeros(T, ns, nbasis, nbasis)
    gp1bp = zeros(T, ns, nbasis, nbasis)
    gp3bp = zeros(T, ns, nbasis, nbasis)
    gr11bp = zeros(T, ns, nbasis, nbasis)
    gr13bp = zeros(T, ns, nbasis, nbasis)
    gr33bp = zeros(T, ns, nbasis, nbasis)
    gw113bp = zeros(T, ns, nbasis, nbasis)
    gw133bp = zeros(T, ns, nbasis, nbasis)
    gw333bp = zeros(T, ns, nbasis, nbasis)

    AveG{T}(
        gn, gp1, gp3, gr11, gr13, gr33, gw113, gw133, gw333,
        gt1, gt3, gu1, gu3, gu33, gu3gt1, gu3gt3, gu33gt1, gu33gt3,
        gnp0, gp1p0, gp3p0, gr11p0, gr13p0, gr33p0, c_tor_par_gp1p0, c_tor_par_gr11p0, c_tor_par_gr13p0,
        gnb0, gp1b0, gp3b0, gr11b0, gr13b0, gr33b0, gw113b0, gw133b0, gw333b0,
        gnbp, gp1bp, gp3bp, gr11bp, gr13bp, gr33bp, gw113bp, gw133bp, gw333bp
    )
end

mutable struct AveWG{T<:Real} <: AbstractAve{T}

    wdgp1p0::Array{T,3}
    wdgr11p0::Array{T,3}
    wdgr13p0::Array{T,3}
    wdgu1::Array{T,3}
    wdgu3::Array{T,3}
    wdgu33::Array{T,3}
    wdgt1::Array{T,3}
    wdgt3::Array{T,3}
    wdgu3gt1::Array{T,3}
    wdgu3gt3::Array{T,3}
    wdgu33gt1::Array{T,3}
    wdgu33gt3::Array{T,3}

    modwdgu1::Array{T,3}
    modwdgu3::Array{T,3}
    modwdgu33::Array{T,3}
    modwdgt1::Array{T,3}
    modwdgt3::Array{T,3}
    modwdgu3gt1::Array{T,3}
    modwdgu3gt3::Array{T,3}
    modwdgu33gt1::Array{T,3}
    modwdgu33gt3::Array{T,3}

    wdgp1b0::Array{T,3}
    wdgr11b0::Array{T,3}
    wdgr13b0::Array{T,3}

    wdgp1bp::Array{T,3}
    wdgr11bp::Array{T,3}
    wdgr13bp::Array{T,3}
end

function AveWG{T}(ns::Int, nbasis::Int) where {T<:Real}
    wdgp1p0 = zeros(T, ns, nbasis, nbasis)
    wdgr11p0 = zeros(T, ns, nbasis, nbasis)
    wdgr13p0 = zeros(T, ns, nbasis, nbasis)
    wdgu1 = zeros(T, ns, nbasis, nbasis)
    wdgu3 = zeros(T, ns, nbasis, nbasis)
    wdgu33 = zeros(T, ns, nbasis, nbasis)
    wdgt1 = zeros(T, ns, nbasis, nbasis)
    wdgt3 = zeros(T, ns, nbasis, nbasis)
    wdgu3gt1 = zeros(T, ns, nbasis, nbasis)
    wdgu3gt3 = zeros(T, ns, nbasis, nbasis)
    wdgu33gt1 = zeros(T, ns, nbasis, nbasis)
    wdgu33gt3 = zeros(T, ns, nbasis, nbasis)

    modwdgu1 = zeros(T, ns, nbasis, nbasis)
    modwdgu3 = zeros(T, ns, nbasis, nbasis)
    modwdgu33 = zeros(T, ns, nbasis, nbasis)
    modwdgt1 = zeros(T, ns, nbasis, nbasis)
    modwdgt3 = zeros(T, ns, nbasis, nbasis)
    modwdgu3gt1 = zeros(T, ns, nbasis, nbasis)
    modwdgu3gt3 = zeros(T, ns, nbasis, nbasis)
    modwdgu33gt1 = zeros(T, ns, nbasis, nbasis)
    modwdgu33gt3 = zeros(T, ns, nbasis, nbasis)

    wdgp1b0 = zeros(T, ns, nbasis, nbasis)
    wdgr11b0 = zeros(T, ns, nbasis, nbasis)
    wdgr13b0 = zeros(T, ns, nbasis, nbasis)

    wdgp1bp = zeros(T, ns, nbasis, nbasis)
    wdgr11bp = zeros(T, ns, nbasis, nbasis)
    wdgr13bp = zeros(T, ns, nbasis, nbasis)

    AveWG{T}(
        wdgp1p0, wdgr11p0, wdgr13p0, wdgu1, wdgu3, wdgu33, wdgt1, wdgt3, wdgu3gt1, wdgu3gt3, wdgu33gt1, wdgu33gt3,
        modwdgu1, modwdgu3, modwdgu33, modwdgt1, modwdgt3, modwdgu3gt1, modwdgu3gt3, modwdgu33gt1, modwdgu33gt3,
        wdgp1b0, wdgr11b0, wdgr13b0,
        wdgp1bp, wdgr11bp, wdgr13bp
    )
end

mutable struct AveKG{K<:Complex} <: AbstractAve{K}

    kpargnp0::Array{K,3}
    kpargp1p0::Array{K,3}
    kpargp3p0::Array{K,3}
    kpargu1::Array{K,3}
    kpargu3::Array{K,3}
    kpargt1::Array{K,3}
    kpargt3::Array{K,3}
    modkpargu1::Array{K,3}
    modkpargu3::Array{K,3}

    kpargp1b0::Array{K,3}
    kpargr11b0::Array{K,3}
    kpargr13b0::Array{K,3}

    kpargnbp::Array{K,3}
    kpargp3bp::Array{K,3}
    kpargp1bp::Array{K,3}
    kpargr11bp::Array{K,3}
    kpargr13bp::Array{K,3}
end

function AveKG{K}(ns::Int, nbasis::Int) where {K<:Complex}
    kpargnp0 = zeros(K, ns, nbasis, nbasis)
    kpargp1p0 = zeros(K, ns, nbasis, nbasis)
    kpargp3p0 = zeros(K, ns, nbasis, nbasis)
    kpargu1 = zeros(K, ns, nbasis, nbasis)
    kpargu3 = zeros(K, ns, nbasis, nbasis)
    kpargt1 = zeros(K, ns, nbasis, nbasis)
    kpargt3 = zeros(K, ns, nbasis, nbasis)
    modkpargu1 = zeros(K, ns, nbasis, nbasis)
    modkpargu3 = zeros(K, ns, nbasis, nbasis)

    kpargp1b0 = zeros(K, ns, nbasis, nbasis)
    kpargr11b0 = zeros(K, ns, nbasis, nbasis)
    kpargr13b0 = zeros(K, ns, nbasis, nbasis)

    kpargnbp = zeros(K, ns, nbasis, nbasis)
    kpargp3bp = zeros(K, ns, nbasis, nbasis)
    kpargp1bp = zeros(K, ns, nbasis, nbasis)
    kpargr11bp = zeros(K, ns, nbasis, nbasis)
    kpargr13bp = zeros(K, ns, nbasis, nbasis)

    AveKG{K}(
        kpargnp0, kpargp1p0, kpargp3p0, kpargu1, kpargu3, kpargt1, kpargt3, modkpargu1, modkpargu3,
        kpargp1b0, kpargr11b0, kpargr13b0,
        kpargnbp, kpargp3bp, kpargp1bp, kpargr11bp, kpargr13bp
    )
end

# Convenience constructor defaulting to ComplexF64
AveKG(ns::Int, nbasis::Int) = AveKG{ComplexF64}(ns, nbasis)


mutable struct AveGrad{T<:Real}  <: AbstractAve{T}
    gradhp1::Array{T,3}
    gradhr11::Array{T,3}
    gradhr13::Array{T,3}
    gradhp1p1::Array{T,3}
    gradhr11p1::Array{T,3}
    gradhr13p1::Array{T,3}
    gradhp1p0::Array{T,3}
    gradhr11p0::Array{T,3}
    gradhr13p0::Array{T,3}

    gradgp1::Array{T,3}
    gradgr11::Array{T,3}
    gradgr13::Array{T,3}
    gradgp1p1::Array{T,3}
    gradgr11p1::Array{T,3}
    gradgr13p1::Array{T,3}
    gradgp1p0::Array{T,3}
    gradgr11p0::Array{T,3}
    gradgr13p0::Array{T,3}
end

function AveGrad{T}(ns::Int, nbasis::Int) where {T<:Real}
    gradhp1 = zeros(T, ns, nbasis, nbasis)
    gradhr11 = zeros(T, ns, nbasis, nbasis)
    gradhr13 = zeros(T, ns, nbasis, nbasis)
    gradhp1p1 = zeros(T, ns, nbasis, nbasis)
    gradhr11p1 = zeros(T, ns, nbasis, nbasis)
    gradhr13p1 = zeros(T, ns, nbasis, nbasis)
    gradhp1p0 = zeros(T, ns, nbasis, nbasis)
    gradhr11p0 = zeros(T, ns, nbasis, nbasis)
    gradhr13p0 = zeros(T, ns, nbasis, nbasis)

    gradgp1 = zeros(T, ns, nbasis, nbasis)
    gradgr11 = zeros(T, ns, nbasis, nbasis)
    gradgr13 = zeros(T, ns, nbasis, nbasis)
    gradgp1p1 = zeros(T, ns, nbasis, nbasis)
    gradgr11p1 = zeros(T, ns, nbasis, nbasis)
    gradgr13p1 = zeros(T, ns, nbasis, nbasis)
    gradgp1p0 = zeros(T, ns, nbasis, nbasis)
    gradgr11p0 = zeros(T, ns, nbasis, nbasis)
    gradgr13p0 = zeros(T, ns, nbasis, nbasis)

    AveGrad{T}(
        gradhp1, gradhr11, gradhr13, gradhp1p1, gradhr11p1, gradhr13p1, gradhp1p0, gradhr11p0, gradhr13p0,
        gradgp1, gradgr11, gradgr13, gradgp1p1, gradgr11p1, gradgr13p1, gradgp1p0, gradgr11p0, gradgr13p0
    )
end

mutable struct AveGradB{T<:Real}  <: AbstractAve{T}
    gradBhp1::Array{T,3}
    gradBhp3::Array{T,3}
    gradBhr11::Array{T,3}
    gradBhr13::Array{T,3}
    gradBhr33::Array{T,3}
    gradBhu1::Array{T,3}
    gradBhu3::Array{T,3}
    gradBhu33::Array{T,3}

    gradBgp1::Array{T,3}
    gradBgp3::Array{T,3}
    gradBgr11::Array{T,3}
    gradBgr13::Array{T,3}
    gradBgr33::Array{T,3}
    gradBgu1::Array{T,3}
    gradBgu3::Array{T,3}
    gradBgu33::Array{T,3}
end

function AveGradB{T}(ns::Int, nbasis::Int) where {T<:Real}
    gradBhp1 = zeros(T, ns, nbasis, nbasis)
    gradBhp3 = zeros(T, ns, nbasis, nbasis)
    gradBhr11 = zeros(T, ns, nbasis, nbasis)
    gradBhr13 = zeros(T, ns, nbasis, nbasis)
    gradBhr33 = zeros(T, ns, nbasis, nbasis)
    gradBhu1 = zeros(T, ns, nbasis, nbasis)
    gradBhu3 = zeros(T, ns, nbasis, nbasis)
    gradBhu33 = zeros(T, ns, nbasis, nbasis)

    gradBgp1 = zeros(T, ns, nbasis, nbasis)
    gradBgp3 = zeros(T, ns, nbasis, nbasis)
    gradBgr11 = zeros(T, ns, nbasis, nbasis)
    gradBgr13 = zeros(T, ns, nbasis, nbasis)
    gradBgr33 = zeros(T, ns, nbasis, nbasis)
    gradBgu1 = zeros(T, ns, nbasis, nbasis)
    gradBgu3 = zeros(T, ns, nbasis, nbasis)
    gradBgu33 = zeros(T, ns, nbasis, nbasis)

    AveGradB{T}(
        gradBhp1, gradBhp3, gradBhr11, gradBhr13, gradBhr33, gradBhu1, gradBhu3, gradBhu33,
        gradBgp1, gradBgp3, gradBgr11, gradBgr13, gradBgr33, gradBgu1, gradBgu3, gradBgu33
    )
end

function reset_ave!(ave::AbstractAve{T}) where {T<:Number}
    for fieldname in fieldnames(typeof(ave))
        fill!(getproperty(ave, fieldname), zero(T))
    end
end

const AvesTuple{T} = Tuple{Ave{T},AveH{T},AveWH{T},AveKH{Complex{T}},AveG{T},AveWG{T},AveKG{Complex{T}},AveGrad{T},AveGradB{T}}

"""
    LSWorkspace{T}()

Reusable buffers for repeated single-ky linear solves (`TJLF.run(input; workspace)` /
`tjlf_LS(...; workspace)`): the dense `iur × iur` eigenmatrices `amat`/`bmat`, the two
copies the LAPACK solve destroys, the inverse-iteration matrix `zmat`, the eigenvector,
and the `Ave*` structs. Without it every call allocates these afresh (five complex
`iur × iur` matrices, ~165 MB at 3 species / NBASIS=32), a heavy load on the garbage
collector when many threads solve at once.

The buffers are sized on first use and resized when `(NS, nbasis)` changes. A workspace
must not be shared between concurrently running tasks: keep one per task. Results are
bit-identical to a call without a workspace (same operations on the same values; only
the storage is reused).
"""
mutable struct LSWorkspace{T<:Real}
    ns::Int
    nbasis::Int
    iur::Int
    amat::Matrix{Complex{T}}
    bmat::Matrix{Complex{T}}
    amat_copy::Matrix{Complex{T}}
    bmat_copy::Matrix{Complex{T}}
    zmat::Matrix{Complex{T}}
    eigenvector::Vector{Complex{T}}
    aves::Union{Nothing,AvesTuple{T}}
end

function LSWorkspace{T}() where {T<:Real}
    K = Complex{T}
    e() = Matrix{K}(undef, 0, 0)
    return LSWorkspace{T}(0, 0, 0, e(), e(), e(), e(), e(), Vector{K}(undef, 0), nothing)
end
LSWorkspace() = LSWorkspace{Float64}()
