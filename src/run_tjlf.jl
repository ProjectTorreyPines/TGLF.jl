"""
    run(input::Union{InputTJLF,InputTGLF}; use_gpu::Bool=false)

Run the full TGLF quasi-linear transport model and return everything it computes
as a `NamedTuple` (unlike [`run_tjlf`](@ref), which returns only the summed
fluxes). An `InputTGLF` is converted to an `InputTJLF{Float64}` first.

parameters:
    input::InputTJLF | InputTGLF            - the local-equilibrium input
    use_gpu::Bool                           - offload the dense eigensolves to CUDA (default false)

outputs (NamedTuple):
    QL_weights       - quasi-linear weights      [field, species, mode, ky, type]
    eigenvalue       - growth rate & frequency   [mode, ky, type]
    QL_flux_out      - summed fluxes             [field, species, type]   (== `run_tjlf` output)
    flux_spectrum    - flux per ky               [field, species, mode, ky, type]
    field_weight_out - mode (eigen) structure    [field, basis, mode, ky]

description:
    Wraps the solve in `with_blas_threads(1)` so the threaded `ky` loop does not
    oversubscribe BLAS. With `USE_TRANSPORT_MODEL=true` (the usual path) it runs
    the full saturation-rule spectrum sum; otherwise it returns a single-`ky`
    linear solution with empty `QL_flux_out`/`flux_spectrum`. As with `run_tjlf`,
    the width and eigenvalue spectra found during the run are saved back into the
    input struct for reuse via the `FIND_WIDTH`/`FIND_EIGEN` flags.

    For repeated single-`ky` solves (`USE_TRANSPORT_MODEL=false`) pass
    `workspace = LSWorkspace()` (one per task) to reuse the large eigenmatrix
    buffers between calls; results are unchanged.

!!! note
    `run` is exported but shadows `Base.run`, so call it as `TJLF.run`.
"""
function run(inputTJLF::InputTJLF; use_gpu::Bool=false, workspace::Union{Nothing,LSWorkspace}=nothing)
    return with_blas_threads(1) do
        _run_impl(inputTJLF; use_gpu=use_gpu, workspace=workspace)
    end
end

function _run_impl(inputTJLF::InputTJLF; use_gpu::Bool=false, workspace::Union{Nothing,LSWorkspace}=nothing)
    use_tm = inputTJLF.USE_TRANSPORT_MODEL
    if use_tm
        apply_presets!(inputTJLF)
        checkInput(inputTJLF)
        outputHermite = gauss_hermite(inputTJLF)
        satParams = get_sat_params(inputTJLF)
        inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
        tm = tjlf_TM(inputTJLF, satParams, outputHermite; use_gpu=use_gpu)
        QL_weights            = tm.QL_weights
        firstPass_eigenvalue  = tm.firstPass_eigenvalue
        secondPass_eigenvalue = tm.secondPass_eigenvalue
        field_weight_out      = tm.field_weight_out
        phi_bar_matrix        = tm.phi_bar_matrix

        vexb_shear_s = inputTJLF.VEXB_SHEAR * inputTJLF.SIGN_IT
        if inputTJLF.ALPHA_QUENCH == 0.0 && vexb_shear_s != 0.0
            eigenvalue_for_flux = secondPass_eigenvalue
        else
            eigenvalue_for_flux = firstPass_eigenvalue
        end

        if inputTJLF.SAT_RULE == 2 || inputTJLF.SAT_RULE == 3
            most_unstable_gamma_first_pass = firstPass_eigenvalue[1, :, 1]
            vzf_out_first_pass, kymax_out_first_pass, jmax_out_first_pass = get_zonal_mixing(inputTJLF, satParams, most_unstable_gamma_first_pass)
            QL_flux_out, flux_spectrum = sum_ky_spectrum(inputTJLF, satParams, eigenvalue_for_flux[:, :, 1], QL_weights;
                                                         vzf_out_param=vzf_out_first_pass,
                                                         kymax_out_param=kymax_out_first_pass,
                                                         jmax_out_param=jmax_out_first_pass)
        else
            QL_flux_out, flux_spectrum = sum_ky_spectrum(inputTJLF, satParams, eigenvalue_for_flux[:, :, 1], QL_weights;
                                                         phi_bar_matrix=phi_bar_matrix)
        end
        return (QL_weights=QL_weights, eigenvalue=firstPass_eigenvalue, QL_flux_out=QL_flux_out,
                flux_spectrum=flux_spectrum, field_weight_out=field_weight_out)
    else
        checkInput(inputTJLF)
        outputHermite = gauss_hermite(inputTJLF)
        satParams = get_sat_params(inputTJLF)
        inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
        ns     = inputTJLF.NS
        nmodes = inputTJLF.NMODES
        nbasis = inputTJLF.NBASIS_MAX
        nmodes_out, gamma_out, freq_out,
            particle_QL_out, energy_QL_out, stress_tor_QL_out, stress_par_QL_out, exchange_QL_out,
            _ft, field_weight_out_3d, _phi = tjlf_LS(inputTJLF, satParams, outputHermite,
                                                      inputTJLF.KY, nbasis, inputTJLF.VEXB_SHEAR, 1;
                                                      use_gpu=use_gpu, workspace=workspace)

        eigenvalue       = zeros(Float64, nmodes, 1, 2)
        QL_weights       = zeros(Float64, 3, ns, nmodes, 1, 5)
        field_weight_out = zeros(ComplexF64, 3, nbasis, nmodes, 1)

        eigenvalue[:, 1, 1]                      .= gamma_out
        eigenvalue[:, 1, 2]                      .= freq_out
        QL_weights[:, :, 1:nmodes_out, 1, 1]    .= particle_QL_out[:, :, 1:nmodes_out]
        QL_weights[:, :, 1:nmodes_out, 1, 2]    .= energy_QL_out[:, :, 1:nmodes_out]
        QL_weights[:, :, 1:nmodes_out, 1, 3]    .= stress_tor_QL_out[:, :, 1:nmodes_out]
        QL_weights[:, :, 1:nmodes_out, 1, 4]    .= stress_par_QL_out[:, :, 1:nmodes_out]
        QL_weights[:, :, 1:nmodes_out, 1, 5]    .= exchange_QL_out[:, :, 1:nmodes_out]
        field_weight_out[:, :, 1:nmodes_out, 1] .= field_weight_out_3d[:, :, 1:nmodes_out]

        QL_flux_out   = Array{Float64}(undef, 0, 0, 0)
        flux_spectrum = Array{Float64}(undef, 0, 0, 0)
        return (QL_weights=QL_weights, eigenvalue=eigenvalue, QL_flux_out=QL_flux_out,
                flux_spectrum=flux_spectrum, field_weight_out=field_weight_out)
    end
end

function run(inputTGLF::InputTGLF; use_gpu::Bool=false)
    inputTJLF = InputTJLF{Float64}(inputTGLF)
    return run(inputTJLF; use_gpu=use_gpu)
end

"""
    run_tjlf(inputTJLF::InputTJLF)

parameters:
    input_tjlfs::InputTJLF                  - InputTJLF struct

outputs:
    outputs                                 - fluxes (field, species, type)

description:
    Runs TJLF on a single InputTJLF struct, during the run, will save the width spectrum and eigenvalue spectrum to the InputTJLF struct.
    If you want to use these widths and eigenvalues in future runs, there is a flag: FIND_WIDTH and FIND_EIGEN that you set to false
"""
function run_tjlf(inputTJLF::InputTJLF; use_gpu::Bool=false)
    return run(inputTJLF; use_gpu=use_gpu).QL_flux_out
end

"""
    run_tjlf(input_tjlfs::Vector{InputTJLF{T}}) where {T<:Real}

parameters:
    input_tjlfs::Vector{InputTJLF{T}}       - vector of InputTJLF structs

outputs:
    outputs                                 - vector of fluxes (field, species, type)

description:
    Runs TJLF on a vector of InputTJLF structs, during the run, will save the width spectrum and eigenvalue spectrum to the InputTJLF struct.
    If you want to use these widths and eigenvalues in future runs, there is a flag: FIND_WIDTH and FIND_EIGEN that you set to false
"""
function run_tjlf(input_tjlfs::Vector{InputTJLF{T}}; use_gpu::Bool=false) where {T<:Real}
    checkInput(input_tjlfs)
    outputs = Vector{Array{T,3}}(undef, length(input_tjlfs))
    # set BLAS=1 before spawning threads; the per-element run() short-circuits
    with_blas_threads(1) do
        Threads.@threads for idx in eachindex(input_tjlfs)
            outputs[idx] = TJLF.run_tjlf(input_tjlfs[idx]; use_gpu=use_gpu)
        end
    end
    return outputs
end

# (field, species, type)
# type: (particle, energy, torodial stress, parallel stress, exchange)

Qe(QL_flux_out::Array{<:Real}) = sum(QL_flux_out[:, 1, 2])

Qi(QL_flux_out::Array{<:Real}) = sum(QL_flux_out[:, 2:end, 2])

Πi(QL_flux_out::Array{<:Real}) = sum(QL_flux_out[:, 2:end, 3])

Γe(QL_flux_out::Array{<:Real}) = sum(QL_flux_out[:, 1, 1])

Γi(QL_flux_out::Array{<:Real}) = [sum(QL_flux_out[:, k, 1]) for k in 2:size(QL_flux_out)[2]]


"""
    pick_device(mode::Symbol)

parameters:
    mode::Symbol - :cpu, :gpu, or :auto

outputs:
    device::Symbol - :cpu or :gpu

description:
    Determines the device to use based on the mode and CUDA availability.
"""
function pick_device(mode::Symbol)
    if mode === :cpu
        return :cpu
    elseif mode === :gpu
        (_cuda_functional() && _cuda_device_count() > 0) || error("GPU requested, but CUDA is not functional")
        return :gpu
    elseif mode === :auto
        return _cuda_functional() ? :gpu : :cpu
    else
        error("Unknown mode: $mode")
    end
end

