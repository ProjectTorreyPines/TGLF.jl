using Test
using TJLF

# saturation rule test
satRuleDirectory = joinpath(@__DIR__, "test_SAT_rules")
tests = readdir(satRuleDirectory)
for dir_name in tests
    if dir_name == ".DS_Store"
        continue
    end
    baseDirectory = joinpath(satRuleDirectory, dir_name)

    @testset "$baseDirectory" begin

        fileDirectory = joinpath(baseDirectory, "out.tglf.QL_flux_spectrum")
        lines = readlines(fileDirectory)
        # 5, 2, 3, 21, 2
        (ntype, nspecies, nfield, nky, nmodes) = parse.(Int32, split(lines[4]))
        ql = Vector{Float64}()
        for line in lines[7:length(lines)]
            line = split(line)
            if any(occursin.(["m", "s"], string(line)))
                continue
            end

            for x in line
                push!(ql, parse(Float64, string(x)))
            end

        end
        # 5, 21, 2, 3, 2
        QLw = reshape(ql, (ntype, nky, nmodes, nfield, nspecies))
        QL_data = permutedims(QLw, (4, 5, 3, 2, 1))
        # QL = AxisArray(QL_data, Axis{:nky}(1:nky), Axis{:nmodes}(1:nmodes), Axis{:nspecies}(1:nspecies), Axis{:nfield}(1:nfield), Axis{:ntype}(("particle","energy","toroidal_stress","parrllel_stress","exchange_stress")))

        particle_QL = QL_data[:, :, :, :, 1]
        energy_QL = QL_data[:, :, :, :, 2]
        toroidal_stress_QL = QL_data[:, :, :, :, 3]
        parallel_stress_QL = QL_data[:, :, :, :, 4]
        exchange_QL = QL_data[:, :, :, :, 5]

        # Read spectral shift and ave_p0 (only needed for SAT0)
        fileDirectory = joinpath(baseDirectory, "out.tglf.spectral_shift_spectrum")
        lines = readlines(fileDirectory)
        kx0_e = Vector{Float64}()
        for line in lines[6:length(lines)]
            push!(kx0_e, parse(Float64, line))
        end

        fileDirectory = joinpath(baseDirectory, "out.tglf.ave_p0_spectrum")
        lines = readlines(fileDirectory)
        ave_p0 = Vector{Float64}()
        for line in lines[4:length(lines)]
            push!(ave_p0, parse(Float64, line))
        end

        # Read scalar saturation parameters
        fileDirectory = joinpath(baseDirectory, "out.tglf.scalar_saturation_parameters")
        lines = readlines(fileDirectory)
        inputComparison = Dict()
        for line in lines
            #### no idea why this is here
            if any(occursin.(["!", "UNITS", "SAT_RULE", "XNU_MODEL", "ETG_FACTOR", "ALPHA_ZF", "RULE"], string(line)))
                continue
            end
            field, value = strip.(split(line, "="))
            if contains(line, "R_unit")
                R_unit = parse(Float64, value)
            else
                inputComparison[string(field)] = parse(Float64, value)
            end
        end

        # Get ky spectrum
        fileDirectory = joinpath(baseDirectory, "out.tglf.ky_spectrum")
        lines = readlines(fileDirectory)
        ky_spect = Array{Float64}(undef, 0)
        for line in lines[3:length(lines)]
            push!(ky_spect, parse(Float64, line))
        end

        # Get eigenvalue spectrum
        fileDirectory = joinpath(baseDirectory, "out.tglf.eigenvalue_spectrum")
        lines = readlines(fileDirectory)
        lines = split(join(lines[3:length(lines)]))
        lines = [parse(Float64, l) for l in lines]

        gamma = []
        freq = []
        for k in 1:nmodes
            push!(gamma, lines[2k-1:2*nmodes:end])
            push!(freq, lines[2k-1:2*nmodes:end])
        end

        gammas = hcat(gamma...)
        R_unit = ones(size(gammas)) * R_unit

        # Get potential spectrum
        fileDirectory = joinpath(baseDirectory, "out.tglf.field_spectrum")
        lines = readlines(fileDirectory)

        columns = split.(lines[2], ",")
        nc = length(columns)

        lines = split(join(lines[6:length(lines)]))

        tmpdict = Dict()
        for (ik, k) in enumerate(columns)
            tmp = []
            for nm in 1:nmodes
                push!(tmp, parse.(Float64, lines[ik-3+nm*nc:nc*nmodes:end]))
            end
            tmpdict[k] = tmp
        end
        potentialTmp = tmpdict[columns[length(columns)]]
        potential = hcat(potentialTmp...)

        fileDirectory = joinpath(baseDirectory, "out.tglf.gbflux")
        lines = readlines(fileDirectory)
        width::Integer = round(length(split(lines[1])) / 4)
        fluxes = transpose(reshape(parse.(Float64, split(lines[1])), (width, 4)))

        #******************************************************************************#************************
        # Read input.tglf
        #******************************************************************************#************************

        inputTJLF = readInput(joinpath(baseDirectory, "input.tglf"))

        #*******************************************************************************************************
        #   start running stuff
        #*******************************************************************************************************

        satParams = get_sat_params(inputTJLF)

        inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
        @test isapprox(ky_spect, inputTJLF.KY_SPECTRUM, rtol=1e-6)

        kx0epy = xgrid_functions_geo(inputTJLF, satParams, Matrix(gammas'))
        @test isapprox(kx0epy, kx0_e, rtol=1e-3)
        @test isapprox(inputComparison["SAT_geo0_out"], satParams.SAT_geo0, rtol=1e-6)
        @test isapprox(inputComparison["SAT_geo1_out"], satParams.SAT_geo1, rtol=1e-6)
        @test isapprox(inputComparison["SAT_geo2_out"], satParams.SAT_geo2, rtol=1e-5)
        @test isapprox(R_unit[1, 1], satParams.R_unit, rtol=1e-6)
        @test isapprox(inputComparison["Bt0_out"], satParams.Bt0, rtol=1e-6)
        @test isapprox(inputComparison["grad_r0_out"], satParams.grad_r0, rtol=1e-6)

        if inputTJLF.VEXB_SHEAR != 0.0
            @test isapprox(inputComparison["B_geo0_out"], satParams.B_geo[1], rtol=1e-6)
        end

        QL_flux_out, flux_out = sum_ky_spectrum(inputTJLF, satParams, Matrix(gammas'), QL_data)
        # julia_sat1 = sum(sum(sat_1["energy_flux_integral"], dims=3)[:,:,1], dims=1)[1,:]
        julia_sat1 = sum(QL_flux_out; dims=1)[1, :, 2]
        expected_sat1 = fluxes[2, :]
        @test isapprox(julia_sat1, expected_sat1, rtol=1e-3)
    end
end


# SAT3 with the spectral peak on a wavenumber of the grid, and a mode that carries no weight at
# the next wavenumber: the weight ratio at the peak must come from the two wavenumbers bracketing
# it, where the mode exists, and the fluxes must be finite
@testset "SAT3 weight ratio at a grid wavenumber" begin
    @test TJLF.linear_interpolation([1.0, 2.0, 3.0], [1.0, 2.0, NaN], 2.0) == 2.0
    @test TJLF.linear_interpolation([1.0, 2.0, 3.0], [1.0, 2.0, 4.0], 2.5) == 3.0

    input = readInput(joinpath(@__DIR__, "test_SAT_rules", "SAT3", "input.tglf"))
    out = TJLF.run(input)
    satParams = TJLF.get_sat_params(input)
    gammas = out.eigenvalue[:, :, 1]
    vzf, _, jmax = TJLF.get_zonal_mixing(input, satParams, gammas[1, :])
    ky = input.KY_SPECTRUM
    @test 1 < jmax < length(ky)
    sum_fluxes(W) = TJLF.sum_ky_spectrum(input, satParams, gammas, W; vzf_out_param=vzf, kymax_out_param=ky[jmax], jmax_out_param=jmax)[1]

    W = copy(out.QL_weights)
    mode = size(W, 3)
    W[:, :, mode, jmax-1:jmax, :] .= W[:, :, 1, jmax-1:jmax, :]
    reference = sum_fluxes(W)
    @test all(isfinite, reference)

    W[:, :, mode, jmax+1, :] .= 0.0
    @test all(isfinite, sum_fluxes(W))
end
