using Test
using TJLF

testTMDirectory = joinpath(@__DIR__, "test_TM")
tests = readdir(testTMDirectory)
for dir_name in tests
    if dir_name == ".DS_Store"
        continue
    end
    baseDirectory = joinpath(testTMDirectory, dir_name)

    @testset "$baseDirectory" begin

        # get particle, energy, exchange fluxes
        fileDirectory = joinpath(baseDirectory, "out.tglf.QL_flux_spectrum")
        lines = readlines(fileDirectory)
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
        QLw = reshape(ql, (ntype, nky, nmodes, nfield, nspecies))
        QL_data = permutedims(QLw, (4, 5, 3, 2, 1))

        # Get eigenvalue spectrum
        fileDirectory = joinpath(baseDirectory, "out.tglf.eigenvalue_spectrum")
        lines = readlines(fileDirectory)
        lines = split(join(lines[3:length(lines)]))
        lines = [parse(Float64, l) for l in lines]

        gamma = Vector{Vector{Float64}}()
        freq = Vector{Vector{Float64}}()
        for k in 1:nmodes
            push!(gamma, lines[2k-1:2*nmodes:end])
            push!(freq, lines[2k:2*nmodes:end])
        end
        gamma = vcat(gamma...)
        freq = vcat(freq...)

        # Get the integral flux
        fileDirectory = joinpath(baseDirectory, "out.tglf.gbflux")
        lines = readlines(fileDirectory)
        width::Integer = round(length(split(lines[1])) / 4)
        fluxesFortran = transpose(reshape(parse.(Float64, split(lines[1])), (width, 4)))

        #******************************************************************************************************
        # Read input.tglf
        #******************************************************************************************************

        inputTJLF = readInput(joinpath(baseDirectory, "input.tglf"))
        inputTJLF2 = readInput(joinpath(baseDirectory, "input.tglf"))

        #*******************************************************************************************************
        #   start running stuff
        #*******************************************************************************************************

        outputHermite = gauss_hermite(inputTJLF)
        satParams = get_sat_params(inputTJLF)
        inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)

        tm = tjlf_TM(inputTJLF, satParams, outputHermite)
        fluxes = tm.QL_weights
        eigenvalue = tm.firstPass_eigenvalue
        field_weight_out = tm.field_weight_out
        gammaJulia = eigenvalue[:, :, 1]
        freqJulia = eigenvalue[:, :, 2]

        # TGLF writes the first nfield fields; compare them element by element (eachindex checks the shapes)
        QL_julia = fluxes[1:nfield, :, :, :, :]
        for i in eachindex(QL_data, QL_julia)
            @test isapprox(QL_data[i], QL_julia[i], atol=1e-6)
        end
        for i in eachindex(gammaJulia)
            @test isapprox(gamma[i], gammaJulia'[i], atol=1e-6)
            @test isapprox(freq[i], freqJulia'[i], atol=1e-6)
        end

        QL_flux_out, flux_out = sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:, :, 1], fluxes)
        juliaEnergy = sum(QL_flux_out; dims=1)[1, :, 2]
        fortranEnergy = fluxesFortran[2, :]
        @test isapprox(juliaEnergy, fortranEnergy, rtol=1e-3)

        if dir_name == "simple_test"
            nky = length(inputTJLF.KY_SPECTRUM)
            @test size(field_weight_out) == (3, inputTJLF.NBASIS_MAX, inputTJLF.NMODES, nky)

            plot_field_out, plot_angle_out, max_plot, nmodes_out = get_wavefunction(inputTJLF, satParams, field_weight_out[:,:,:,1])
            @test nmodes_out == inputTJLF.NMODES
            @test size(plot_field_out) == (inputTJLF.NMODES, 3, max_plot)
            @test length(plot_angle_out) == max_plot
            @test any(!iszero, plot_field_out)
        end


        #*******************************************************************************************************0
        #   test that providing widths yields same results (NO LONGER VALID AFTER NEW QUICK/ITERATIVE EIGENSOLVER METHOD)
        #*******************************************************************************************************

        # inputTJLF2.KY_SPECTRUM .= inputTJLF.KY_SPECTRUM
        # inputTJLF2.WIDTH_SPECTRUM .= inputTJLF.WIDTH_SPECTRUM
        # inputTJLF2.EIGEN_SPECTRUM .= eigenvalue[inputTJLF.NMODES,:,1]
        # inputTJLF2.FIND_WIDTH = false

        # fluxes2, eigenvalue2 = tjlf_TM(inputTJLF2, satParams, outputHermite)

        # for i in eachindex(fluxes)
        #     @test isapprox(fluxes[i], fluxes2[i],atol=1e-8)
        # end

    end
end