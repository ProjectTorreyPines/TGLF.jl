"""
    tjlf_LS(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T},ky::T,nbasis::Int,vexb_shear_s::T;kx0_e::T = NaN,gamma_reference_kx0::Vector{T} = T[],freq_reference_kx0::Vector{T} = T[],outputGeo::Union{OutputGeometry{T},Missing} = missing) where T <: Real

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    satParams::SaturationParameters{T}  - SaturationParameters struct constructed in tjlf_geometry.jl
    outputHermite::OutputHermite{T}     - OutputHermite struct constructed in tjlf_hermite.jl
    ky::T                               - ky value
    nbasis::Int                         - number of basis for matrix dimension
    vexb_shear_s::T                     - e x b shear value (=VEXB_SHEAR*SIGN_IT)
    kx0_e::T = NaN                      - kx0_e value calculated on second pass with eigen values from first pass
    gamma_reference_kx0                 - gamma vector if second pass
    freq_reference_kx0                  - freq vector if second pass

outputs:
    nmodes_out                          - number of most unstable modes calculated
    gamma_out                           - array of growth rate eigenvalue calculated
    freq_out                            - array of frequency eigenvalue calculated
    particle_QL_out                     - particle QL flux
    energy_QL_out                       - energy QL flux
    stress_tor_QL_out                   - torodial stress QL flux
    stress_par_QL_out                   - parallel stress QL flux
    exchange_QL_out                     - exchange QL flux

description:
    TGLF Linear Stability driver computes all quasilinear quantities for a single ky
"""
function construct_linear_map(A, B, sigma)
    P = lu(A - sigma*B)    
    LinearMap{eltype(A)}((y,x) -> ldiv!(y, P, B*x),size(A, 1))
end

"""
    fit_workspace!(ws::LSWorkspace{T}, inputs::InputTJLF{T}, nbasis::Int)

Size (or resize) the workspace buffers for `inputs.NS` species and `nbasis` basis functions.
A no-op when they already fit. Returns `ws`.
"""
function fit_workspace!(ws::LSWorkspace{T}, inputs::InputTJLF{T}, nbasis::Int) where {T<:Real}
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    iur = (ns - ns0 + 1) * 15 * nbasis
    if ws.aves === nothing || ws.ns != ns || ws.nbasis != nbasis
        ws.aves = new_aves(inputs, nbasis)
        ws.ns = ns
        ws.nbasis = nbasis
    end
    if ws.iur != iur || size(ws.amat, 1) != iur
        K = Complex{T}
        ws.amat = Matrix{K}(undef, iur, iur)
        ws.bmat = Matrix{K}(undef, iur, iur)
        ws.amat_copy = Matrix{K}(undef, iur, iur)
        ws.bmat_copy = Matrix{K}(undef, iur, iur)
        ws.zmat = Matrix{K}(undef, iur, iur)
        ws.eigenvector = Vector{K}(undef, iur)
        ws.iur = iur
    end
    return ws
end

function tjlf_LS(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T},
            ky::T,
            nbasis::Int,
            vexb_shear_s::T,
            ky_index::Int;
            kx0_e = T(NaN),
            gamma_reference_kx0::Vector{T} = T[],
            freq_reference_kx0::Vector{T} = T[],
            outputGeo::Union{OutputGeometry{T},Missing} = missing,
            find_eigenvector::Bool = false,
            workspace::Union{Nothing,LSWorkspace{T}} = nothing,
            aves = workspace === nothing ? new_aves(inputs, nbasis) : fit_workspace!(workspace, inputs, nbasis).aves::AvesTuple{T},
            use_gpu::Bool = false) where T <: Real

    epsilon1 = 1.0e-12
    nmodes_in = inputs.NMODES
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    # compute the rank of the eigenmatrix iur
    nroot = 15
    iur = (ns - ns0 + 1) * nroot * nbasis

    alpha_quench_in = inputs.ALPHA_QUENCH
    new_eikonal_in = inputs.NEW_EIKONAL
    filter_in = inputs.FILTER

    R_unit = satParams.R_unit
    q_unit = satParams.q_unit

    new_geometry = false ######################## hardcoded for now ########################
    new_width = true ######################## hardcoded for now ########################
    # check co-dependencies
    if (new_geometry)
        new_width = true
    end
    if (new_width)
        new_matrix = true
    end

    if (new_eikonal_in)
        if (new_geometry)
            error("REMEMBER TO UPDATE R and Q_UNIT")
            # MILLER geometry
            igeo = 1
            if (igeo == 1)
                # miller_geo(inputs)
                # mercier_luc(inputs)
            # FOURIER geometry
            elseif(igeo==2)
                error("sorry fourier geometry not implemented :(")
                fourier_geo(inputs)
            # ELITE geometry
            elseif(igeo==3)
                error("sorry fourier geometry not implemented :(")
                ELITE_geo(inputs)
            end
            # compute the eikonal functions for general geometry (igeo>0)
            # if(igeo > 0) mercier_luc(inputs) end
        end
        new_geometry = false

        #  load the x-grid eikonal functions v_QL_out,b0x
        if(ismissing(outputGeo))
            outputGeo = xgrid_functions_geo(inputs, satParams, outputHermite, ky, ky_index; kx0_e)
        end
    end  #new_eikonal_in

    new_matrix = true ######################## hardcode for now ########################
    if (new_matrix)
        ave, aveH, aveWH, aveKH,
        aveG, aveWG, aveKG,
        aveGrad, aveGradB = get_matrix(inputs, outputGeo, outputHermite, ky, nbasis, ky_index; aves)
    end

    K = Complex{T}
    if workspace === nothing
        amat = Matrix{K}(undef, iur, iur)
        bmat = Matrix{K}(undef, iur, iur)
    else
        # reuse the caller's buffers (sized by fit_workspace! in the kwarg default above);
        # tjlf_eigensolver overwrites every entry of both matrices
        fit_workspace!(workspace, inputs, nbasis)
        amat = workspace.amat
        bmat = workspace.bmat
    end
    #  solver for linear eigenmodes of tglf equations
    eigenvalues, v = tjlf_eigensolver(inputs,outputGeo,satParams,ave,aveH,aveWH,aveKH,aveG,aveWG,aveKG,aveGrad,aveGradB, nbasis,ky, amat,bmat,ky_index,find_eigenvector; use_gpu=use_gpu, workspace=workspace)

    rr = real.(eigenvalues)
    ri = imag.(eigenvalues)

    # filter out numerical instabilities that sometimes occur with high mode frequency
    if filter_in > 0.0
        max_freq = 2 * abs(ave.wdh[1, 1]) / R_unit
        for is = ns0:ns
            rlts = inputs.RLTS[is]
            rlns = inputs.RLNS[is]
            as = inputs.AS[is]
            zs = inputs.ZS[is]
            test = abs(as * zs * (aveH.hp3p0[is, 1, 1] * rlns + 1.5 * (aveH.hr13p0[is, 1, 1] - aveH.hp3p0[is, 1, 1]) * rlts))
            max_freq = max(max_freq, test)
        end
        max_freq *= filter_in * abs(ky)
        # if imaginary part > max freq, flip sign of real part
        rr .*= ifelse.((rr .> 0.0) .& (abs.(ri) .> max_freq), -1, 1)
    end

    jmax = zeros(Int, nmodes_in)
    if size(gamma_reference_kx0,1) == 0
        gamma_out = zeros(T, nmodes_in)
        freq_out = zeros(T, nmodes_in)
    else
        gamma_out = similar(gamma_reference_kx0)
        freq_out = similar(freq_reference_kx0)
    end

    if (inputs.IBRANCH == 0)
        di = zeros(Int, size(rr))
        de = zeros(Int, size(rr))
        nmodes_out = nmodes_in
        # sort the unstable modes into electron and ion frequencies
        mi = 0
        me = 0
        for j1 = eachindex(rr)
            if (rr[j1] > epsilon1)
                if (ri[j1] > 0.0)
                    # note that ri = -freq, rr = gamma
                    mi = mi + 1
                    di[mi] = j1
                else
                    me = me + 1
                    de[me] = j1
                end
            end
        end

        # find the most unstable mode for each branch
        if (me > 0)
            zgamax = 0.0
            for iroot = 1:me
                if (rr[de[iroot]] > zgamax)
                    zgamax = rr[de[iroot]]
                    jmax[1] = de[iroot]
                end
            end
            gamma_out[1] = rr[jmax[1]]
            freq_out[1] = -ri[jmax[1]]
        end
        if (mi > 0)
            zgamax = 0.0
            for iroot = 1:mi
                if (rr[di[iroot]] > zgamax)
                    zgamax = rr[di[iroot]]
                    jmax[2] = di[iroot]
                end
            end
            gamma_out[2] = rr[jmax[2]]
            freq_out[2] = -ri[jmax[2]]
        end

    elseif(inputs.IBRANCH==-1)
        # find the top nmodes most unstable modes
        ### put the unstable modes in ascending order by growthrate
        jmax = sortperm(rr) #### sort_eigenvalues(nmodes_in,jmax), i believe this does the same thing
        jmax .= ifelse.(rr[jmax] .> epsilon1, jmax, 0)
        reverse!(jmax)
        nmodes_out = 0
        for j1 = 1:min(nmodes_in, size(rr, 1))
            if (jmax[j1] != 0)
                nmodes_out = nmodes_out + 1
                gamma_out[j1] = rr[jmax[j1]]
                freq_out[j1] = -ri[jmax[j1]]
            end
        end
    end

    # apply quench rule
    if (alpha_quench_in != 0.0)
        for j1 = 1:nmodes_in
            gamma_out[j1] = get_gamma_net(inputs, vexb_shear_s, gamma_out[j1])
        end
    # use spectral shift model for second pass
    elseif (vexb_shear_s != 0.0)
        gamma_out .= gamma_reference_kx0
        freq_out .= freq_reference_kx0
    end

    phi_bar_out = fill(NaN, nmodes_out)

    # get the fluxes for the most unstable modes
    if (inputs.IFLUX || find_eigenvector)
        #  initalize output to zero
        phi_bar_out = zeros(T, nmodes_out)
        a_par_bar_out = zeros(T, nmodes_out)
        b_par_bar_out = zeros(T, nmodes_out)
        v_bar_out = zeros(T, nmodes_out)
        ne_te_phase_out = zeros(T, nmodes_out)

        field_weight_out = zeros(K, 3, nbasis, nmodes_out)

        particle_QL_out = zeros(T, 3, ns, nmodes_out)
        energy_QL_out = zeros(T, 3, ns, nmodes_out)
        stress_par_QL_out = zeros(T, 3, ns, nmodes_out)
        stress_tor_QL_out = zeros(T, 3, ns, nmodes_out)
        exchange_QL_out = zeros(T, 3, ns, nmodes_out)

        N_QL_out = zeros(T, ns, nmodes_out)
        T_QL_out = zeros(T, ns, nmodes_out)
        U_QL_out = zeros(T, ns, nmodes_out)
        Q_QL_out = zeros(T, ns, nmodes_out)
        N_bar_out = zeros(T, ns, nmodes_out)
        T_bar_out = zeros(T, ns, nmodes_out)
        U_bar_out = zeros(T, ns, nmodes_out)
        Q_bar_out = zeros(T, ns, nmodes_out)
        Ns_Ts_phase_out = zeros(T, ns, nmodes_out)

        wd_bar_out = zeros(T, nmodes_out)
        b0_bar_out = zeros(T, nmodes_out)
        modB_bar_out = zeros(T, nmodes_out)
        v_QL_out = zeros(T, nmodes_out)
        a_par_QL_out = zeros(T, nmodes_out)
        b_par_QL_out = zeros(T, nmodes_out)
        kx_bar_out = zeros(T, nmodes_out)
        kpar_bar_out = zeros(T, nmodes_out)

        # used for computing eigenvector
        if inputs.SMALL != 0.0
            if workspace === nothing
                zmat = similar(amat)
                eigenvector = Vector{K}(undef,iur)
            else
                zmat = workspace.zmat
                eigenvector = workspace.eigenvector
            end
        end
        for imax = 1:nmodes_out
            if (jmax[imax] > 0)
                if inputs.SMALL != 0.0
                    # calculate eigenvector
                    eigenvector .= inputs.SMALL
                    # alpha and beta from eigensolver.jl
                    @. zmat = amat - (inputs.SMALL + rr[jmax[imax]] + im * ri[jmax[imax]]) * bmat

                    if use_gpu && K === ComplexF64
                        # Inverse-iteration LU solve on the GPU (CUSOLVER getrf/getrs),
                        # mirroring the eigenvalue solve. Moves this O(n^3) factorisation
                        # off the CPU; overwrites zmat.
                        eigenvector = TJLF._gpu_lu_solve!(zmat, eigenvector)
                    else
                        # in-place inverse iteration: lu!(zmat) factors in place (no
                        # matrix copy), ldiv! solves in place. zmat is rebuilt fresh
                        # each iteration above, so destroying it here is safe. This
                        # avoids the ~2.24 MiB/solve allocation of `zmat \ rhs`, which
                        # under multithreading triggered heavy GC contention.
                        ldiv!(lu!(zmat), eigenvector)
                    end
                else
                    if inputs.FIND_EIGEN || isnan(v[1,1])
                        nev1 = size(amat)[1]
                        L = construct_linear_map(sparse(amat), sparse(bmat), eigenvalues[jmax[imax]])
                        _, vec, _ = KrylovKit.eigsolve(L, nev1, 1, :LM)
                        eigenvector = vec[1]
                    else
                        eigenvector = v[:, jmax[imax]]
                    end
                end

                Ns_Ts_phase,
                Ne_Te_phase,
                N_weight,
                T_weight,
                U_weight,
                Q_weight,
                wd_bar,
                b0_bar,
                modB_bar,
                v_weight,
                a_par_weight,
                b_par_weight,
                kx_bar,
                kpar_bar,
                field_weight_QL_out,
                particle_weight,
                energy_weight,
                stress_par_weight,
                stress_tor_weight,
                exchange_weight = get_QL_weights(inputs, ave, aveH, ky, nbasis, eigenvalues[jmax[imax]], eigenvector)

                #### probably outputs
                wd_bar_out[imax] = wd_bar
                b0_bar_out[imax] = b0_bar
                modB_bar_out[imax] = modB_bar
                v_QL_out[imax] = v_weight
                a_par_QL_out[imax] = a_par_weight
                b_par_QL_out[imax] = b_par_weight
                kx_bar_out[imax] = kx_bar
                kpar_bar_out[imax] = kpar_bar/(R_unit*q_unit*inputs.WIDTH_SPECTRUM[ky_index])

                field_weight_out[:, :, imax] .= field_weight_QL_out

                particle_QL_out[:, :, imax] .= particle_weight
                energy_QL_out[:, :, imax] .= energy_weight
                stress_par_QL_out[:, :, imax] .= stress_par_weight
                stress_tor_QL_out[:, :, imax] .= stress_tor_weight
                exchange_QL_out[:, :, imax] .= exchange_weight

                N_QL_out[:, imax] .= N_weight
                T_QL_out[:, imax] .= T_weight
                U_QL_out[:, imax] .= U_weight
                Q_QL_out[:, imax] .= Q_weight
                Ns_Ts_phase_out[:, imax] .= Ns_Ts_phase

                ne_te_phase_out[imax] = Ne_Te_phase

                if (abs(v_QL_out[imax]) < epsilon1)
                    v_bar_out[imax] = 0.0
                    phi2_bar = 0.0
                else
                    kyi = ky
                    v_bar_out[imax] = get_intensity(inputs, ave, satParams, outputGeo.kx0_e, R_unit, kyi, gamma_out[imax]) ############### can use some cleaning
                    phi2_bar = v_bar_out[imax] / v_QL_out[imax]
                end

                phi_bar_out[imax] = phi2_bar

            end
        end
        a_par_bar_out .= phi_bar_out .* a_par_QL_out
        b_par_bar_out .= phi_bar_out .* b_par_QL_out

        phi_bar_out_diagonal = Diagonal(phi_bar_out)

        N_bar_out[ns0:ns, :] = N_QL_out[ns0:ns, :] * phi_bar_out_diagonal
        T_bar_out[ns0:ns, :] = T_QL_out[ns0:ns, :] * phi_bar_out_diagonal
        U_bar_out[ns0:ns, :] = U_QL_out[ns0:ns, :] * phi_bar_out_diagonal
        Q_bar_out[ns0:ns, :] = Q_QL_out[ns0:ns, :] * phi_bar_out_diagonal

        # check for inward ballooing
        sum_modB_bar = sum(v_bar_out .* modB_bar_out)
        sum_v_bar = sum(v_bar_out)

        ft_test = 0.0
        if (sum_v_bar > epsilon1)
            ft_test = sum_modB_bar / sum_v_bar
        end
        modB_min = abs(minimum(satParams.B_geo))
        ft_test = ft_test / modB_min

        reduce = 1.0
        if false # && inputs.SAT_RULE==0 && inputs.NBASIS_MAX!=inputs.NBASIS_MIN && gamma_nb_min_out<gamma_out[1] && gamma_nb_min_out<gamma_cutoff
            error("This is in the Fortran, but gamma_nb_min_out should always equal gamma_out[1] so the condition should never be true")
        end
        phi_bar_out .*= reduce # this does nothing for the reason above ^^^

        phi_bar_output = reshape(phi_bar_out,(1,1,nmodes_out))

        return nmodes_out, gamma_out, freq_out,
        particle_QL_out, energy_QL_out, stress_tor_QL_out, stress_par_QL_out, exchange_QL_out,
        ft_test, field_weight_out, phi_bar_output
    end

    phi_bar_output= reshape(phi_bar_out,(1,1,nmodes_out))
    particle_QL_out = fill(NaN, (3, ns, nmodes_in))
    energy_QL_out = fill(NaN, (3, ns, nmodes_in))
    stress_par_QL_out = fill(NaN, (3, ns, nmodes_in))
    stress_tor_QL_out = fill(NaN, (3, ns, nmodes_in))
    exchange_QL_out = fill(NaN, (3, ns, nmodes_in))
    field_weight_out = fill(NaN, (3, nbasis, nmodes_in))

    return nmodes_out, gamma_out, freq_out,
    particle_QL_out, energy_QL_out, stress_tor_QL_out, stress_par_QL_out, exchange_QL_out,
    NaN, field_weight_out, phi_bar_output
    # return  gamma_out,
    #         freq_out,
    #         v_QL_out,
    #         a_par_QL_out,
    #         b_par_QL_out,
    #         phi_bar_out,
    #         a_par_bar_out,
    #         b_par_bar_out,
    #         v_bar_out,
    #         ne_te_phase_out,
    #         field_weight_out,
    #         particle_QL_out,
    #         energy_QL_out,
    #         stress_par_QL_out,
    #         stress_tor_QL_out,
    #         exchange_QL_out,
    #         N_QL_out,
    #         T_QL_out,
    #         U_QL_out,
    #         Q_QL_out,
    #         N_bar_out,
    #         T_bar_out,
    #         U_bar_out,
    #         Q_bar_out,
    #         Ns_Ts_phase_out

end

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

"""
    get_wavefunction(inputs::InputTJLF, satParams::SaturationParameters, field_weight_out::Array{ComplexF64,3})

Reconstruct the linear mode structure (wavefunction) along the extended
ballooning angle from the Hermite-basis field weights produced by [`run`](@ref)
(its `field_weight_out`).

parameters:
    inputs::InputTJLF                       - the input used for the run
    satParams::SaturationParameters         - saturation/geometry params (`get_sat_params`)
    field_weight_out::Array{ComplexF64,3}   - mode field weights [field, basis, mode] from `run`

outputs:
    plot_field_out   - fields on the angle grid   [mode, field, point]
    plot_angle_out   - ballooning angle θ per grid point
    max_plot         - number of angle-grid points
    nmodes_out       - number of modes returned

description:
    Sums the basis-`j` Hermite functions weighted by `field_weight_out` onto a
    fixed extended-angle grid, for each mode and each active field
    (electrostatic, plus `USE_BPER`/`USE_BPAR` magnetic components). Useful for
    inspecting mode localization and parity.
"""
function get_wavefunction(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, field_weight_out::Array{ComplexF64, 3}) where T<:Real 
    ms = 128
    npi = 9
    np = Int(ms/8)
    nb = inputs.NBASIS_MAX
    igeo = 1 # Hard-Coded for now as with the previous LS functions:
    max_plot = Int(18*ms/8+1) # 289 length vector
    maxmodes = inputs.NMODES
    xp = fill(NaN, max_plot)
    plot_angle_out = fill(NaN, max_plot)
    plot_field_out = zeros(ComplexF64, maxmodes, 3, max_plot)
    hp = fill(NaN, (nb, max_plot))
    nbasis = inputs.NBASIS_MAX
    nmodes_out = inputs.NMODES
    if (igeo == 0)
        dx = npi*2.0*pi/(max_plot-1) # /(max_plot-1)
        for i = 1:max_plot
            xp[i] = -npi*pi+(i-1)*dx
            plot_angle_out[i] = xp[i]
        end
    else
        # satParams.y is the 1-based port of Fortran y(0:ms): y[m] = y(m-1), so
        # Fortran y(ms) is y[ms+1] here (y[ms+1] = Ly, the full field-line length).
        dx = 2.0*pi/(satParams.y[ms+1]*inputs.WIDTH) # y comes from s_grid. width_in is not a tjlf struct thing; they mean width.
        j0 = npi*np+1
        xp[j0] = 0.0 
        plot_angle_out[j0] = 0.0
        j = 0
        k = 0
        imax = np*npi
        for i = 1:imax
            j += 1
            if (j > 2*np)
                j -= 2*np
                k += 1
            end # I needed to change the ms to ms+1 as the max dimension of y is ms+1 (see tjlf_geometry.jl)
            # Fortran y(4*j)/t_s(4*j) port to 1-based y[4*j+1]/theta[4*j+1]. Using
            # y[4*j] shifted the positive-theta samples by one cell (the negative side
            # was exact), making the wavefunction spuriously asymmetric and inflating
            # the tearing-parity metric |phi(theta)-phi(-theta)| for marginal modes.
            xp[j0+i] = dx*(k*satParams.y[ms+1] + satParams.y[4*j+1])
            xp[j0-i] = -dx*((k+1)*satParams.y[ms+1] - satParams.y[ms+1-4*j])
            plot_angle_out[j0+i] = -(k*satParams.theta[ms+1] + satParams.theta[4*j+1]) #t_s comes from sgrid module. t_s is theta in satParams
            plot_angle_out[j0-i] = (k+1)*satParams.theta[ms+1] - satParams.theta[ms+1-4*j]
        end
    end
    hp0 = sqrt(2)/pi^0.25
    for i = 1:max_plot
        hp[1, i] = hp0*exp(-xp[i]*xp[i]/2.0)
        hp[2, i] = xp[i]*sqrt(2)*hp0*exp(-xp[i]*xp[i]/2.0)
        if (nbasis > 2) # This needs some revision most likely. 
            for j = 3:nbasis
                hp[j, i] = xp[i]*sqrt(2.0/(j-1))*hp[j-1, i] - sqrt((j-2)/(j-1))*hp[j-2, i]
            end # end of j
        end # end of if nbasis > 2
    end # end of i 
    # compute the fields of the theta-grid
    for n = 1:nmodes_out
        for k = 1:3
            for i = 1:max_plot
                if ((k == 1) || (k == 2 && inputs.USE_BPER) || (k == 3 && inputs.USE_BPAR))
                    for j = 1:nbasis
                        plot_field_out[n, k, i] = plot_field_out[n, k, i] + field_weight_out[k, j, n]*hp[j, i]
                    end
                end
            end # end of i
        end # end of k
    end # end of n
    return plot_field_out, plot_angle_out, max_plot, nmodes_out #nfields_out?
end

"""
    get_intensity(inputs::InputTJLF{T}, ave::Ave{T}, satParams::SaturationParameters{T}, kx0_e::T, R_unit::T, kp::T, gp::T) where T<:Real

description:
    helper function to get intensity coefficent given the saturation rule
"""
function get_intensity(inputs::InputTJLF{T}, ave::Ave{T}, satParams::SaturationParameters{T}, kx0_e::T, R_unit::T, kp::T, gp::T) where {T<:Real}

    nmodes_in = inputs.NMODES
    sat_rule_in = inputs.SAT_RULE
    etg_factor_in = inputs.ETG_FACTOR
    alpha_quench_in = inputs.ALPHA_QUENCH

    pols = (ave.p0[1, 1] / abs(inputs.AS[1] * inputs.ZS[1]^2))^2
    ks = kp * √(inputs.TAUS[1] * inputs.MASS[2]) / abs(inputs.ZS[1])
    measure = √(inputs.TAUS[1] * inputs.MASS[2])

    if (sat_rule_in == 0)
        igeo = 1 ################### MILLER GEO for now
        if (igeo == 0)
            if (nmodes_in <= 2) #this fit is for nmodes_in=2
                cnorm = 30.40 * pols
                exponent1 = 1.657
            else #this fit is for nmodes_in=4
                cnorm = 24.58 * pols
                exponent1 = 1.761
            end

            if (ks > 1.0)
                cnorm = cnorm / (ks)^etg_factor_in
            end
            c1 = 0.0

        elseif (igeo >= 1)
            if (nmodes_in <= 2) #this fit is for nmodes_in=2
                cnorm = 32.48 * pols
                exponent1 = 1.547
                c1 = 0.534
            else #this fit is for nmodes_in=4
                cnorm = 30.03 * pols
                exponent1 = 1.66
                c1 = 0.234
            end

            if (ks > 1.0)
                cnorm = cnorm / (ks)^etg_factor_in
            end
        end

        wd0 = ks * √(inputs.TAUS[1] / inputs.MASS[2]) / R_unit  #renomalized for scale invariance
        gnet = gp / wd0
        intensity = cnorm * (wd0^2) * (gnet^exponent1 + c1 * gnet) / (kp^4)

        if (alpha_quench_in == 0.0 && abs(kx0_e) > 0.0)
            intensity = intensity / (1.0 + 0.56 * kx0_e^2)^2
            intensity = intensity / (1.0 + (1.15 * kx0_e)^4)^2
        end

        # geometry factor from get_sat_params: 1.0 for GYRO units, 0.946/qrat_geo(0) otherwise
        intensity = intensity * satParams.SAT_geo0 * measure
    elseif (sat_rule_in >= 1)
        intensity = 1.0
    end

    return intensity
end

#--------------------------------------------------------------

function get_gamma_net(inputs::InputTJLF{T}, vexb_shear_s::T, gp::T) where {T<:Real}
    alpha_quench_in = inputs.ALPHA_QUENCH
    kappa_loc = inputs.KAPPA_LOC ##### only true for MILLER
    alpha_exb = 0.3
    igeo = 1 ####### MILLER GEOMETRY for now

    if (igeo == 1)
        alpha_exb = 0.3 * √(kappa_loc)
    end
    get_gamma_net = max(gp - abs(alpha_exb * alpha_quench_in * vexb_shear_s), 0.0)

    return get_gamma_net
end

#--------------------------------------------------------------

"""
    get_QL_weights(inputs::InputTJLF{T}, ave::Ave{T}, aveH::AveH{T}, ky::T, nbasis::Int, eigenvalue::K, v::Vector{K}) where T<:Real where K<:Complex

description:
    helper function to compute the quasilinear weights for a single eigenmode with eigenvector v.
    All of the QL weights are normalized to phi_norm
"""
function get_QL_weights(inputs::InputTJLF{T}, ave::Ave{T}, aveH::AveH{T},
    ky::T, nbasis::Int, eigenvalue::K, v::Vector{K}) where {T<:Real} where {K<:Complex}

    epsilon1 = 1.e-12
    sat_rule_in = inputs.SAT_RULE
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    nroot = 15 ### hardcoded
    iur = (ns - ns0 + 1) * nroot * nbasis

    taus = inputs.TAUS
    mass = inputs.MASS
    as = inputs.AS
    zs = inputs.ZS
    vs = .√(taus ./ mass)
    vpar = inputs.VPAR
    vpar_shear = inputs.VPAR_SHEAR

    vpar_model_in = inputs.VPAR_MODEL
    alpha_mach_in = inputs.ALPHA_MACH
    sign_It_in = inputs.SIGN_IT
    alpha_p_in = inputs.ALPHA_P
    freq_QL = im * eigenvalue

    n = zeros(K, nbasis, ns)
    u_par = zeros(K, nbasis, ns)
    p_par = zeros(K, nbasis, ns)
    p_tot = zeros(K, nbasis, ns)
    q_par = zeros(K, nbasis, ns)
    q_tot = zeros(K, nbasis, ns)

    for is = ns0:ns
        j = (is - ns0) * nroot * nbasis
        for i = 1:nbasis
            n[i, is] = v[j+i]
            u_par[i, is] = v[j+nbasis+i]
            p_par[i, is] = v[j+nbasis*2+i]
            p_tot[i, is] = v[j+nbasis*3+i]
            q_par[i, is] = v[j+nbasis*4+i]
            q_tot[i, is] = v[j+nbasis*5+i]
            if (nroot > 6)
                n[i, is] = n[i, is] - v[j+nbasis*6+i] + v[j+nbasis*12+i]
                u_par[i, is] = u_par[i, is] - v[j+nbasis*7+i]
                p_par[i, is] = p_par[i, is] - v[j+nbasis*8+i] + v[j+nbasis*13+i]
                p_tot[i, is] = p_tot[i, is] - v[j+nbasis*9+i] + v[j+nbasis*14+i]
                q_par[i, is] = q_par[i, is] - v[j+nbasis*10+i]
                q_tot[i, is] = q_tot[i, is] - v[j+nbasis*11+i]
            end
        end
    end


    # compute vnorm
    vnorm = 0.0
    if (ns <= 2)
        vnorm = real((adjoint(v) * v))
    else #weight vnorm by equililibrium densities
        j = 1
        for i = 1:iur
            if (i > j * nbasis * nroot)
                j = j + 1
            end
            vnorm = vnorm + real(v[i]*conj(v[i])) * abs(inputs.AS[j] * inputs.ZS[j])
        end
        vnorm = vnorm / abs(inputs.AS[1] * inputs.ZS[1])   #normalize to electron charge density
    end

    # compute the electromagnetic potentials
    betae_s = inputs.BETAE ##### not true for 'GENE' units
    betae_psi = 0.0
    if (inputs.USE_BPER)
        betae_psi = 0.5 * betae_s / (ky * ky)
    end
    betae_sig = 0.0
    if (inputs.USE_BPAR)
        betae_sig = 0.5 * betae_s
    end

    phi = zeros(K, nbasis)
    psi = zeros(K, nbasis)
    bsig = zeros(K, nbasis)
    U0 = sum((alpha_mach_in * sign_It_in) .* vpar .* zs .^ 2 .* as ./ taus) ### defined in startup.f90
    @views phi .= sum(ave.p0inv * (n[:, ns0:ns] * Diagonal((as.*zs)[ns0:ns])), dims=2)
    if (inputs.USE_BPER)
        psi .= (betae_psi .* sum(ave.b0inv * (u_par[:, ns0:ns] * Diagonal((as.*zs.*vs)[ns0:ns])), dims=2))
        if (vpar_model_in == 0)
            @views phi .= phi .+ (U0 * betae_psi) .* sum(ave.bpinv * (u_par[:, ns0:ns] * Diagonal((as.*zs.*vs)[ns0:ns])), dims=2)
            @views psi .= psi .- (U0 * betae_psi) .* sum(ave.bpinv * (n[:, ns0:ns] * Diagonal((as.*zs)[ns0:ns])), dims=2)
        end
    end
    if (inputs.USE_BPAR)
        @views bsig .= -betae_sig .* (1.5 .* p_tot[:, ns0:ns] .- 0.5 .* p_par[:, ns0:ns]) * (as.*taus)[ns0:ns]
    end

    # add the adiabatic terms to the total moments
    @views n[:, ns0:ns] .= n[:, ns0:ns] .- (phi .* transpose((zs./taus)[ns0:ns])) #### outer product to make matrix, idk why its this order tbh -DSUN
    @views p_par[:, ns0:ns] .= p_par[:, ns0:ns] .- (phi .* transpose((zs./taus)[ns0:ns]))
    @views p_tot[:, ns0:ns] .= p_tot[:, ns0:ns] .- (phi .* transpose((zs./taus)[ns0:ns]))

    # compute phi_norm, psi_norm, sig_norm
    phi_norm = real(adjoint(phi) * phi)
    psi_norm = real(adjoint(psi) * psi)
    bsig_norm = real(adjoint(bsig) * bsig)
    if (phi_norm < epsilon1)
        phi_norm = epsilon1
    end

    #save the field weights
    field_weight_QL_out = Matrix{K}(undef, 3, nbasis)
    field_weight_QL_out[1, :] .= phi .* (im / √(phi_norm))
    field_weight_QL_out[2, :] .= psi .* (im / √(phi_norm))
    field_weight_QL_out[3, :] .= bsig .* (im / √(phi_norm))

    #compute <phi|*|phi> averages
    phi_wd_phi = quadform(phi, ave.wdh)
    phi_b0_phi = quadform(phi, ave.b0)
    phi_modB_phi = quadform(phi, ave.c_par_par)
    phi_kx_phi = quadform(phi, ave.kx)
    phi_kpar_phi = im * quadform(phi, ave.kpar)

    wd_bar = real(phi_wd_phi) / phi_norm
    b0_bar = real(phi_b0_phi) / phi_norm
    modB_bar = abs(real(phi_modB_phi) / phi_norm)
    kx_bar = real(phi_kx_phi) / phi_norm
    kpar_bar = real(phi_kpar_phi) / phi_norm


    # fill the stress moments
    stress_par = zeros(K, nbasis, ns, 2)
    stress_per = zeros(K, nbasis, ns, 2)

    stress_correction = fill(one(T), ns - ns0 + 1)
    if (sat_rule_in == 0)
        @views wp = (ky * abs(alpha_p_in)) .* aveH.hp1[ns0:ns, 1, 1] .* vpar_shear ./ vs
        stress_correction .= (imag(freq_QL) .+ 2.0 .* wp) ./ (imag(freq_QL) .+ wp)
    end
    stress_correction = Diagonal(stress_correction)

    @views stress_par[:, ns0:ns, 1] .= u_par[:, ns0:ns] * stress_correction
    @views stress_par[:, ns0:ns, 2] .= p_par[:, ns0:ns] * stress_correction
    @views stress_per[:, ns0:ns, 1] .= ((im * ky) .* ave.kx) * (1.5 .* p_tot[:, ns0:ns] .- 0.5 .* p_par[:, ns0:ns])   # (is,j) x (j,i)
    @views stress_per[:, ns0:ns, 2] .= ((im * ky) .* ave.kx) * (1.5 .* q_tot[:, ns0:ns] .- 0.5 .* q_par[:, ns0:ns])


    # compute the quasilinear weights for the fluxes
    particle_weight = zeros(T, 3, ns)
    energy_weight = zeros(T, 3, ns)
    stress_par_weight = zeros(T, 3, ns)
    stress_tor_weight = zeros(T, 3, ns)
    exchange_weight = zeros(T, 3, ns)

    ### real() with an im in it is funky
    #### CHECK THIS PLEASE
    @views particle_weight[1, ns0:ns] .= .-vec(imag.(adjoint(phi) * n[:, ns0:ns]))
    @views energy_weight[1, ns0:ns] .= .-vec(imag.(adjoint(phi) * p_tot[:, ns0:ns]))
    @views stress_par_weight[1, ns0:ns] .= .-vec(imag.(adjoint(phi) * (ave.c_par_par * stress_par[:, ns0:ns, 1])))
    @views stress_tor_weight[1, ns0:ns] .= .-vec(imag.(adjoint(phi) * (ave.c_tor_par * stress_par[:, ns0:ns, 1]
                                                                 .+ ave.c_tor_per * stress_per[:, ns0:ns, 1])))
    @views exchange_weight[1, ns0:ns] .= .-vec(imag.(freq_QL .* adjoint(phi) * (n[:, ns0:ns])) * Diagonal(zs[ns0:ns]))


    if (inputs.USE_BPER)
        @views particle_weight[2, ns0:ns] .= vec(imag.(adjoint(psi) * u_par[:, ns0:ns]) * Diagonal(vs[ns0:ns]))
        @views energy_weight[2, ns0:ns] .= vec(imag.(adjoint(psi) * q_tot[:, ns0:ns]) * Diagonal(vs[ns0:ns]))
        @views exchange_weight[2, ns0:ns] .= vec(imag.(freq_QL .* adjoint(psi) * u_par[:, ns0:ns]) * Diagonal((zs.*vs)[ns0:ns]))
        @views stress_par_weight[2, ns0:ns] .= vec(imag.(adjoint(psi) * (ave.c_par_par * stress_par[:, ns0:ns, 2])))
        @views stress_tor_weight[2, ns0:ns] .= vec(imag.(adjoint(psi) * (ave.c_tor_par * stress_par[:, ns0:ns, 2]
                                                                  .+ ave.c_tor_per * stress_per[:, ns0:ns, 2])))
    end
    if (inputs.USE_BPAR)
        @views particle_weight[3, ns0:ns] .= .-vec(imag.(adjoint(bsig) * (1.5 .* p_tot[:, ns0:ns] .- 0.5 .* p_par[:, ns0:ns]) * Diagonal((taus./zs)[ns0:ns])))
        @views exchange_weight[3, ns0:ns] .= .-vec(imag.(conj(freq_QL) * adjoint(bsig) * (1.5 .* p_tot[:, ns0:ns] .- 0.5 .* p_par[:, ns0:ns])) * Diagonal(taus[ns0:ns]))
    end

    @views particle_weight[:, ns0:ns] .= ky .* (particle_weight[:, ns0:ns] * Diagonal(as[ns0:ns])) ./ phi_norm
    @views energy_weight[:, ns0:ns] .= (1.5 * ky) .* (energy_weight[:, ns0:ns] * Diagonal((as.*taus)[ns0:ns])) ./ phi_norm
    @views stress_par_weight[:, ns0:ns] .= ky .* (stress_par_weight[:, ns0:ns] * Diagonal((as.*mass.*vs)[ns0:ns])) ./ phi_norm
    @views stress_tor_weight[:, ns0:ns] .= ky .* sign_It_in .* (stress_tor_weight[:, ns0:ns] * Diagonal((mass.*as.*vs)[ns0:ns])) ./ phi_norm
    @views exchange_weight[:, ns0:ns] .= (exchange_weight[:, ns0:ns] * Diagonal(as[ns0:ns])) ./ phi_norm


    #  add the vpar shifts to the total  moments
    if (vpar_model_in == 0)
        vpar_s = (alpha_mach_in * sign_It_in) .* vpar
        @views n[:, ns0:ns] = psi .* transpose((vpar_s.*zs./taus)[ns0:ns]) ### outer product is slightly weird
        @views u_par[:, ns0:ns] = phi .* transpose((.-(vpar_s ./ vs).*(zs./taus))[ns0:ns])
        @views p_par[:, ns0:ns] = psi .* transpose((vpar_s.*(zs./taus))[ns0:ns])
        @views p_tot[:, ns0:ns] = psi .* transpose((vpar_s.*(zs./taus))[ns0:ns])
        @views q_par[:, ns0:ns] = phi .* transpose((-3 .* (vpar_s./vs).*(zs./taus))[ns0:ns])
        @views q_tot[:, ns0:ns] = phi .* transpose((-(5 / 3).*(vpar_s./vs).*(zs./taus))[ns0:ns])
    end



    #### outputs
    # compute the density and temperature amplitude weights
    N_weight = zeros(T, ns)
    T_weight = zeros(T, ns)
    U_weight = zeros(T, ns)
    Q_weight = zeros(T, ns)
    temp = Matrix{K}(undef, nbasis, ns)

    temp .= p_tot .- n
    N_weight .= vec(sum(abs2.(n), dims=1)) ./ phi_norm
    T_weight .= vec(sum(abs2.(temp), dims=1)) ./ phi_norm
    U_weight .= vec(sum(abs2.(u_par), dims=1)) ./ phi_norm
    Q_weight .= vec(sum(abs2.(q_tot), dims=1)) ./ phi_norm
    v_weight = vnorm / phi_norm
    a_par_weight = psi_norm / phi_norm
    b_par_weight = bsig_norm / phi_norm


    #compute electron density-temperature phase
    @views Ne_Te_cos = real(adjoint(n[:, 1]) * temp[:, 1])
    @views Ne_Te_sin = imag(adjoint(n[:, 1]) * temp[:, 1])
    Ne_Te_phase = atan(Ne_Te_sin, Ne_Te_cos)

    #compute species density-temperature phase
    Ns_Ts_phase = zeros(T, ns)
    Ns_Ts_cos = zeros(T, ns)
    Ns_Ts_sin = zeros(T, ns)

    temp2 = @views @. conj(n[:, ns0:ns]) .* temp[:, ns0:ns]
    @views Ns_Ts_cos .= vec(sum(real(temp2), dims=1))
    @views Ns_Ts_sin .= vec(sum(imag(temp2), dims=1))
    @. Ns_Ts_phase = atan(Ns_Ts_sin, Ns_Ts_cos)

    return Ns_Ts_phase, Ne_Te_phase,
    N_weight, T_weight, U_weight, Q_weight,
    wd_bar, b0_bar, modB_bar, v_weight, a_par_weight, b_par_weight, kx_bar, kpar_bar,
    field_weight_QL_out, particle_weight, energy_weight, stress_par_weight, stress_tor_weight, exchange_weight

end

function quadform(phi::Vector{<:Complex}, A::Matrix{<:Real})
    sum(adjoint(phi[i]) * A[i, j] * phi[j] for j in eachindex(phi), i in eachindex(phi))
end
