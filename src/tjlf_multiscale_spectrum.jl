"""
    get_rho_ion(inputs::InputTJLF{T}) where T<:Real

description:
    single definition of the ion gyroradius used by the saturation rules,
    mirroring the module-level rho_ion of the Fortran (tglf_startup.f90):
    charge-weighted average over ions when USE_AVE_ION_GRID is set,
    otherwise the plain species-2 value
"""
function get_rho_ion(inputs::InputTJLF{T}) where T<:Real
    if !inputs.USE_AVE_ION_GRID
        return √(inputs.TAUS[2]*inputs.MASS[2]) / abs(inputs.ZS[2])
    end
    rho_ion = 0.0
    charge = 0.0
    for is in 2:inputs.NS
        if inputs.ZS[is]*inputs.AS[is]/abs(inputs.ZS[1]*inputs.AS[1]) > 0.1
            charge += inputs.ZS[is]*inputs.AS[is]
            rho_ion += inputs.AS[is]*√(inputs.TAUS[is]*inputs.MASS[is])
        end
    end
    if charge!=0.0 rho_ion = rho_ion/charge end
    return rho_ion
end

"""
    get_zonal_mixing(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, most_unstable_gamma::AbstractArray{T}) where T<:Real

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    satParams::SaturationParameters{T}  - SaturationParameters struct constructed in tjlf_geometry.jl
    most_unstable_gamma::Vector{T}      - vector of most unstable gammas (mode = 1) along the ky spectrum

outputs:
    vzf_mix - zonal flow mixing rate
    kymax_mix - ky value at the calculated vzf_mix
    jmax_mix - index of ky_spect array where vzf_mix is calculated

description:
    finds the maximum of gamma/ky spectrum at low-k values by going through the ky_spect
    array, then after finding the max, interpolate the value to improve accuracy
"""
function get_zonal_mixing(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, most_unstable_gamma::Vector{T}) where T<:Real

    sat_rule_in = inputs.SAT_RULE
    ky_spect = inputs.KY_SPECTRUM
    alpha_zf = inputs.ALPHA_ZF

    rho_ion = get_rho_ion(inputs)

    # find the local maximum of most_unstable_gamma/ky_spect with the largest most_unstable_gamma/ky_spect^2
    kycut = 0.8/rho_ion
    kymin = 0.0

    testmax = 0.0
    jmax_mix = 1
    use_kymin = false

    if(alpha_zf < 0.0) 
        use_kymin = true
        kymin = 0.173 * √(2.0) / rho_ion 
    end
    # saturation rules
    if sat_rule_in==2 || sat_rule_in==3
        grad_r0 = satParams.grad_r0
        kycut = grad_r0 * kycut
        kymin = grad_r0 * kymin
    end

    # find the low and high ky peaks of gamma/ky

    # go through ky_spect except for last value
    # if next val is >= kymin and curr val is less than kycut save it
    # update testmax and max index if necessary

    j1 = nothing
    jmin = 0
    for j in 1:length(ky_spect)-1
        if ky_spect[j] < kymin
            jmin = j
        end
        if((ky_spect[j+1] >= kymin) && (ky_spect[j] <= kycut))
            # save index in case no max
            j1=j
            testmax1 = most_unstable_gamma[j]/ky_spect[j]
            # find maxes
            if(testmax1 > testmax)
                testmax = testmax1
                jmax_mix = j
            end
        end
    end

    if isnothing(j1)
        error("Check that ky spectrum has entries in the range kymin <= ky <= kycut.")
    end
    # if testmax is not updated a single time
    if(testmax==0.0) jmax_mix=j1 end

    # no unstable modes in range set kymax index to end of range
    # this is cut of at j1 since a maximum may not exist in the low-k range
    kymax1 = ky_spect[jmax_mix]
    gammamax1 = most_unstable_gamma[jmax_mix]




    # linearly interpolate to find a more accurate low-k maximum gamma/ky
    if(kymax1<kymin)
        kymax1 = kymin
        gammamax1 = (most_unstable_gamma[1]
            + (most_unstable_gamma[2]-most_unstable_gamma[1])*(kymin-ky_spect[1])/(ky_spect[2]-ky_spect[1]))
    end
    # determine kymax1 and gammamax1 bounded by the tree points f0,f1,f2
    # use a quadratic fit: f = a + b x + c x^2    to f = gamma/ky centered at jmax1
    # scale it to be quadratic where x goes from 0 to 1
    if(jmax_mix>1 && jmax_mix < j1)
        jmax1 = jmax_mix
        f0 = most_unstable_gamma[jmax1-1] / ky_spect[jmax1-1]
        f1 = most_unstable_gamma[jmax1] / ky_spect[jmax1]
        f2 = most_unstable_gamma[jmax1+1] / ky_spect[jmax1+1]
        deltaky = ky_spect[jmax1+1] - ky_spect[jmax1-1]
        x1 = (ky_spect[jmax1] - ky_spect[jmax1-1])/deltaky
        a = f0
        b = (f1 - f0*(1-x1*x1)-f2*x1*x1)/(x1-x1*x1)
        c = f2 - f0 - b

        # if f0>=f1 then f1 is not a local maximum (Fortran tie-break: .ge.)
        if(f0 >= f1)
            kymax1 = ky_spect[jmax1-1]
            gammamax1 = f0*kymax1

            #interpolate to find the value of gammamax1 at kymin
            if(kymax1 < kymin)
                kymax1 = kymin
                xmin = (kymin - ky_spect[jmax1-1])/deltaky
                gammamax1 = (a + b*xmin + c*xmin*xmin)*kymin
            end
        end

        # if f0<f1 then f1>f2 due to the maximum search
        # use the quadratic fit to refine the local maximum:
        if(f0<f1)
            xmax = -b/(2.0*c)
            xmin = 0.0

            if(ky_spect[jmax1-1]<kymin)
                xmin = (kymin - ky_spect[jmax1-1])/deltaky
            end

            # if xmax >= 1    use f2 as the maximum
            if(xmax >= 1.0)
                kymax1 = ky_spect[jmax1+1]
                gammamax1 = f2*kymax1
            elseif(xmax < xmin)
                # use the quadratic fit to determine gammamax1 at kymin
                if(xmin > 0.0)
                    kymax1 = kymin
                    gammamax1 = (a + b*xmin + c*xmin*xmin)*kymin
                # if xmax<=0 use f0 as the maximum
                else
                    kymax1 = ky_spect[jmax1-1]
                    gammamax1 = f0*kymax1
                end

            # the conditions f0<f1<f2 and xmin<xmax<1 are satisfied
            # use the quadratic fit to determine gammamax1 and kymax1
            else
                kymax1 = ky_spect[jmax1-1] + deltaky*xmax
                gammamax1 = (a + b*xmax + c*xmax*xmax)*kymax1
            end #xmax tests

        end #f0 < f1

    end    # jmax_mix > 1
    
    # Match Fortran behavior: if gammamax1 is zero, return vzf_mix = 0.0
    if gammamax1 == 0.0
        vzf_mix = 0.0
    elseif testmax == 0.0
        vzf_mix = 0.0
    else
        vzf_mix = gammamax1/kymax1
    end
    kymax_mix = kymax1
    ### commented out in original Fortran code
    # jmax_mix = jmax1

    return vzf_mix, kymax_mix, jmax_mix

end

"""
    mode_transition_function(x, y1, y2, x_ITG, x_TEM)

description:
    helper function that returns either the smaller value if x is less than, larger value if x is greater than, 
    or a linear interpolation of the two if x is within the bounds
"""
function mode_transition_function(x, y1, y2, x_ITG, x_TEM)

    if (x < x_ITG)
        y = y1
    elseif (x > x_TEM)
        y = y2
    else
        y = y1*((x_TEM - x) / (x_TEM - x_ITG)) + y2*((x - x_ITG) / (x_TEM - x_ITG))
    end

    return y
end

"""
    linear_interpolation(x::AbstractVector, y::AbstractVector, x0)

description:
    helper function that returns either the linear interpolated y value
    at a x0 value given the x and y arrays
"""
function linear_interpolation(x::AbstractVector, y::AbstractVector, x0)
    # strict comparison as in the Fortran: at x0 == x[i] the pair is (i-1, i), the same pair
    # that the callers test or fill before interpolating
    i = 1
    while (x[i] < x0)
        i += 1
    end
    return ((y[i] - y[i-1]) * x0 + (x[i] * y[i-1] - x[i-1]*y[i])) / (x[i] - x[i-1])

end

"""
    intensity_sat(inputs::InputTJLF{T},satParams::SaturationParameters{T},gamma_matrix::Array{T},QL_weights::Array{T,5},expsub::T=2.0,return_phi_params::Bool=false) where T<:Real
    
parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    satParams::SaturationParameters{T}  - SaturationParameters struct constructed in tjlf_geometry.jl
    gamma_matrix::Matrix{T}             - matrix of gamma (mode, ky)
    QL_weights::Array{T,5}              - array of quasilinear weights (field, species, mode, ky, type)
                                          type: (particle, energy, torodial stress, parallel stress, exchange)
    expsub::T (opt)                     - float for the exponent
    return_phi_params::T (opt)          - boolean flag to change how much is output

outputs:
  If return_phi_params=false (default): tuple (phinorm, QLA_P, QLA_E, QLA_O)
    phinorm           - intensity of the saturation rule (nky, nmodes)
    QLA_P             - QLA particle value (nmodes)
    QLA_E             - QLA energy value (nmodes)
    QLA_O             - QLA stress / exchange value (nmodes)

  If return_phi_params=true: NamedTuple with fields needed to reconstruct
  |phi|^2(kx, ky, mode) from the Staebler SAT kx-Lorentzian shape
    phinorm   [nky, nmodes] - SAT intensity at the kx peak (kx = kx0_e)
    kx_width  [nky]         - Lorentzian kx width (kx rms width)
    gammaeff  [nky, nmodes] - effective growth rate
    kx0_e     [nky]         - spectral shift (peak of kx distribution)
    ax, ay, exp_ax          - SAT Lorentzian coefficients
    (SAT2 adds d1, d2, kycut, b3 geometry parameters)

description:
    takes in the input.tglf file, ky values, gammas, quasilinear weights, and
    some other optional parameters and returns the intensity value and QLA
    parameter values dependent on the 3 saturation rules
    called get_multiscale_spectrum in the Fortran
"""
function intensity_sat(
    inputs::InputTJLF{T},
    satParams::SaturationParameters{T},
    gamma_matrix::Matrix{T},
    QL_weights::Array{T,5},
    expsub=T(2.0),
    return_phi_params::Bool=false;
    vzf_out_param=T(NaN),
    kymax_out_param=T(NaN),
    jmax_out_param::Int=0) where T<:Real

    ############ figure out how to make this prettier
    units_in = inputs.UNITS
    if inputs.SAT_RULE == 2 || inputs.SAT_RULE == 3
        inputs.UNITS = "CGYRO"
        units_in = "CGYRO"
    end

    sat_rule_in = inputs.SAT_RULE
    rlnp_cutoff = inputs.RLNP_CUTOFF
    beta_loc = inputs.BETA_LOC
    ns = inputs.NS
    rmaj_loc = inputs.RMAJ_LOC
    p_prime_loc = inputs.P_PRIME_LOC
    alpha_zf = inputs.ALPHA_ZF
    alpha_quench = inputs.ALPHA_QUENCH

    mass_2 = inputs.MASS[2]
    taus_1 = inputs.TAUS[1]
    nmodes = inputs.NMODES
    rho_ion = get_rho_ion(inputs)
    small = 10^-10
    ky_spect = inputs.KY_SPECTRUM
    nky = length(ky_spect)

    SAT_geo0_out = satParams.SAT_geo0
    SAT_geo1_out = satParams.SAT_geo1
    SAT_geo2_out = satParams.SAT_geo2
    Bt0_out = satParams.Bt0
    b_geo0_out = satParams.B_geo[1]
    grad_r0_out = satParams.grad_r0

    most_unstable_gamma = gamma_matrix[1, :] # SAT1 and SAT2 use the growth rates of the most unstable modes

    # Use passed parameters when available, otherwise calculate them
    if !isnan(vzf_out_param) && !isnan(kymax_out_param)
        vzf_out = vzf_out_param
        kymax_out = kymax_out_param
        if jmax_out_param > 0
            jmax_out = jmax_out_param
        else
            # Find jmax_out from kymax_out
            kymax_tolerance = 1e-10
            jmax_out = findfirst(x -> abs(x - kymax_out) < kymax_tolerance, ky_spect)
            if jmax_out === nothing
                jmax_out = argmin(abs.(ky_spect .- kymax_out))
            end
        end
    else
        vzf_out, kymax_out, jmax_out = get_zonal_mixing(inputs, satParams, most_unstable_gamma)
    end

    # Now calculate spectral shift with the correct parameters
    kx0_e = xgrid_functions_geo(inputs, satParams, gamma_matrix; vzf_out_param=vzf_out, kymax_out_param=kymax_out)

    # model fit parameters
    # Miller geometry values igeo=1
    if rlnp_cutoff > 0.0
        if beta_loc == 0.0
            dlnpdr = 0.0
            ptot = 0.0
            for i in 1:ns
                ptot = ptot + inputs.AS[i]*inputs.TAUS[i]
                dlnpdr = dlnpdr + inputs.AS[i]*inputs.TAUS[i]*(inputs.RLNS[i]+inputs.RLTS[i])
            end
            ### kwargs["RMAJ_LOC"] used for rmaj_input
            dlnpdr = rmaj_loc*dlnpdr/max(ptot,0.01)
        else
            ##### need to define rmin_input
            dlnpdr = (-p_prime_loc*(8.0π/beta_loc)*
                (inputs.RMIN_LOC/inputs.Q_LOC)*inputs.RMAJ_LOC)
        end

        if dlnpdr > rlnp_cutoff
            dlnpdr = rlnp_cutoff
        end
        if dlnpdr < 4.0
            dlnpdr = 4.0
        end

    else
        dlnpdr = 12.0
    end


    ### coefficients for SAT_RULE = 1
    czf = abs(alpha_zf)

    if (inputs.ETG_FACTOR < 0)
        etg_stiff = abs(inputs.ETG_FACTOR)
      else
        etg_stiff = 1.0
    end



    cnorm = 14.29
    cz1=0.48*czf
    cz2=1.0*czf
    cky=3.0
    sqcky=√(cky)
    cnorm = 14.29
    etg_streamer=1.05
    # if(USE_SUB1)
    #     cnorm=12.12
    #     expsub=1
    # end
    if(alpha_quench != 0.0) etg_streamer=2.1 end
    kyetg = etg_streamer / rho_ion   # fixed to ion gyroradius (Fortran: etg_streamer/rho_ion)

    measure = √(taus_1 * mass_2)
    # if(igeo == 0)then # s-alpha
    #     cnorm=14.63
    #     cz1=0.90*czf
    #     cz2=1.0*czf
    # end

    ### coefficents for SAT_RULE = 2
    if(sat_rule_in == 2 || sat_rule_in == 3)
        b0 = 0.76
        b1 = 1.22
        b2 = 3.74
        if(nmodes > 1) b2 = 3.55 end
        b3 = 1.0

        d1 = (Bt0_out/b_geo0_out)^4
        d1 = d1/grad_r0_out
        Gq = b_geo0_out/grad_r0_out
        d2 = b3/(Gq^2)
        cnorm = b2*(12.0/dlnpdr)
        kyetg = 1000.0
        cky = 3.0
        sqcky = √(cky)
        kycut = b0*kymax_out
        cz1 = 0.0
        cz2 = 1.05*czf
        measure = 1.0/kymax_out
    end
    ### coefficents for SAT_RULE = 3
    if(sat_rule_in==3)
        kmax = kymax_out
        gmax = vzf_out * kymax_out
        kmin = 0.685 * kmax
        aoverb = - 1.0 / (2 * kmin)
        coverb = - 0.751 * kmax
        kT = 1.0/rho_ion
        k0 = 0.6 * kmin
        kP = 2.0 * kmin
        c_1 = - 2.42
        x_ITG = 0.8
        x_TEM = 1.0
        Y_ITG = 3.3 * (gmax^2) / (kmax^5)
        Y_TEM = 12.7 * (gmax^2) / (kmax^4)
       
        scal = 0.82 # Q(SAT3 GA D) / (2 * QLA(ITG,Q) * Q(SAT2 GA D))

        Ys = Vector{T}(undef, nmodes)
        xs = Vector{T}(undef, nmodes)

        for k in 1:nmodes


            sum_W_i = zeros(T, length(QL_weights[1, 1, k, :, 2]))

            for is in 2:size(QL_weights)[2] # sum over ion species, requires electrons to be species 1



                sum_W_i .= sum_W_i .+ QL_weights[1, is, k, :, 2]
            # (field, species, mode, ky, type) type = 2 is energy flux
            end
            
            # check for singularities in weight ratio near kmax
            if(kmax<=ky_spect[1])
                x = 0.5
            else
            ### why is sum_W_i an array now???
                i = 1
                while(ky_spect[i] < kmax)
                    i += 1
                end
                if(sum_W_i[i]==0.0 || sum_W_i[i-1]==0.0)
                    x = 0.5
                else
                    abs_W_ratio = abs.(QL_weights[1,1,k,:,2]./sum_W_i)
                    x = linear_interpolation(ky_spect, abs_W_ratio, kmax)
                end
            end
            Y = mode_transition_function(x, Y_ITG, Y_TEM, x_ITG, x_TEM)
            
            xs[k] = x
            Ys[k] = Y
        end
    end


    # coefficients for spectral shift model for ExB shear
    ax=0.0
    ay=0.0
    exp_ax = 1
    if(alpha_quench == 0.0)
        #spectral shift model parameters
        ax = 1.15
        ay = 0.56
        exp_ax = 4

        if(sat_rule_in==2 || sat_rule_in==3)
            ax = 1.21
            ay = 1.0
            exp_ax = 2
        end
    end

    # if(!first_pass)  # second pass for spectral shift model
    gamma_net = Vector{T}(undef,nky)
    for i in 1:nky
        kx = kx0_e[i]
        if(sat_rule_in==2 || sat_rule_in==3)
            ky0 = ky_spect[i]
            if(ky0<kycut)
                kx_width = kycut/grad_r0_out
            else
                kx_width = kycut/grad_r0_out + b1*(ky0 - kycut)*Gq
            end
            kx = kx*ky0/kx_width
        end
        gamma_net[i] = most_unstable_gamma[i]/(1.0 + abs(ax*kx)^exp_ax)
    end

    if(sat_rule_in==1)
        vzf_out, kymax_out, jmax_out = get_zonal_mixing(inputs, satParams, gamma_net)
    else
        vzf_out_fp = vzf_out # used for recreation of first pass for SAT3
        vzf_out = vzf_out*gamma_net[jmax_out]/max(most_unstable_gamma[jmax_out], small)
    end
    # end   # second pass complete

    gammamax1= vzf_out*kymax_out
    kymax1 = kymax_out
    jmax1 = jmax_out
    vzf1 = vzf_out

    gamma_mix1 = Vector{T}(undef,nky)
    gamma = Vector{T}(undef,nky)

    if(sat_rule_in==1)
        gamma .= ifelse.(ky_spect.<kymax1,
                    max.(gamma_net .- (cz1.*(kymax1 .- ky_spect).*vzf1),0.0),
                    (cz2*gammamax1) .+ etg_stiff*max.(gamma_net .- (cz2.*vzf1.*ky_spect),0.0))
    elseif(sat_rule_in==2 || sat_rule_in==3)
        gamma .= ifelse.(ky_spect.<kymax1,
                    gamma_net,  # Use spectral-shifted values for low-k  
                    gammamax1 .+ etg_stiff*max.(gamma_net .- cz2*vzf1.*ky_spect, 0.0))  # Use gamma_net (not most_unstable_gamma) for high-k
    end
    gamma_mix1 .= gamma


    # if(USE_MIX)
        # mix over ky > kymax with integration weight = sqcky*ky0**2/(ky0**2 + cky*(ky-ky0)**2)
        #### unsure of this index
        mixnorm1 = 0.0
        for j in jmax1+2:nky
            gamma_ave = 0.0
            ky0 = ky_spect[j]
            mixnorm1 = (ky0*
                        (atan(sqcky*(ky_spect[nky]/ky0 - 1.0))
                        - atan(sqcky*(ky_spect[jmax1+1]/ky0 - 1.0)))
            )
            for i in jmax1+1:nky-1
                ky1 = ky_spect[i]
                ky2 = ky_spect[i+1]
                mix1 = (ky0*
                        (atan(sqcky*(ky2/ky0 - 1.0))
                        - atan(sqcky*(ky1/ky0-1.0)))
                )
                delta = (gamma[i+1]-gamma[i])/(ky2-ky1)
                mix2 = (ky0*mix1
                    + (ky0*ky0/(2.0*sqcky))
                    * (log(cky*(ky2-ky0)^2 + ky0^2)
                    -  log(cky*(ky1-ky0)^2 + ky0^2))
                )
                gamma_ave = gamma_ave + (gamma[i] - ky1*delta)*mix1 + delta*mix2
            end
            gamma_mix1[j] = gamma_ave/mixnorm1
        end
    # end

    # recreate first pass for ExB rule in SAT3
    if(sat_rule_in==3)
        # if(!first_pass)
            gamma_fp = similar(ky_spect)
            gamma = similar(ky_spect)

            gamma .= ifelse.(ky_spect.<kymax1,
                             most_unstable_gamma,  # Use original first-pass values for SAT3 recreation  
                            ((gammamax1 * ifelse(vzf_out != 0, vzf_out_fp/vzf_out, 1.0)) .+ max.(most_unstable_gamma .- (cz2.*vzf_out_fp.*ky_spect),0.0)))

            gamma_fp .= gamma
        # end


        # if(USE_MIX)
            for j in jmax1+2:nky
                gamma_ave = 0.0
                ky0 = ky_spect[j]
                mixnorm = (ky0*
                            (atan(sqcky*(ky_spect[nky]/ky0 - 1.0))
                            - atan(sqcky*(ky_spect[jmax1+1]/ky0 - 1.0)))
                )

                for i in jmax1+1:nky-1
                    ky1 = ky_spect[i]
                    ky2 = ky_spect[i+1]
                    mix1 = (ky0*
                            (atan(sqcky*(ky2/ky0-1.0))
                            - atan(sqcky*(ky1/ky0-1.0)))
                    )
                    delta = (gamma[i+1]-gamma[i])/(ky2-ky1)
                    mix2 = (ky0*mix1
                            + (ky0*ky0/(2.0*sqcky))*(log(cky*(ky2-ky0)^2+ky0^2)
                            - log(cky*(ky1-ky0)^2+ky0^2))
                    )
                    gamma_ave = gamma_ave + (gamma[i] - ky1*delta)*mix1 + delta*mix2
                end
                gamma_fp[j] = gamma_ave/mixnorm
            end
        # end
    # else
    #     gamma_fp = gamma_mix1
    end



    # generate SAT2 potential for SAT3 to connect to for electron scale

    # zero-initialized like the Fortran: when the connecting-quadratic denominator
    # vanishes YTs stays 0.0 and the 1e-5 floor below takes over
    YTs = zeros(T, nmodes)
    if(sat_rule_in==3)
        if(ky_spect[nky] >= kT)
            dummy_interp = zeros(T, size(ky_spect)...)
            k = 1
            while ky_spect[k] < kT
                k += 1
            end
            for l in k-1:k
                gamma0 = most_unstable_gamma[l]  # Use original eigenvalues for gamma0
	            ky0 = ky_spect[l]
	            kx = kx0_e[l]

                if(ky0 < kycut)
                    kx_width = kycut/grad_r0_out
                    sat_geo_factor = SAT_geo0_out*d1*SAT_geo1_out
                else
                    kx_width = kycut/grad_r0_out + b1*(ky0 - kycut)*Gq
                    sat_geo_factor = (SAT_geo0_out*
                                (d1*SAT_geo1_out*kycut +
                                (ky0 - kycut)*d2*SAT_geo2_out)/ky0
                    )
                end
	            kx = kx*ky0/kx_width
	            gammaeff = 0.0
	            if(gamma0 > small) gammaeff = gamma_fp[l] end
	            # potentials without multimode and ExB effects, added later
	            dummy_interp[l] = scal*measure*cnorm*(gammaeff/(kx_width*ky0))^2
                if(units_in != "GYRO") dummy_interp[l] = sat_geo_factor*dummy_interp[l] end
            end
	        YT = linear_interpolation(ky_spect, dummy_interp, kT)
	        YTs .= YT
	    else
	        if(aoverb*(kP^2)+kP+coverb-((kP-kT)*(2*aoverb*kP+1))!=0)
	            for l in 1:nmodes
	                YTs[l] = (Ys[l]*
                            (
                                ((aoverb*(k0^2) + k0+coverb)
                                /(
                                    aoverb*(kP^2)
                                    + kP
                                    + coverb
                                    -((kP - kT)*(2*aoverb*kP+1))
                                    )
                                )^abs(c_1))
                    )
                end
            end
        end
	end

    # intensity model

    field_spectrum_out = Matrix{T}(undef, nky,nmodes)
    gammaeff_out = Matrix{T}(undef, nky,nmodes)
    kx_width_out = Vector{T}(undef, nky)
    sat_geo_factor_out = Vector{T}(undef, nky)

    for j in 1:nky
        gamma0 = most_unstable_gamma[j]  # Use original eigenvalue for gamma0 reference
        ky0 = ky_spect[j]
        kx = kx0_e[j]
        if(sat_rule_in==1)
            sat_geo_factor = SAT_geo0_out
            kx_width = ky0
        elseif(sat_rule_in==2 || sat_rule_in==3)
            if(ky0 < kycut)
                kx_width = kycut/ grad_r0_out
                sat_geo_factor = SAT_geo0_out * d1 * SAT_geo1_out
            else
                kx_width = kycut/grad_r0_out + b1*(ky0 - kycut)*Gq
                sat_geo_factor = (SAT_geo0_out*
                                    (d1*SAT_geo1_out*kycut +
                                    (ky0 - kycut)*d2*SAT_geo2_out)/ky0)
            end
            kx = kx*ky0/kx_width
        end
        kx_width_out[j] = kx_width
        sat_geo_factor_out[j] = sat_geo_factor

        if(sat_rule_in==1 || sat_rule_in==2)
            for i in 1:nmodes
                gammaeff = 0.0
                if(gamma0>small)
                    gammaeff = gamma_mix1[j]*(gamma_matrix[i,j]/gamma0)^expsub
                end
                if(ky0>kyetg) gammaeff = gammaeff * √(ky0/kyetg) end

                field_spectrum_out[j,i] = measure*cnorm*((gammaeff/(kx_width*ky0))/(1.0+ay*kx^2))^2
                                
                if(units_in != "GYRO")
                    field_spectrum_out[j,i] = sat_geo_factor*field_spectrum_out[j,i]
                end
                gammaeff_out[j, i] = gammaeff
            end

        elseif(sat_rule_in==3) # SAT3
            if(gamma_fp[j]==0)
                Fky=0.0
            else
		        Fky =  ( (gamma_mix1[j] / gamma_fp[j])^2 ) / ((1.0 + ay*(kx^2))^2)
            end
            for i in 1:nmodes
                field_spectrum_out[j,i] = 0.0
                gammaeff = 0.0
                if(gamma0 > small) # B: Based on the reduction of the matrices field_spectrum_out and gamma_matrix, I'm going to presume that they are in Fortran's [2, :, :] and [1, :, :] space respectively.
                    if (ky0 <= kT)
		                if(kP>=kT)
                           
                            field_spectrum_out[j,i] = 0.0
                        elseif(ky0 <= kP) # initial quadratic
                            sig_ratio = (aoverb * (ky0^2) + ky0 + coverb) / (aoverb * (k0^2) + k0 + coverb)
                            field_spectrum_out[j,i] = Ys[i] * (sig_ratio^c_1) * Fky * (gamma_matrix[i,j]/gamma0)^(2 * expsub)
                        else # connecting quadratic
                            if YTs[i] == 0.0
                               
                                YTs[i] = 1e-5
                            end
                            doversig0 = ( (Ys[i] / YTs[i]) ^ (1.0/abs(c_1))
                                        - (aoverb*kP^2 + kP + coverb - (kP-kT)*(2*aoverb*kP + 1))
                                        /(aoverb*k0^2 + k0 + coverb))
                            doversig0 = doversig0 * (1.0/((kP-kT)^2))
                            eoversig0 = -2 * doversig0 * kP + ((2 * aoverb * kP + 1)/(aoverb * (k0^2) + k0 + coverb))
                            foversig0 = ((Ys[i] / YTs[i])^(1.0/abs(c_1))) - eoversig0 * kT - doversig0 * (kT^2)
                            sig_ratio = doversig0*(ky0^2) + eoversig0*ky0 + foversig0
                            field_spectrum_out[j,i] = Ys[i] * (sig_ratio ^ c_1) * Fky * (gamma_matrix[i,j]/gamma0)^(2*expsub)
                        end

                    else # SAT2 for electron scale

                        gammaeff = gamma_mix1[j]*(gamma_matrix[i,j]/gamma0)^expsub
                        if(ky0 > kyetg) gammaeff = gammaeff*√(ky0/kyetg) end

                        field_spectrum_out[j,i] = (scal*measure*cnorm*
                                        ((gammaeff/(kx_width*ky0))/(1.0+ay*kx^2))^2)
                        if(units_in != "GYRO")
                            field_spectrum_out[j,i] = sat_geo_factor*field_spectrum_out[j,i]
                        end

                    end
                end
                gammaeff_out[j, i] = gammaeff
            end
        end
    end

	#SAT3 QLA here
    if(sat_rule_in==3)
        QLA_P = zeros(T, nmodes)
	    QLA_E = zeros(T, nmodes)
	    for k in 1:nmodes
	        # factor of 2 included for real symmetry
            QLA_P[k] = 2 * mode_transition_function(xs[k], 1.1, 0.6, x_ITG, x_TEM)
            QLA_E[k] = 2 * mode_transition_function(xs[k], .75, 0.6, x_ITG, x_TEM)
        end
	    QLA_O = fill(2.0 * 0.8, nmodes)
    else
	   QLA_P = fill(1.0,nmodes)
	   QLA_E = fill(1.0,nmodes)
	   QLA_O = fill(1.0,nmodes)
    end

    phinorm = field_spectrum_out
    # so the normal behavior doesn't change,
    if return_phi_params
        out = (;
            phinorm      = phinorm,        # [nky, nmodes] SAT intensity at kx peak (= kx0_e)
            kx_width     = kx_width_out,   # [nky] kx_model (kx rms width)
            gammaeff     = gammaeff_out,   # [nky, nmodes] effective growthrate
            kx0_e        = kx0_e,          # [nky] spectral shift in kx
            ax           = ax,             # SAT1 (cx), SAT2 (alpha_x)
            ay           = ay,             # SAT1 (cy)
            exp_ax       = exp_ax,         # SAT2 (sigma_x)
        )
        if sat_rule_in == 2
            # add bonus geometry params,
            out = (; out..., d1=d1, d2=d2, kycut=kycut, b3=b3)
        end
    else
        out = phinorm, QLA_P, QLA_E, QLA_O  # SAT123 intensity and QLA params
    end

    return out

end

"""
    flux_integrals(inputs::InputTJLF, QL::Array{T,5}, QL_flux_out::Array{T,3},q_low_out::Matrix{T},i::Int,ky::T,dky0::T,dky1::T) where T <: Real

description:
    helper function that calculates the flux integrals given the index,
    and ky/QL values among other things
"""
######### i think dky0 is sometimes a different type for some reason #########
function flux_integrals(inputs::InputTJLF, QL::Array{T,5}, QL_flux_out::Array{T,3},q_low_out::Matrix{T},
    i::Int,ky::T,dky0::T,dky1::T) where T <: Real

    # Compute the flux integrals
    # particle
    QL_flux_out[:,:,1] .+= (dky0 .* (i==1 ? 0 : sum(QL[:,:,:,i-1,1],dims=3))
                             .+ dky1 .* sum(QL[:,:,:,i,1],dims=3))
    # energy
    QL_flux_out[:,:,2] .+= (dky0 .* (i==1 ? 0 : sum(QL[:,:,:,i-1,2],dims=3))
                             .+ dky1 .* sum(QL[:,:,:,i,2],dims=3))
    # torodial stress
    QL_flux_out[:,:,3] .+= (dky0 .* (i==1 ? 0 : sum(QL[:,:,:,i-1,3],dims=3))
                             .+ dky1 .* sum(QL[:,:,:,i,3],dims=3))
    # torodial stress
    QL_flux_out[:,:,4] .+= (dky0 .* (i==1 ? 0 : sum(QL[:,:,:,i-1,4],dims=3))
                             .+ dky1 .* sum(QL[:,:,:,i,4],dims=3))
    # exchange
    QL_flux_out[:,:,5] .+= (dky0 .* (i==1 ? 0 : sum(QL[:,:,:,i-1,5],dims=3))
                             .+ dky1 .* sum(QL[:,:,:,i,5],dims=3))

    if ky * inputs.TAUS[1] * inputs.MASS[2] <= 1
        # sum of first and second field (phi, apar) of energy flux; an
        # electrostatic-only QL tensor (nf = 1) has no apar contribution
        if size(QL_flux_out, 1) >= 2
            q_low_out .= QL_flux_out[1,:,2] .+ QL_flux_out[2,:,2]
        else
            q_low_out .= QL_flux_out[1,:,2]
        end
    end

    return QL_flux_out, q_low_out
end

"""
    sum_ky_spectrum(inputs::InputTJLF{T},satParams::SaturationParameters{T},gamma_matrix::Matrix{T},QL::Array{T,5})where T <: Real

parameters:
    inputs              - InputTJLF struct constructed using the input.TGLF file
    gamma_matrix        - matrix of gamma (mode, ky)
    QL_weights          - split into separate types of QL weights (field, species, mode, ky, type)
                          type: (particle, energy, torodial stress, parallel stress, exchange)
    (optional)          - a lot of optional parameters that I don't use -DSUN

outputs:
    QL_flux_out         - flux integral of the QL weights (field, species, type)

    takes in the input.tglf file, ky values, gammas, quasilinear weights, calls
    intensity_sat() which returns the QL intensity values, then calls
    flux_integrals() which numerically integrates the fluxes
"""
function sum_ky_spectrum(
    inputs::InputTJLF{T},
    satParams::SaturationParameters{T},
    gamma_matrix::Matrix{T},
    QL_weights::Array{T,5};
    vzf_out_param=T(NaN),
    kymax_out_param=T(NaN),
    jmax_out_param::Int=0,
    phi_bar_matrix::Matrix{T}=Matrix{T}(undef, 0, 0)
)where T <: Real

    sat_rule_in = inputs.SAT_RULE
    # Size the output by the QL tensor itself, not by fixed constants:
    # QLNN bundles trained with drop_bpar hand in nf < 3, and QLNN packs the
    # canonical (e, D, C) species set which may differ from inputs.NS.
    # Native TJLF callers always build (3, NS, ...) so this is a no-op for them.
    nf = size(QL_weights, 1) # get the number of fields
    ns = size(QL_weights, 2) # get the number of species
    nm = inputs.NMODES # get the number of modes
    ky_spect = inputs.KY_SPECTRUM
    nky = length(ky_spect)

    flux_spectrum = similar(QL_weights)

    # Multiply QL weights with desired intensity
    if sat_rule_in >= 1 && sat_rule_in <= 3
        intensity_factor, QLA_P, QLA_E, QLA_O = intensity_sat(inputs, satParams, gamma_matrix, QL_weights; vzf_out_param=vzf_out_param, kymax_out_param=kymax_out_param, jmax_out_param=jmax_out_param)
        # Ql size (nf,ns,nm,nky,ntype)
        # QLA_P and QLA_E are vectors of size (nm)
        # intensity factor size (nky, nm)
        flux_spectrum[:,:,:,:,1] = QL_weights[:,:,:,:,1] .* reshape((QLA_P .* intensity_factor'),(1,1,nm,nky)) # particle
        flux_spectrum[:,:,:,:,2] = QL_weights[:,:,:,:,2] .* reshape((QLA_E .* intensity_factor'),(1,1,nm,nky)) # energy
        flux_spectrum[:,:,:,:,3] = QL_weights[:,:,:,:,3] .* reshape((QLA_O .* intensity_factor'),(1,1,nm,nky)) # toroidal stress
        flux_spectrum[:,:,:,:,4] = QL_weights[:,:,:,:,4] .* reshape((QLA_O .* intensity_factor'),(1,1,nm,nky)) # parallel stress
        flux_spectrum[:,:,:,:,5] = QL_weights[:,:,:,:,5] .* reshape((QLA_O .* intensity_factor'),(1,1,nm,nky)) # exchange
    elseif sat_rule_in == 0
        @assert !isempty(phi_bar_matrix) "phi_bar_matrix must be provided for SAT_RULE=0"
        phi = reshape(phi_bar_matrix, 1, 1, nm, nky, 1)
        flux_spectrum .= QL_weights .* phi
    else
        throw(error("sat_rule_in must be 0,1,2,or 3, not $sat_rule_in"))
    end
    

    # outputs
    q_low_out = zeros(T, (ns, nm))
    QL_flux_out = zeros(T, nf, ns, 5)
    dky0 = zero(T)
    dky1 = zero(T)
    ky0 = zero(T)
    for i in eachindex(ky_spect)
        ky = ky_spect[i]
        ky1 = ky
        if i == 1
            dky1 = ky1
        elseif inputs.KYGRID_MODEL != 0
            dky = log(ky1 / ky0) / (ky1 - ky0)
            dky1 = ky1 * (1.0 - ky0 * dky)
            dky0 = ky0 * (ky1 * dky - 1.0)
        end

        QL_flux_out, q_low_out = flux_integrals(inputs, flux_spectrum, QL_flux_out, q_low_out,
                                                i,ky,dky0,dky1)

        ky0 = ky1
    end

    return QL_flux_out, flux_spectrum

end
