"""
    tjlf_eigensolver(inputs::InputTJLF{T},outputGeo::OutputGeometry{T},satParams::SaturationParameters{T},ave::Ave{T},aveH::AveH{T},aveWH::AveWH{T},aveKH::AveKH,aveG::AveG{T},aveWG::AveWG{T},aveKG::AveKG,aveGrad::AveGrad{T},aveGradB::AveGradB{T},nbasis::Int, ky::T,amat::Matrix{K},bmat::Matrix{K},ky_index::Int) where T<:Real where K<:Complex

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    ave::Ave{T}                         - structs created in tjlf_matrix.jl
    nbasis::Int                         - used to determine size of the matrix
    ky::T                               - current ky value
    amat/bmat::Matrix{K}                - matrix to be populated
    ky_index::Int                       - index used for multithreading

outputs:
    eigenvalues::Complex                - eigenvalues of the matrices

description:
    uses the structs calculated in tjlf_matrix.jl to populate matrix amat and bmat, solves the generalized eigenvalue problem for only the eigenvalues, returns those eigenvalues
"""
function tjlf_eigensolver(inputs::InputTJLF{T},outputGeo::OutputGeometry{T},satParams::SaturationParameters{T},
                        ave::Ave{T},aveH::AveH{T},aveWH::AveWH{T},aveKH::AveKH,
                        aveG::AveG{T},aveWG::AveWG{T},aveKG::AveKG,
                        aveGrad::AveGrad{T},aveGradB::AveGradB{T},
                        nbasis::Int, ky::T,
                        amat::Matrix{K},bmat::Matrix{K},
                        ky_index::Int,
                        find_eigenvector::Bool; use_gpu::Bool=false,
                        workspace::Union{Nothing,LSWorkspace{T}}=nothing) where T<:Real where K<:Complex

    ft = outputGeo.fts[1]  # electrons
    ft2 = ft^2
    ft3 = ft^3
    ft4 = ft^4
    ft5 = ft^5

    betae_s = inputs.BETAE ##### not true for 'GENE' units
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    nroot = 15 #### hardcoded

    vpar = inputs.VPAR
    mass = inputs.MASS
    rlns = inputs.RLNS
    taus = inputs.TAUS
    zs = inputs.ZS
    as = inputs.AS


    vs2 = √(taus[2] / mass[2])
    vs1 = √(taus[1] / mass[1])

    width_in = inputs.WIDTH_SPECTRUM[ky_index]
    use_bpar_in = inputs.USE_BPAR
    use_bper_in = inputs.USE_BPER
    vpar_model_in = inputs.VPAR_MODEL
    damp_psi_in = inputs.DAMP_PSI
    damp_sig_in = inputs.DAMP_SIG
    linsker_factor_in = inputs.LINSKER_FACTOR
    gchat_in = inputs.GCHAT
    park_in = inputs.PARK
    rmaj_input = inputs.RMAJ_LOC #### might be different for different geometries
    alpha_p_in = inputs.ALPHA_P
    alpha_mach_in = inputs.ALPHA_MACH
    sign_it_in = inputs.SIGN_IT
    xnu_factor_in = inputs.XNU_FACTOR
    wdia_trapped_in = inputs.WDIA_TRAPPED
    gradB_factor_in = inputs.GRADB_FACTOR
    xnue_s = inputs.XNUE
    xnu_model_in = inputs.XNU_MODEL
    ghat_in = inputs.GHAT

    R_unit = satParams.R_unit
    B_unit = satParams.B_unit
    q_unit = satParams.q_unit

    U0 = sum((alpha_mach_in*sign_it_in).*vpar.*zs.^2 .* as./taus) ### defined in startup.f90

    #*************************************************************
    # START
    #*************************************************************

    k_par0 = park_in/(R_unit*q_unit*width_in)
    w_d0 = ky/R_unit
    w_cd = -gchat_in*w_d0
    w_s = -ky/B_unit

    if(!use_bper_in)
        hnb0 = 0.0
        hp1b0 = 0.0
        hp3b0 = 0.0
        hr11b0 = 0.0
        hr13b0 = 0.0
        hr33b0 = 0.0
        hw113b0 = 0.0
        hw133b0 = 0.0
        hw333b0 = 0.0
        kpar_hp1b0 = zero(K)
        kpar_hr11b0 = zero(K)
        kpar_hr13b0 = zero(K)
        wdhp1b0 = 0.0
        wdhr11b0 = 0.0
        wdhr13b0 = 0.0
        gnb0 = 0.0
        gp1b0 = 0.0
        gp3b0 = 0.0
        gr11b0 = 0.0
        gr13b0 = 0.0
        gr33b0 = 0.0
        gw113b0 = 0.0
        gw133b0 = 0.0
        gw333b0 = 0.0
        kpar_gp1b0 = zero(K)
        kpar_gr11b0 = zero(K)
        kpar_gr13b0 = zero(K)
        wdgp1b0 = 0.0
        wdgr11b0 = 0.0
        wdgr13b0 = 0.0
    end
    if(!use_bpar_in)
        h10n = 0.0
        h10p1 = 0.0
        h10p3 = 0.0
        h10r13 = 0.0
        h10r33 = 0.0
        g10n = 0.0
        g10p1 = 0.0
        g10p3 = 0.0
        g10r13 = 0.0
        g10r33 = 0.0
    end
    if(vpar_model_in!=0)
        hnbp = 0.0
        hp1bp = 0.0
        hp3bp = 0.0
        hr11bp = 0.0
        hr13bp = 0.0
        hr33bp = 0.0
        hw113bp = 0.0
        hw133bp = 0.0
        hw333bp = 0.0
        wdhp1bp = 0.0
        wdhr11bp = 0.0
        wdhr13bp = 0.0
        kpar_hnbp = zero(K)
        kpar_hp1bp = zero(K)
        kpar_hp3bp = zero(K)
        kpar_hr11bp = zero(K)
        kpar_hr13bp = zero(K)
        gnbp = 0.0
        gp1bp = 0.0
        gp3bp = 0.0
        gr11bp = 0.0
        gr13bp = 0.0
        gr33bp = 0.0
        gw113bp = 0.0
        gw133bp = 0.0
        gw333bp = 0.0
        wdgp1bp = 0.0
        wdgr11bp = 0.0
        wdgr13bp = 0.0
        kpar_gnbp = zero(K)
        kpar_gp1bp = zero(K)
        kpar_gp3bp = zero(K)
        kpar_gr11bp = zero(K)
        kpar_gr13bp = zero(K)
    end

    betae_psi = 0.0
    damp_psi = 0.0 ##### this is just never assigned a different value....
    if(use_bper_in)
        if(nbasis==2)
            betae_psi = 0.5*betae_s/(ky^2+(damp_psi_in*vs2/(q_unit*width_in))^2)
        else
            betae_psi = 0.5*betae_s/ky^2
        end
    end
    betae_sig = zero(T)
    damp_sig = 0.0 ##### this is just never assigned a different value....
    if(use_bpar_in)
        if(nbasis==2)
            betae_sig = 0.5*betae_s/(ky^2 + (damp_sig_in*vs2/(q_unit*width_in))^2)
        else
            betae_sig = 0.5*betae_s/(ky^2)
        end
    end

    linsker = 0.5*linsker_factor_in
    if(nbasis==1) linsker=0.0 end

    am = 1.0
    bm = 0.0
    if(linsker!=0.0)
        am = am/2
        bm = bm/2
    end


    #  GLF toroidal closure coefficients
    v1_r = uv_constants.v[1]
    v1_i = uv_constants.v[2]
    v2_r = uv_constants.v[3]
    v2_i = uv_constants.v[4]
    v3_r = uv_constants.v[5]
    v3_i = uv_constants.v[6]
    v4_r = uv_constants.v[7]
    v4_i = uv_constants.v[8]
    v5_r = uv_constants.v[9]
    v5_i = uv_constants.v[10]
    v6_r = uv_constants.v[11]
    v6_i = uv_constants.v[12]
    v7_r = uv_constants.v[13]
    v7_i = uv_constants.v[14]
    v8_r = uv_constants.v[15]
    v8_i = uv_constants.v[16]
    v9_r = uv_constants.v[17]
    v9_i = uv_constants.v[18]
    v10_r = uv_constants.v[19]
    v10_i = uv_constants.v[20]

    vb1_r = uv_constants.vb[1]
    vb1_i = uv_constants.vb[2]
    vb2_r = uv_constants.vb[3]
    vb2_i = uv_constants.vb[4]
    vb3_r = uv_constants.vb[5]
    vb3_i = uv_constants.vb[6]
    vb4_r = uv_constants.vb[7]
    vb4_i = uv_constants.vb[8]
    vb5_r = uv_constants.vb[9]
    vb5_i = uv_constants.vb[10]
    vb6_r = uv_constants.vb[11]
    vb6_i = uv_constants.vb[12]
    vb7_r = uv_constants.vb[13]
    vb7_i = uv_constants.vb[14]
    vb8_r = uv_constants.vb[15]
    vb8_i = uv_constants.vb[16]
    vb9_r = uv_constants.vb[17]
    vb9_i = uv_constants.vb[18]
    vb10_r = uv_constants.vb[19]
    vb10_i = uv_constants.vb[20]

    # GLF parallel closure coefficients
    bpar_HP = 3.0 + (32.0 - 9.0π)/(3.0π - 8.0)
    bper_DH = 1.0
    dper_DH = √(π/2.0)
    dpar_HP = 2.0*√(2.0π)/(3.0π - 8.0)

    b1 = bpar_HP
    d1 = dpar_HP
    b3 = bper_DH
    d3 = dper_DH
    b33 = (b1 - b3)/3.0
    d33 = (d1 - d3)/3.0

    # include R(theta)/R0 factor like gyro convetions. Note that sign_Bt_in is in ave_c_tor_par
    vpar_shear = Vector{T}(undef, ns)
    vpar = Vector{T}(undef, ns)
    for is = 1:ns
        vpar_shear_in = inputs.VPAR_SHEAR[is]
        vpar_in = inputs.VPAR[is]

        vpar_shear[is] = vpar_shear_in          *(alpha_p_in   *sign_it_in *ave.c_tor_par[1,1]/rmaj_input)
        if(vpar_model_in==0) vpar[is] = (vpar_in*alpha_mach_in*sign_it_in) *ave.c_tor_par[1,1]/rmaj_input end
    end

    #*************************************************************
    #  compute electron-ion collsion model coefficients
    #*************************************************************
    xnu_model = inputs.XNU_MODEL
    zeff_in = inputs.ZEFF
    k1=0.0
    k2=0.0
    k3=0.0
    k4=0.0
    k5=0.0
    if(xnu_model<=1)
        k1 = 0.601248 + 1.12838*zeff_in
        k2 = 0.797885 + 1.12838*zeff_in
        k3 = 1.795240 + 2.25676*zeff_in
    end
    xnu_p1_1    = -(4/5)*k2
    xnu_u_u_1   = -(2/3)*k1
    xnu_u_q3_1  = -(2/5)*k2 + k1
    xnu_q1_u_1  = -(4/5)*k2
    xnu_q1_q1_1 = -(16/35)*k3
    xnu_q1_q3_1 =  (6/5)*k2 + (12/35)*k3
    xnu_q3_u_1  = -(4/9)*k2
    xnu_q3_q3_1 =  (2/3)*k2 - (4/15)*k3

    if(park_in==0.0)
        xnu_u_u_1   = 0.0
        xnu_u_q3_1  = 0.0
        xnu_q1_u_1  = 0.0
        xnu_q1_q1_1 = 0.0
        xnu_q1_q3_1 = 0.0
        xnu_q3_u_1  = 0.0
        xnu_q3_q3_1 = 0.0
    end
    if(nroot>6)

        xnu_hat = xnue_s/(ky*taus[1]/R_unit)

        # model for effective ion wavenumber averaged with ion charge densities zs*as
        # ki = k_theta*sqrt(sum_s(rho_s**2*as*zs))

        ki = 0.0
        charge_tot = 0.0
        for is = 2:ns
            ki = ki + taus[is]*mass[is]*as[is]*zs[is]
            charge_tot = charge_tot + as[is]*zs[is]
        end
        ki = √(ki/charge_tot)*ky
        ks0 = ky*√(taus[1]*mass[2])


        xnu_phi_b = 0.0
        if(xnu_model==1) xnu_phi_b = 1.0 end
        if(xnu_phi_b==0.0)
            # boundary collision model without phi terms
            gradne = rlns[1]*R_unit
            gradne_s = max(gradne+10.8, 1.8)
            xnu_c = gradne_s*(1.5*(1.0-tanh((gradne_s/12.6)^2))+0.13)

            xnu_a = ki*(max(0.36+0.10*gradne,0.0) + xnu_c*(ki/ks0)*(1.0-tanh(ki/0.55)))
            xnu_b = 3.1/(1.0+(2.1*ki + 8.0*ki^2)*xnu_hat)
        else
            # boundary collsion model with phi terms
            xnu_a = 0.41 + ((ki/ks0)^1.7)*0.70*(1.0 + 1.4*ki/0.38)/(1.0 + (ki/0.38)^4)
            xnu_b = 3.79/(1.0 + 4.63*ki*xnu_hat)
        end
        xnu_bndry = (1.0 - ft2)*xnu_a*xnu_b

        xnu_n_b     = xnu_factor_in*xnu_bndry*xnu_p1_1
        xnu_p3_b    = xnu_n_b
        xnu_p1_b    = xnu_n_b
        xnu_u_b     = xnu_factor_in*xnu_bndry*xnu_q1_q1_1
        xnu_q3_b    = xnu_u_b
        xnu_q1_b    = xnu_u_b

        if(park_in==0.0)
            xnu_u_b = 0.0
            xnu_q1_b=0.0
            xnu_q3_b =0.0
        end
    end

    cnuei = 0.0
    if(xnu_model>=2) cnuei = xnue_s end
    kparvthe = abs(k_par0)*vs1/√(2.0)
    kparvthe = max(kparvthe,1.0E-10)

    #*************************************************************
    # full velocity space terms
    #*************************************************************
    #### from tglf_modules.f90
    nuei_c1 = [0.4919, 0.7535, 0.6727, 0.8055, 1.627,
                2.013, 0.4972, 0.7805, 1.694, 3.103,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    k1 = cnuei*(nuei_c1[1] +zeff_in*nuei_c1[2])
    k2 = cnuei*(nuei_c1[3] +zeff_in*nuei_c1[4])
    k3 = cnuei*(nuei_c1[5] +zeff_in*nuei_c1[6])
    k4 = cnuei*(nuei_c1[7] +zeff_in*nuei_c1[8])
    k5 = cnuei*(nuei_c1[9] +zeff_in*nuei_c1[10])
    c01 = (4/5)*k4
    c02 = (2/5)*k2 - k1
    c03 = (2/3)*k1
    c04 = (4/15)*k3 - (4/3)*k2+(5/3)*k1
    c05 = (16/35)*k5
    c06 = 0.0
    c07 = 0.0
    c08 = 0.0
    c09 = 0.0
    c010 = 0.0

    nuei_p1_p1_1 = c01
    nuei_p1_p3_1 = -nuei_p1_p1_1
    nuei_u_q3_1  = c02
    nuei_u_u_1   = c03 - (5/3)*nuei_u_q3_1
    nuei_q3_q3_1 = c04
    nuei_q3_u_1  = (10/9)*c02 - (5/3)*nuei_q3_q3_1
    nuei_q3_u_1  = nuei_q3_u_1 + (5/3)*nuei_u_u_1
    nuei_q3_q3_1 = nuei_q3_q3_1 + (5/3)*nuei_u_q3_1
    nuei_q1_q1_1 = c05
    nuei_q1_q3_1 = -(9/5)*nuei_q1_q1_1
    nuei_q1_u_1  = (9/5)*nuei_q3_u_1
    nuei_q1_q3_1 = nuei_q1_q3_1 + (9/5)*nuei_q3_q3_1

    nuei_p1_p1_t = c01
    nuei_p1_p3_t = -ft2*nuei_p1_p1_t
    nuei_u_q3_t  = c02
    nuei_u_u_t   = c03 -(5.0/3.0)*nuei_u_q3_t
    nuei_q3_q3_t = c04
    nuei_q3_u_t  = (10.0/9.0)*c02 -(5.0/3.0)*nuei_q3_q3_t
    nuei_q3_u_t  = nuei_q3_u_t + (5.0/3.0)*nuei_u_u_t
    nuei_q3_q3_t = nuei_q3_q3_t + (5.0/3.0)*nuei_u_q3_t
    nuei_q1_q1_t = c05
    nuei_q1_q3_t = -(9.0/5.0)*ft2*nuei_q1_q1_t
    nuei_q1_u_t  = (9.0/5.0)*ft2*nuei_q3_u_t
    nuei_q1_q3_t = nuei_q1_q3_t + (9.0/5.0)*ft2*nuei_q3_q3_t

    an = 0.75
    ap3 = 1.25
    ap1 = 2.25
    bn = ft
    bp3 = ft
    bp1 = ft3
    cb3=0.0
    cb5 =0.0

    _cb1_arg = kparvthe*cnuei*(1.0 + 0.82*zeff_in)
    cb1 = _cb1_arg > 0 ? 0.163*√(_cb1_arg) : zero(_cb1_arg)
    if(xnu_model_in==3)
        _cb1_arg2 = cnuei*(1.0 + 0.82*zeff_in)
        if(wdia_trapped_in==0.0)
            cb1 = _cb1_arg2 > 0 ? 0.50*(kparvthe^0.34)*(_cb1_arg2)^0.66 : zero(_cb1_arg2)
        else
            cb1 = _cb1_arg2 > 0 ? 0.315*(kparvthe^0.34)*(_cb1_arg2)^0.66 : zero(_cb1_arg2)
        end
    end
    cb1 = cb1*xnu_factor_in
    cb2 = cb1
    cb4 = cb1

    # even trapped region terms
    nuei_n_n   = (1.0 - ft2)*cb1
    nuei_n_p3  = (1.0 - ft2)*cnuei*cb3
    nuei_n_p1  = 0.0
    nuei_n_p3  = nuei_n_p3 - ft2*nuei_n_p1
    nuei_n_n   = nuei_n_n - nuei_n_p3 - ft2*nuei_n_p1
    nuei_p3_n  = (2.0/3.0)*(1.0 - ft2)*cnuei*cb3
    nuei_p3_p3 = (1.0 -ft2)*cb2
    nuei_p3_p1 = 0.0
    nuei_p3_p3 = nuei_p3_p3 - ft2*nuei_p3_p1
    nuei_p3_n  = nuei_p3_n - nuei_p3_p3 -ft2*nuei_p3_p1
    nuei_p3_n  = nuei_p3_n + nuei_n_n
    nuei_p3_p3 = nuei_p3_p3 + nuei_n_p3
    nuei_p3_p1 = nuei_p3_p1 + nuei_n_p1
    nuei_p1_n  = 0.0
    nuei_p1_p3 = 0.0
    nuei_p1_p1 = (1.0 - ft2)*cb4
    nuei_p1_p3 = nuei_p1_p3 - ft2*nuei_p1_p1
    nuei_p1_n  = nuei_p1_n - nuei_p1_p3 -ft2*nuei_p1_p1
    nuei_p1_n  = nuei_p1_n + ft2*nuei_p3_n
    nuei_p1_p3 = nuei_p1_p3 + ft2*nuei_p3_p3
    nuei_p1_p1 = nuei_p1_p1 + ft2*nuei_p3_p1

    cb5 = 0.0
    cb6 = 0.0
    cb7 = 0.0
    cb8 = 0.0
    nuei_u_u = (1.0 -ft2)*cnuei*cb5
    nuei_u_q3 = (1.0 -ft2)*cnuei*cb7

    nuei_u_q1 = 0.0
    nuei_u_q3 = nuei_u_q3 - (9.0/5.0)*ft2*nuei_u_q1
    nuei_u_u = nuei_u_u -(5.0/3.0)*nuei_u_q3 -3.0*ft2*nuei_u_q1
    nuei_q3_u = (10.0/9.0)*(1.0 -ft2)*cnuei*cb7
    nuei_q3_q3 = (1.0 - ft2)*cnuei*cb6

    nuei_q3_q1 = 0.0
    nuei_q3_q3 = nuei_q3_q3 -(9.0/5.0)*ft2*nuei_q3_q1
    nuei_q3_u = nuei_q3_u - (5.0/3.0)*nuei_q3_q3 -3.0*ft2*nuei_q3_q1
    nuei_q3_u = nuei_q3_u + (5.0/3.0)*nuei_u_u
    nuei_q3_q3 = nuei_q3_q3 +(5.0/3.0)*nuei_u_q3
    nuei_q3_q1 = nuei_q3_q1 +(5.0/3.0)*nuei_u_q1
    nuei_q1_u = 0.0
    nuei_q1_q3 = 0.0
    nuei_q1_q1 = (1.0 -ft2)*cnuei*cb8
    nuei_q1_q3 = nuei_q1_q3 - (9.0/5.0)*ft2*nuei_q1_q1
    nuei_q1_u = nuei_q1_u - (5.0/3.0)*nuei_q1_q3 - 3.0*ft2*nuei_q1_q1
    nuei_q1_u = nuei_q1_u +(9.0/5.0)*ft2*nuei_q3_u
    nuei_q1_q3 = nuei_q1_q3 +(9.0/5.0)*ft2*nuei_q3_q3
    nuei_q1_q1 = nuei_q1_q1 + (9.0/5.0)*ft2*nuei_q3_q1

    #*************************************************************
    # done with electron-ion collision model
    #*************************************************************



    #*************************************************************
    # start of loop over species is,js for amat
    #*************************************************************
    # amat = Matrix{ComplexF64}(undef, iur, iur)
    # bmat = Matrix{ComplexF64}(undef, iur, iur)
    if(nroot>6)
        v = Vector{T}(undef, 20)
        vb = Vector{T}(undef, 20)
    end
    for is = ns0:ns
        rlnsIS = inputs.RLNS[is]
        rltsIS = inputs.RLTS[is]
        tausIS = inputs.TAUS[is]
        massIS = inputs.MASS[is]
        vsIS = √(tausIS / massIS)
        zsIS = inputs.ZS[is]

        ### i hate it here (this is necessary)
        ft = outputGeo.fts[is]  # electrons
        ft2 = ft^2
        ft3 = ft^3
        ft4 = ft^4
        ft5 = ft^5

        if(nroot>6)
            index = floor(Int, 20*ft) + 1
            if(index==21) index=20 end
            df = (ft-uv_constants.fm[index])/(uv_constants.fm[index+1]-uv_constants.fm[index])

            @views @. v = uv_constants.vm[:,index] + uv_constants.vm_diff[:,index] * df
            @views @. vb = uv_constants.vbm[:,index] + uv_constants.vbm_diff[:,index] * df

            u1_r = v[1]
            u1_i = v[2]
            u2_r = v[3]
            u2_i = v[4]
            u3_r = v[5]
            u3_i = v[6]
            u4_r = v[7]
            u4_i = v[8]
            u5_r = v[9]
            u5_i = v[10]
            u6_r = v[11]
            u6_i = v[12]
            u7_r = v[13]
            u7_i = v[14]
            u8_r = v[15]
            u8_i = v[16]
            u9_r = v[17]
            u9_i = v[18]
            u10_r = v[19]
            u10_i = v[20]

            ub1_r = vb[1]
            ub1_i = vb[2]
            ub2_r = vb[3]
            ub2_i = vb[4]
            ub3_r = vb[5]
            ub3_i = vb[6]
            ub4_r = vb[7]
            ub4_i = vb[8]
            ub5_r = vb[9]
            ub5_i = vb[10]
            ub6_r = vb[11]
            ub6_i = vb[12]
            ub7_r = vb[13]
            ub7_i = vb[14]
            ub8_r = vb[15]
            ub8_i = vb[16]
            ub9_r = vb[17]
            ub9_i = vb[18]
            ub10_r = vb[19]
            ub10_i = vb[20]

            u2_r = u2_r*ft2
            u2_i = u2_i*ft2
            u3_r = u3_r/ft2
            u3_i = u3_i/ft2
            u5_r = u5_r*ft2
            u5_i = u5_i*ft2
            u7_r = u7_r*ft2
            u7_i = u7_i*ft2
            u9_r = u9_r/ft2
            u9_i = u9_i/ft2
            ub2_r = ub2_r*ft2
            ub2_i = ub2_i*ft2
            ub3_r = ub3_r/ft2
            ub3_i = ub3_i/ft2
            ub5_r = ub5_r*ft2
            ub5_i = ub5_i*ft2
            ub7_r = ub7_r*ft2
            ub7_i = ub7_i*ft2
            ub9_r = ub9_r/ft2
            ub9_i = ub9_i/ft2
        end

        #*************************************************************
        # start of loop over basis ib,jb for amat
        #*************************************************************
        for js = ns0:ns
            tausJS = inputs.TAUS[js]
            massJS = inputs.MASS[js]
            vsJS = √(tausJS / massJS)
            zsJS = inputs.ZS[js]
            asJS = inputs.AS[js]

            for ib = 1:nbasis
                for jb = 1:nbasis

                    #*************************************************************
                    # collision  model
                    #*************************************************************
                    xnuei = 0.0
                    if(is==1) xnuei = xnue_s end
                    xnuion = 0.0
                    d_ab = 0.0
                    d_ij = 0.0
                    d_ee = 0.0
                    d_ab_psi = 0.0
                    d_11 = 0.0
                    d_1 = 0.0
                    if(is==1)              d_1=1.0 end
                    if(is==1 && js==1)     d_11=1.0 end
                    if(is==js)             d_ij=1.0 end
                    if(ib==jb && is==js)   d_ab=1.0 end
                    if(is==1 && d_ab==1.0) d_ee=1.0 end
                    if(ib==jb)             d_ab_psi=1.0 end

                    b1 = bpar_HP
                    d1 = dpar_HP
                    b3 = bper_DH
                    d3 = dper_DH
                    bs = 20.0*√(tausIS*massIS)*ky/abs(zsIS)
                    b33 = (b1 - b3)/(3.0)
                    d33 = (d1 - d3)/(3.0)

                    hn = aveH.hnp0[is,ib,jb]
                    hp1 = aveH.hp1p0[is,ib,jb]
                    hp3 = aveH.hp3p0[is,ib,jb]
                    hr11 = aveH.hr11p0[is,ib,jb]
                    hr13 = aveH.hr13p0[is,ib,jb]
                    hr33 = aveH.hr33p0[is,ib,jb]
                    c_tor_par_hp1 = aveH.c_tor_par_hp1p0[is,ib,jb]/rmaj_input
                    c_tor_par_hr11 = aveH.c_tor_par_hr11p0[is,ib,jb]/rmaj_input
                    c_tor_par_hr13 = aveH.c_tor_par_hr13p0[is,ib,jb]/rmaj_input

                    if(use_bper_in)
                        hnb0 = aveH.hnb0[is,ib,jb]
                        hp1b0 = aveH.hp1b0[is,ib,jb]
                        hp3b0 = aveH.hp3b0[is,ib,jb]
                        hr11b0 = aveH.hr11b0[is,ib,jb]
                        hr13b0 = aveH.hr13b0[is,ib,jb]
                        hr33b0 = aveH.hr33b0[is,ib,jb]
                        hw113b0 = aveH.hw113b0[is,ib,jb]
                        hw133b0 = aveH.hw133b0[is,ib,jb]
                        hw333b0 = aveH.hw333b0[is,ib,jb]
                        if(vpar_model_in==0)
                            hnbp = aveH.hnbp[is,ib,jb]
                            hp1bp = aveH.hp1bp[is,ib,jb]
                            hp3bp = aveH.hp3bp[is,ib,jb]
                            hr11bp = aveH.hr11bp[is,ib,jb]
                            hr13bp = aveH.hr13bp[is,ib,jb]
                            hr33bp = aveH.hr33bp[is,ib,jb]
                            hw113bp = aveH.hw113bp[is,ib,jb]
                            hw133bp = aveH.hw133bp[is,ib,jb]
                            hw333bp = aveH.hw333bp[is,ib,jb]
                        end
                    end
                    if(use_bpar_in)
                        h10n = 1.5*(hnb0-hp3b0)
                        h10p1 = 2.5*hp1b0 - 1.5*hr13b0
                        h10p3 = 2.5*hp3b0 - 1.5*hr33b0
                        h10r13 = 3.5*hr13b0 - 1.5*hw133b0
                        h10r33 = 3.5*hr33b0 - 1.5*hw333b0
                    end
                    hu1 = aveH.hu1[is,ib,jb]
                    hu3 = aveH.hu3[is,ib,jb]
                    ht1 = aveH.ht1[is,ib,jb]
                    ht3 = aveH.ht3[is,ib,jb]
                    wdhp1p0 = aveWH.wdhp1p0[is,ib,jb]
                    wdhr11p0 = aveWH.wdhr11p0[is,ib,jb]
                    wdhr13p0 = aveWH.wdhr13p0[is,ib,jb]
                    if(use_bper_in)
                        wdhp1b0 = aveWH.wdhp1b0[is,ib,jb]
                        wdhr11b0 = aveWH.wdhr11b0[is,ib,jb]
                        wdhr13b0 = aveWH.wdhr13b0[is,ib,jb]
                        if(vpar_model_in==0)
                            wdhp1bp = aveWH.wdhp1bp[is,ib,jb]
                            wdhr11bp = aveWH.wdhr11bp[is,ib,jb]
                            wdhr13bp = aveWH.wdhr13bp[is,ib,jb]
                        end
                    end
                    wdhu1 = aveWH.wdhu1[is,ib,jb]
                    wdhu3 = aveWH.wdhu3[is,ib,jb]
                    wdhu3ht1 = aveWH.wdhu3ht1[is,ib,jb]
                    wdhu3ht3 = aveWH.wdhu3ht3[is,ib,jb]
                    wdhu33 = aveWH.wdhu33[is,ib,jb]
                    wdhu33ht1 = aveWH.wdhu33ht1[is,ib,jb]
                    wdhu33ht3 = aveWH.wdhu33ht3[is,ib,jb]
                    modwdhu1 = aveWH.modwdhu1[is,ib,jb]
                    modwdhu3 = aveWH.modwdhu3[is,ib,jb]
                    modwdhu3ht1 = aveWH.modwdhu3ht1[is,ib,jb]
                    modwdhu3ht3 = aveWH.modwdhu3ht3[is,ib,jb]
                    modwdhu33 = aveWH.modwdhu33[is,ib,jb]
                    modwdhu33ht1 = aveWH.modwdhu33ht1[is,ib,jb]
                    modwdhu33ht3 = aveWH.modwdhu33ht3[is,ib,jb]

                    hv1r = (v1_r-vb1_r)*ave.modwdh[ib,jb] + vb1_r*modwdhu3*(3/5)
                    hv2r = (v2_r-vb2_r)*ave.modwdh[ib,jb] + vb2_r*modwdhu3*(3/5)
                    hv3r = (v3_r-vb3_r)*ave.modwdh[ib,jb] + vb3_r*modwdhu33*(3/5)
                    hv4r = (v4_r-vb4_r)*ave.modwdh[ib,jb] + vb4_r*modwdhu33*(3/5)
                    hv5r = (v5_r-vb5_r)*ave.modwdh[ib,jb] + vb5_r*modwdhu3*(3/5)
                    hv6r = (v6_r-vb6_r)*ave.modwdh[ib,jb] + vb6_r*modwdhu3*(3/5)
                    hv7r = (v7_r-vb7_r)*ave.modwdh[ib,jb] + vb7_r*modwdhu3*(3/5)
                    hv8r = (v8_r-vb8_r)*ave.modwdh[ib,jb] + vb8_r*modwdhu33*(3/5)
                    hv9r = (v9_r-vb9_r)*ave.modwdh[ib,jb] + vb9_r*modwdhu33*(3/5)
                    hv10r = (v10_r-vb10_r)*ave.modwdh[ib,jb] + vb10_r*modwdhu33*(3/5)
                    hv1rht1 = (v1_r-vb1_r)*aveWH.modwdht1[is,ib,jb] + vb1_r*modwdhu3ht1*(3/5)
                    hv2rht3 = (v2_r-vb2_r)*aveWH.modwdht3[is,ib,jb] + vb2_r*modwdhu3ht3*(3/5)
                    hv3rht1 = (v3_r-vb3_r)*aveWH.modwdht1[is,ib,jb] + vb3_r*modwdhu33ht1*(3/5)
                    hv4rht3 = (v4_r-vb4_r)*aveWH.modwdht3[is,ib,jb]  + vb4_r*modwdhu33ht3*(3/5)
                    hv6rhu1 = (v6_r*modwdhu1)
                    hv7rhu3 = (v7_r*modwdhu3)
                    hv9rhu1 = (v9_r*modwdhu1)
                    hv10rhu3= (v10_r*modwdhu3)

                    hv1i = (v1_i-vb1_i)*ave.wdh[ib,jb] + vb1_i*wdhu3*(3/5)
                    hv2i = (v2_i-vb2_i)*ave.wdh[ib,jb] + vb2_i*wdhu3*(3/5)
                    hv3i = (v3_i-vb3_i)*ave.wdh[ib,jb] + vb3_i*wdhu33*(3/5)
                    hv4i = (v4_i-vb4_i)*ave.wdh[ib,jb] + vb4_i*wdhu33*(3/5)
                    hv5i = (v5_i-vb5_i)*ave.wdh[ib,jb] + vb5_i*wdhu3*(3/5)
                    hv6i = (v6_i-vb6_i)*ave.wdh[ib,jb] + vb6_i*wdhu3*(3/5)
                    hv7i = (v7_i-vb7_i)*ave.wdh[ib,jb] + vb7_i*wdhu3*(3/5)
                    hv8i = (v8_i-vb8_i)*ave.wdh[ib,jb] + vb8_i*wdhu33*(3/5)
                    hv9i = (v9_i-vb9_i)*ave.wdh[ib,jb] + vb9_i*wdhu33*(3/5)
                    hv10i = (v10_i-vb10_i)*ave.wdh[ib,jb] + vb10_i*wdhu33*(3/5)
                    hv1iht1 = (v1_i-vb1_i)*aveWH.wdht1[is,ib,jb] + vb1_i*wdhu3ht1*(3/5)
                    hv2iht3 = (v2_i-vb2_i)*aveWH.wdht3[is,ib,jb] + vb2_i*wdhu3ht3*(3/5)
                    hv3iht1 = (v3_i-vb3_i)*aveWH.wdht1[is,ib,jb] + vb3_i*wdhu33ht1*(3/5)
                    hv4iht3 = (v4_i-vb4_i)*aveWH.wdht3[is,ib,jb] + vb4_i*wdhu33ht3*(3/5)
                    hv6ihu1 = (v6_i*wdhu1)
                    hv7ihu3 = (v7_i*wdhu3)
                    hv9ihu1 = (v9_i*wdhu1)
                    hv10ihu3= (v10_i*wdhu3)

                    kpar_hnp0 = k_par0*aveKH.kparhnp0[is,ib,jb]
                    kpar_hp1p0 = k_par0*aveKH.kparhp1p0[is,ib,jb]
                    kpar_hp3p0 = k_par0*aveKH.kparhp3p0[is,ib,jb]
                    if(use_bper_in)
                        kpar_hp1b0 = k_par0*aveKH.kparhp1b0[is,ib,jb]
                        kpar_hr11b0 = k_par0*aveKH.kparhr11b0[is,ib,jb]
                        kpar_hr13b0 = k_par0*aveKH.kparhr13b0[is,ib,jb]
                        if(vpar_model_in==0)
                            kpar_hnbp = k_par0*aveKH.kparhnbp[is,ib,jb]
                            kpar_hp1bp = k_par0*aveKH.kparhp1bp[is,ib,jb]
                            kpar_hp3bp = k_par0*aveKH.kparhp3bp[is,ib,jb]
                            kpar_hr11bp = k_par0*aveKH.kparhr11bp[is,ib,jb]
                            kpar_hr13bp = k_par0*aveKH.kparhr13bp[is,ib,jb]
                        end
                    end
                    kpar_hu1 = aveKH.kparhu1[is,ib,jb]
                    kpar_hu3 = aveKH.kparhu3[is,ib,jb]
                    kpar_hb1 = b1*ave.kpar_eff[is,ib,jb]
                    kpar_hb3 = b3*ave.kpar_eff[is,ib,jb]
                    kpar_hb33 = b33*ave.kpar_eff[is,ib,jb]
                    kpar_hb1ht1 = b1*aveKH.kparht1[is,ib,jb]
                    kpar_hb3ht3 = b3*aveKH.kparht3[is,ib,jb]
                    kpar_hb33ht1 = b33*aveKH.kparht1[is,ib,jb]
                    modkpar_hd1 = d1*ave.modkpar_eff[is,ib,jb]
                    modkpar_hd3 = d3*ave.modkpar_eff[is,ib,jb]
                    modkpar_hd33 = d33*ave.modkpar_eff[is,ib,jb]
                    modkpar_hd1hu1 = d1*aveKH.modkparhu1[is,ib,jb]
                    modkpar_hd3hu3 = d3*aveKH.modkparhu3[is,ib,jb]
                    modkpar_hd33hu1 = d33*aveKH.modkparhu1[is,ib,jb]
                    grad_hu1 = 0.0
                    grad_hu3 = 0.0
                    dhr13 = -b3*(aveKH.kparht3[is,ib,jb] - aveKH.kparht1[is,ib,jb]/3 - aveKH.kparhu3[is,ib,jb] + aveKH.kparhu1[is,ib,jb]/3)

                    if(linsker==0.0)
                        gradhp1 = 0.0
                        gradhr11 = 0.0
                        gradhr13 = 0.0
                        gradhp1p1 = 0.0
                        gradhr11p1 = 0.0
                        gradhr13p1 = 0.0
                    else

                        gradhp1 = linsker*aveGrad.gradhp1p0[is,ib,jb]
                        gradhr11 = linsker*aveGrad.gradhr11p0[is,ib,jb]
                        gradhr13 = linsker*aveGrad.gradhr13p0[is,ib,jb]
                        gradhp1p1 = linsker*aveGrad.gradhp1p1[is,ib,jb]
                        gradhr11p1 = linsker*aveGrad.gradhr11p1[is,ib,jb]
                        gradhr13p1 = linsker*aveGrad.gradhr13p1[is,ib,jb]
                    end

                    if(nroot>6)
                        gn = aveG.gnp0[is,ib,jb]
                        gp1 = aveG.gp1p0[is,ib,jb]
                        gp3 = aveG.gp3p0[is,ib,jb]
                        gr11 = aveG.gr11p0[is,ib,jb]
                        gr13 = aveG.gr13p0[is,ib,jb]
                        gr33 = aveG.gr33p0[is,ib,jb]
                        c_tor_par_gp1 = aveG.c_tor_par_gp1p0[is,ib,jb]/rmaj_input
                        c_tor_par_gr11 = aveG.c_tor_par_gr11p0[is,ib,jb]/rmaj_input
                        c_tor_par_gr13 = aveG.c_tor_par_gr13p0[is,ib,jb]/rmaj_input
                        if(use_bper_in)
                            gnb0 = aveG.gnb0[is,ib,jb]
                            gp1b0 = aveG.gp1b0[is,ib,jb]
                            gp3b0 = aveG.gp3b0[is,ib,jb]
                            gr11b0 = aveG.gr11b0[is,ib,jb]
                            gr13b0 = aveG.gr13b0[is,ib,jb]
                            gr33b0 = aveG.gr33b0[is,ib,jb]
                            gw113b0 = aveG.gw113b0[is,ib,jb]
                            gw133b0 = aveG.gw133b0[is,ib,jb]
                            gw333b0 = aveG.gw333b0[is,ib,jb]
                            if(vpar_model_in==0)
                                gnbp = aveG.gnbp[is,ib,jb]
                                gp1bp = aveG.gp1bp[is,ib,jb]
                                gp3bp = aveG.gp3bp[is,ib,jb]
                                gr11bp = aveG.gr11bp[is,ib,jb]
                                gr13bp = aveG.gr13bp[is,ib,jb]
                                gr33bp = aveG.gr33bp[is,ib,jb]
                                gw113bp = aveG.gw113bp[is,ib,jb]
                                gw133bp = aveG.gw133bp[is,ib,jb]
                                gw333bp = aveG.gw333bp[is,ib,jb]
                            end
                        end
                        if(use_bpar_in)
                            g10n = 1.5*(gnb0-gp3b0)
                            g10p1 = 2.5*gp1b0 - 1.5*gr13b0
                            g10p3 = 2.5*gp3b0 - 1.5*gr33b0
                            g10r13 = 3.5*gr13b0 - 1.5*gw133b0
                            g10r33 = 3.5*gr33b0 - 1.5*gw333b0
                        end

                        gu1 = aveG.gu1[is,ib,jb]
                        gu3 = aveG.gu3[is,ib,jb]
                        gt1 = aveG.gt1[is,ib,jb]
                        gt3 = aveG.gt3[is,ib,jb]

                        wdgp1p0 = aveWG.wdgp1p0[is,ib,jb]
                        wdgr11p0 = aveWG.wdgr11p0[is,ib,jb]
                        wdgr13p0 = aveWG.wdgr13p0[is,ib,jb]
                        if(use_bper_in)
                            wdgp1b0 = aveWG.wdgp1b0[is,ib,jb]
                            wdgr11b0 = aveWG.wdgr11b0[is,ib,jb]
                            wdgr13b0 = aveWG.wdgr13b0[is,ib,jb]
                            if(vpar_model_in==0)
                                wdgp1bp = aveWG.wdgp1bp[is,ib,jb]
                                wdgr11bp = aveWG.wdgr11bp[is,ib,jb]
                                wdgr13bp = aveWG.wdgr13bp[is,ib,jb]
                            end
                        end

                        wdgu1 = aveWG.wdgu1[is,ib,jb]
                        wdgu3 = aveWG.wdgu3[is,ib,jb]
                        wdgu33 = aveWG.wdgu33[is,ib,jb]
                        wdgu3gt1 = aveWG.wdgu3gt1[is,ib,jb]
                        wdgu3gt3 = aveWG.wdgu3gt3[is,ib,jb]
                        wdgu33gt1 = aveWG.wdgu33gt1[is,ib,jb]
                        wdgu33gt3 = aveWG.wdgu33gt3[is,ib,jb]
                        modwdgu1 = aveWG.modwdgu1[is,ib,jb]
                        modwdgu3 = aveWG.modwdgu3[is,ib,jb]
                        modwdgu33 = aveWG.modwdgu33[is,ib,jb]
                        modwdgu3gt1 = aveWG.modwdgu3gt1[is,ib,jb]
                        modwdgu3gt3 = aveWG.modwdgu3gt3[is,ib,jb]
                        modwdgu33gt1 = aveWG.modwdgu33gt1[is,ib,jb]
                        modwdgu33gt3 = aveWG.modwdgu33gt3[is,ib,jb]
                        gu1r = (u1_r-ub1_r)*ave.modwdg[ib,jb]+ub1_r*modwdgu3*(3/5)
                        gu2r = (u2_r-ub2_r)*ave.modwdg[ib,jb]+ub2_r*modwdgu3*(3/5)
                        gu3r = (u3_r-ub3_r)*ave.modwdg[ib,jb]+ub3_r*modwdgu33*(3/5)
                        gu4r = (u4_r-ub4_r)*ave.modwdg[ib,jb]+ub4_r*modwdgu33*(3/5)
                        gu5r = (u5_r-ub5_r)*ave.modwdg[ib,jb]+ub5_r*modwdgu3*(3/5)
                        gu6r = (u6_r-ub6_r)*ave.modwdg[ib,jb]+ub6_r*modwdgu3*(3/5)
                        gu7r = (u7_r-ub7_r)*ave.modwdg[ib,jb]+ub7_r*modwdgu3*(3/5)
                        gu8r = (u8_r-ub8_r)*ave.modwdg[ib,jb]+ub8_r*modwdgu33*(3/5)
                        gu9r = (u9_r-ub9_r)*ave.modwdg[ib,jb]+ub9_r*modwdgu33*(3/5)
                        gu10r = (u10_r-ub10_r)*ave.modwdg[ib,jb]+ub10_r*modwdgu33*(3/5)
                        gu1rgt1 = (u1_r-ub1_r)*aveWG.modwdgt1[is,ib,jb]+ub1_r*modwdgu3gt1*(3/5)
                        gu2rgt3 = (u2_r-ub2_r)*aveWG.modwdgt3[is,ib,jb]+ub2_r*modwdgu3gt3*(3/5)
                        gu3rgt1 = (u3_r-ub3_r)*aveWG.modwdgt1[is,ib,jb]+ub3_r*modwdgu33gt1*(3/5)
                        gu4rgt3 = (u4_r-ub4_r)*aveWG.modwdgt3[is,ib,jb]+ub4_r*modwdgu33gt3*(3/5)
                        gu6rgu1 = (u6_r*modwdgu1)
                        gu7rgu3 = (u7_r*modwdgu3)
                        gu9rgu1 = (u9_r*modwdgu1)
                        gu10rgu3= (u10_r*modwdgu3)
                        gu1i = (u1_i-ub1_i)*ave.wdg[ib,jb]+ub1_i*wdgu3*(3/5)
                        gu2i = (u2_i-ub2_i)*ave.wdg[ib,jb]+ub2_i*wdgu3*(3/5)
                        gu3i = (u3_i-ub3_i)*ave.wdg[ib,jb]+ub3_i*wdgu33*(3/5)
                        gu4i = (u4_i-ub4_i)*ave.wdg[ib,jb]+ub4_i*wdgu33*(3/5)
                        gu5i = (u5_i-ub5_i)*ave.wdg[ib,jb]+ub5_i*wdgu3*(3/5)
                        gu6i = (u6_i-ub6_i)*ave.wdg[ib,jb]+ub6_i*wdgu3*(3/5)
                        gu7i = (u7_i-ub7_i)*ave.wdg[ib,jb]+ub7_i*wdgu3*(3/5)
                        gu8i = (u8_i-ub8_i)*ave.wdg[ib,jb]+ub8_i*wdgu33*(3/5)
                        gu9i = (u9_i-ub9_i)*ave.wdg[ib,jb]+ub9_i*wdgu33*(3/5)
                        gu10i = (u10_i-ub10_i)*ave.wdg[ib,jb]+ub10_i*wdgu33*(3/5)
                        gu1igt1 = (u1_i-ub1_i)*aveWG.wdgt1[is,ib,jb]+ub1_i*wdgu3gt1*(3/5)
                        gu2igt3 = (u2_i-ub2_i)*aveWG.wdgt3[is,ib,jb]+ub2_i*wdgu3gt3*(3/5)
                        gu3igt1 = (u3_i-ub3_i)*aveWG.wdgt1[is,ib,jb]+ub3_i*wdgu33gt1*(3/5)
                        gu4igt3 = (u4_i-ub4_i)*aveWG.wdgt3[is,ib,jb]+ub4_i*wdgu33gt3*(3/5)
                        gu6igu1 = (u6_i*wdgu1)
                        gu7igu3 = (u7_i*wdgu3)
                        gu9igu1 = (u9_i*wdgu1)
                        gu10igu3= (u10_i*wdgu3)

                        kpar_gnp0 = k_par0*aveKG.kpargnp0[is,ib,jb]
                        kpar_gp1p0 = k_par0*aveKG.kpargp1p0[is,ib,jb]
                        kpar_gp3p0 = k_par0*aveKG.kpargp3p0[is,ib,jb]
                        if(use_bper_in)
                            kpar_gp1b0 = k_par0*aveKG.kpargp1b0[is,ib,jb]
                            kpar_gr11b0 = k_par0*aveKG.kpargr11b0[is,ib,jb]
                            kpar_gr13b0 = k_par0*aveKG.kpargr13b0[is,ib,jb]
                            if(vpar_model_in==0)
                                kpar_gnbp = k_par0*aveKG.kpargnbp[is,ib,jb]
                                kpar_gp1bp = k_par0*aveKG.kpargp1bp[is,ib,jb]
                                kpar_gp3bp = k_par0*aveKG.kpargp3bp[is,ib,jb]
                                kpar_gr11bp = k_par0*aveKG.kpargr11bp[is,ib,jb]
                                kpar_gr13bp = k_par0*aveKG.kpargr13bp[is,ib,jb]
                            end
                        end
                        kpar_gu1 = aveKG.kpargu1[is,ib,jb]
                        kpar_gu3 = aveKG.kpargu3[is,ib,jb]
                        kpar_gb1 = ft2*b1*ave.kpar_eff[is,ib,jb]
                        kpar_gb3 = ft2*b3*ave.kpar_eff[is,ib,jb]
                        kpar_gb33 = b33*ave.kpar_eff[is,ib,jb]
                        kpar_gb1gt1 = ft2*b1*aveKG.kpargt1[is,ib,jb]
                        kpar_gb3gt3 = ft2*b3*aveKG.kpargt3[is,ib,jb]
                        kpar_gb33gt1 = b33*aveKG.kpargt1[is,ib,jb]
                        modkpar_gd1 = ft*d1*ave.modkpar_eff[is,ib,jb]
                        modkpar_gd3 = ft*d3*ave.modkpar_eff[is,ib,jb]
                        modkpar_gd33 = (d33/ft)*ave.modkpar_eff[is,ib,jb]
                        modkpar_gd1gu1 = ft*d1*aveKG.modkpargu1[is,ib,jb]
                        modkpar_gd3gu3 = ft*d3*aveKG.modkpargu3[is,ib,jb]
                        modkpar_gd33gu1 = (d33/ft)*aveKG.modkpargu1[is,ib,jb]

                        grad_gu1 = 0.0
                        grad_gu3 = 0.0
                        dgr13 = -b3*(ft2*aveKG.kpargt3[is,ib,jb] - aveKG.kpargt1[is,ib,jb]/3
                                    - ft2*aveKG.kpargu3[is,ib,jb] + aveKG.kpargu1[is,ib,jb]/3)

                        if(nbasis==1 || linsker==0.0)

                            gradgp1=0.0
                            gradgr11=0.0
                            gradgr13=0.0
                            gradgp1p1=0.0
                            gradgr11p1=0.0
                            gradgr13p1=0.0
                        else

                            gradgp1  = linsker*aveGrad.gradgp1p0[is,ib,jb]
                            gradgr11 = linsker*aveGrad.gradgr11p0[is,ib,jb]
                            gradgr13 = linsker*aveGrad.gradgr13p0[is,ib,jb]
                            gradgp1p1 = linsker*aveGrad.gradgp1p1[is,ib,jb]
                            gradgr11p1 = linsker*aveGrad.gradgr11p1[is,ib,jb]
                            gradgr13p1 = linsker*aveGrad.gradgr13p1[is,ib,jb]
                        end
                    end  # nroot>6

                    w_d1 = -ghat_in*w_d0
                    k_par1 = k_par0
                    if(js!=is)
                        w_d1 = 0.0
                        k_par1 = 0.0
                    end
                    modw_d1 = abs(w_d1)
                    modk_par1 = abs(k_par1)
                    w_dh= w_d1*ave.wdh[ib,jb]
                    w_dg= w_d1*ave.wdg[ib,jb]
                    k_par = k_par1*ave.kpar_eff[is,ib,jb]
                    k_par_psi = k_par0*ave.kpar_eff[is,ib,jb]
                    gradB1=k_par1*gradB_factor_in
                    if(nbasis==1 || gradB1==0.0)
                        gradB=0.0
                        gradBhp1=0.0
                        gradBhp3=0.0
                        gradBhr11=0.0
                        gradBhr13=0.0
                        gradBhr33=0.0
                        gradBhu1=0.0
                        gradBhu3=0.0
                        gradBhu33=0.0
                        if(nroot>6)
                            gradBgp1=0.0
                            gradBgp3=0.0
                            gradBgr11=0.0
                            gradBgr13=0.0
                            gradBgr33=0.0
                            gradBgu1=0.0
                            gradBgu3=0.0
                            gradBgu33=0.0
                        end
                    else
                        @warn "NOT TESTED eigensolve.jl ln 950"
                        gradB = gradB1*aveGradB.gradB[ib,jb]
                        gradBhp1=gradB1*aveGradB.gradBhp1[is,ib,jb]
                        gradBhp3=gradB1*aveGradB.gradBhp3[is,ib,jb]
                        gradBhr11=gradB1*aveGradB.gradBhr11[is,ib,jb]
                        gradBhr13=gradB1*aveGradB.gradBhr13[is,ib,jb]
                        gradBhr33=gradB1*aveGradB.gradBhr33[is,ib,jb]
                        gradBhu1=gradB1*aveGradB.gradBhu1[is,ib,jb]
                        gradBhu3=gradB1*aveGradB.gradBhu3[is,ib,jb]
                        gradBhu33=gradB1*aveGradB.gradBhu33[is,ib,jb]
                        if(nroot>6)
                            gradBgp1=gradB1*aveGradB.gradBgp1[is,ib,jb]
                            gradBgp3=gradB1*aveGradB.gradBgp3[is,ib,jb]
                            gradBgr11=gradB1*aveGradB.gradBgr11[is,ib,jb]
                            gradBgr13=gradB1*aveGradB.gradBgr13[is,ib,jb]
                            gradBgr33=gradB1*aveGradB.gradBgr33[is,ib,jb]
                            gradBgu1=gradB1*aveGradB.gradBgu1[is,ib,jb]
                            gradBgu3=gradB1*aveGradB.gradBgu3[is,ib,jb]
                            gradBgu33=gradB1*aveGradB.gradBgu33[is,ib,jb]
                        end
                    end

                    M_i = zsIS*vsIS/tausIS
                    J_j = zsJS*asJS*vsJS
                    E_i = zsIS/tausIS
                    N_j = zsJS*asJS




#*************************************************************
# matrix in order n, u, p1, p3, q1,q3
#*************************************************************





#*************************************************************
# n_u equ #1
#*************************************************************
                    ia0 = (is-ns0)*nroot*nbasis
                    ja0 = (js-ns0)*nroot*nbasis

                    ia = ib + ia0
                    #****************
                    #  untrapped terms
                    #****************
                    ja = jb + ja0

                    phi_A = N_j*im*w_s*(rlnsIS*hn + rltsIS*1.5*(hp3-hn))
                    if(vpar_model_in==0)
                        phi_A = phi_A + N_j*E_i*kpar_hnp0*vpar[is]
                    end
                    phi_B = -hn*E_i*N_j
                    sig_A = zero(K)
                    sig_B = zero(K)
                    psi_A = zero(K)
                    psi_B = zero(K)
                    phi_AU = zero(K)
                    phi_BU = zero(K)
                    psi_AN = zero(K)
                    psi_BN = zero(K)
                    if(use_bpar_in)
                        sig_A = -(betae_sig*(asJS*tausJS*zsIS/massIS) *
                            (im*w_s*(rlnsIS*h10n + rltsIS*1.5*(h10p3-h10n))))
                        sig_B = (betae_sig*h10n*asJS*tausJS*zsIS*zsIS
                            /(tausIS*massIS))
                        sig_A = sig_A - damp_sig*sig_B
                    end
                    if(use_bper_in)
                        psi_A = -betae_psi*J_j*im*w_s*vpar_shear[is]*hp1b0
                        if(vpar_model_in==0)
                            psi_A = psi_A -betae_psi*J_j*im*w_cd*vpar[is]*wdhp1b0
                            psi_B = betae_psi*M_i*J_j*vpar[is]*hp1b0/vsIS
                            phi_AU = betae_psi*U0*J_j*(im*w_s*(rlnsIS*hnbp + rltsIS*1.5*(hp3bp-hnbp))
                                + E_i*kpar_hnbp*vpar[is])
                            phi_BU = -betae_psi*U0*E_i*J_j*vpar[is]*hp1bp
                            psi_AN = betae_psi*U0*N_j*im*(w_cd*vpar[is]*wdhp1bp + w_s*vpar_shear[is]*hp1bp)
                            psi_BN = -betae_psi*U0*M_i*N_j*vpar[is]*hp1bp/vsIS
                        end
                    end
                    amat[ia,ja] = phi_A + psi_AN
                    bmat[ia,ja] = d_ab + phi_B + psi_BN

                    ja = nbasis+jb + ja0
                    amat[ia,ja] = psi_A + phi_AU -k_par*vsIS + am*gradB*vsIS
                    bmat[ia,ja] = psi_B + phi_BU

                    ja = 2*nbasis+jb + ja0
                    amat[ia,ja] = -0.5*sig_A -0.5*im*w_dh*tausIS/zsIS
                    bmat[ia,ja] = -0.5*sig_B

                    ja = 3*nbasis+jb + ja0
                    amat[ia,ja] = 1.5*sig_A -1.5*im*w_dh*tausIS/zsIS
                    bmat[ia,ja] = 1.5*sig_B

                    ja = 4*nbasis+jb + ja0
                    amat[ia,ja] = 0.0
                    bmat[ia,ja] = 0.0

                    ja = 5*nbasis+jb + ja0
                    amat[ia,ja] = 0.0
                    bmat[ia,ja] = 0.0


                    if(nroot>6)
                        #****************
                        #  n_u ghost terms
                        #****************
                        ja = 6*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(phi_A + psi_AN)
                        bmat[ia,ja] = -1.0*(phi_B + psi_BN)

                        ja = 7*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(psi_A + phi_AU)
                        bmat[ia,ja] = -1.0*(psi_B + phi_BU)

                        ja = 8*nbasis+jb + ja0
                        amat[ia,ja] = -0.5*(-1.0*sig_A)
                        bmat[ia,ja] = -0.5*(-1.0*sig_B)

                        ja = 9*nbasis+jb + ja0
                        amat[ia,ja] = 1.5*(-1.0*sig_A)
                        bmat[ia,ja] = 1.5*(-1.0*sig_B)

                        ja = 10*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 11*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        #****************
                        #  n_u trapped particle terms
                        #****************
                        ja = 12*nbasis+jb + ja0
                        amat[ia,ja] = phi_A + psi_AN
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = 13*nbasis+jb + ja0
                        amat[ia,ja] = -0.5*sig_A
                        bmat[ia,ja] = -0.5*sig_B

                        ja = 14*nbasis+jb + ja0
                        amat[ia,ja] = 1.5*sig_A
                        bmat[ia,ja] = 1.5*sig_B

                    end #  nroot>6

#*************************************************************
# u_par_u equ #2
#*************************************************************
                    ia = nbasis+ib + ia0

                    phi_A = N_j*im*w_s*vpar_shear[is]*hp1/vsIS
                    phi_B = zero(K)
                    if(vpar_model_in==0)
                        phi_A = phi_A  +N_j*im*w_cd*wdhp1p0*vpar[is]/vsIS
                            + d_1*(nuei_u_u_1+nuei_u_q3_1*5.0/3.0)*hp1*E_i*N_j*vpar[is]/vsIS
                        phi_B = -E_i*N_j*hp1*vpar[is]/vsIS
                    end
                    sig_A = zero(K)
                    sig_B = zero(K)
                    psi_A = zero(K)
                    psi_B = zero(K)
                    phi_AU = zero(K)
                    phi_BU = zero(K)
                    psi_AN = zero(K)
                    psi_BN = zero(K)
                    if(use_bper_in)
                        psi_A = -betae_psi*J_j*vsIS*im*w_s*(rlnsIS*hp1b0+1.5*rltsIS*(hr13b0-hp1b0))
                        psi_B = betae_psi*M_i*J_j*hp1b0
                        psi_A = psi_A - damp_psi*psi_B
                        if(vpar_model_in==0)
                            psi_A = psi_A  -betae_psi*J_j*M_i*kpar_hp1b0*vpar[is]
                            phi_AU = betae_psi*U0*J_j*im*(w_s*vpar_shear[is]*hp1bp +w_cd*wdhp1bp*vpar[is])/vsIS + d_1*(nuei_u_u_1+nuei_u_q3_1*5.0/3.0)*hp1*E_i*betae_psi*U0*J_j*vpar[is]/vsIS
                            phi_BU = -betae_psi*U0*E_i*J_j*hp1bp*vpar[is]/vsIS
                            psi_AN = betae_psi*U0*N_j*(vsIS*im*w_s*(rlnsIS*hp1bp+1.5*rltsIS*(hr13bp-hp1bp)) + M_i*kpar_hp1bp*vpar[is])
                            psi_BN = -betae_psi*U0*M_i*N_j*hp1bp
                        end
                    end

                    #****************
                    # u_par_u  untrapped terms
                    #****************
                    ja = jb + ja0
                    amat[ia,ja] = phi_A + psi_AN
                    bmat[ia,ja] = phi_B + psi_BN

                    ja = nbasis + jb + ja0
                    amat[ia,ja] =  (psi_A +phi_AU -d_ee*nuei_u_u_1
                        + xnuei*(d_ab*xnu_u_u_1 - d_ij*hu3*xnu_u_q3_1))
                    bmat[ia,ja] = d_ab + psi_B +phi_BU

                    ja = 2*nbasis+jb + ja0
                    amat[ia,ja] =  -(k_par - k_par1*gradhp1p1)*vsIS + (am+bm*0.5)*gradB*vsIS
                    bmat[ia,ja] = 0.0

                    ja = 3*nbasis+jb + ja0
                    amat[ia,ja] =  - bm*1.5*gradB*vsIS
                    bmat[ia,ja] = 0.0

                    ja = 4*nbasis+jb + ja0
                    amat[ia,ja] = -0.5*im*w_dh*tausIS/zsIS
                    bmat[ia,ja] = 0.0

                    ja = 5*nbasis+jb + ja0
                    amat[ia,ja] = -1.5*im*w_dh*tausIS/zsIS -d_ee*nuei_u_q3_1 +xnuei*d_ab*xnu_u_q3_1
                    bmat[ia,ja] = 0.0

                    if(nroot>6)
                        #****************
                        #  u_par_u ghost terms
                        #****************
                        ja = 6*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(phi_A + psi_AN)
                        bmat[ia,ja] = -1.0*(phi_B + psi_BN)

                        ja = 7*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(psi_A + phi_AU)
                        bmat[ia,ja] = -1.0*(psi_B + phi_BU)

                        ja = 8*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 9*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 10*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 11*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        #****************
                        # u_par_u trapped particle terms
                        #****************
                        ja = 12*nbasis+jb + ja0
                        amat[ia,ja] = phi_A + psi_AN
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = 13*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 14*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
                    end  #  nroot>6

#*************************************************************
# p_par_u equ #3
#*************************************************************
                    ia = 2*nbasis+ib + ia0

                    phi_A = N_j*im*w_s*(rlnsIS*hp1 + rltsIS*1.5*(hr13-hp1))
                    if(vpar_model_in==0)
                        phi_A = phi_A +N_j*E_i*kpar_hp1p0*vpar[is]
                    end
                    phi_B = -hp1*E_i*N_j
                    sig_A = zero(K)
                    sig_B = zero(K)
                    psi_A = zero(K)
                    psi_B = zero(K)
                    phi_AU = zero(K)
                    phi_BU = zero(K)
                    psi_AN = zero(K)
                    psi_BN = zero(K)
                    if(use_bpar_in)
                        sig_A = -betae_sig*(asJS*tausJS*zsIS/massIS)*(im*w_s*(rlnsIS*h10p1 + rltsIS*1.5*(h10r13-h10p1)))
                        sig_B = betae_sig*h10p1*asJS*tausJS*zsIS*zsIS/(tausIS*massIS)
                        sig_A = sig_A - damp_sig*sig_B
                    end
                    if(use_bper_in)
                        psi_A = -betae_psi*J_j*im*w_s*vpar_shear[is]*hr11b0
                        if(vpar_model_in==0)
                            psi_A = psi_A -betae_psi*J_j*im*w_cd*vpar[is]*wdhr11b0
                            psi_B = betae_psi*M_i*J_j*vpar[is]*hr11b0/vsIS
                            phi_AU = betae_psi*U0*J_j*(im*w_s*(rlnsIS*hp1bp + rltsIS*1.5*(hr13bp-hp1bp))
                                    +E_i*kpar_hp1bp*vpar[is])
                            phi_BU = -betae_psi*U0*hp1bp*E_i*J_j
                            psi_AN = betae_psi*U0*N_j*im*(w_cd*vpar[is]*wdhr11bp+w_s*vpar_shear[is]*hr11bp)
                            psi_BN = -betae_psi*U0*M_i*N_j*vpar[is]*hr11bp/vsIS
                        end
                    end
                    #****************
                    #  p_par_u  untrapped terms
                    #****************
                    ja = jb + ja0
                    amat[ia,ja] = (phi_A + psi_AN
                            + 2*tausIS*(modw_d1*hv1rht1/abs(zsIS) + w_d1*im*hv1iht1/zsIS)
                            + 2*tausIS*(modw_d1*hv2rht3/abs(zsIS) +w_d1*im*hv2iht3/zsIS)
                            -xnuei*d_ij*xnu_p1_1*(ht1-ht3))
                    bmat[ia,ja] = phi_B + psi_BN

                    ja = nbasis+jb + ja0
                    amat[ia,ja] = k_par1*grad_hu1*vsIS + psi_A + phi_AU
                    bmat[ia,ja] = psi_B + phi_BU

                    ja = 2*nbasis+jb + ja0
                    amat[ia,ja] =  (-0.5*sig_A
                            -im*w_d1*(tausIS/zsIS)*(0.5*wdhu1+1.5*wdhu3)
                            -2.0*tausIS*(modw_d1*hv1r/abs(zsIS) +w_d1*im*hv1i/zsIS)
                            -d_ee*nuei_p1_p1_1
                            +xnuei*d_ab*xnu_p1_1)
                    bmat[ia,ja] = d_ab -0.5*sig_B

                    ja = 3*nbasis+jb + ja0
                    amat[ia,ja] = (1.5*sig_A -2.0*tausIS*
                            (modw_d1*hv2r/abs(zsIS) +w_d1*im*hv2i/zsIS)
                            -d_ee*nuei_p1_p3_1 -xnuei*d_ab*xnu_p1_1)
                    bmat[ia,ja] = 1.5*sig_B

                    ja = 4*nbasis+jb + ja0
                    amat[ia,ja] = -k_par*vsIS + (am +bm)*gradB*vsIS
                    bmat[ia,ja] = 0.0

                    ja = 5*nbasis+jb + ja0
                    amat[ia,ja] = -bm*3.0*gradB*vsIS
                    bmat[ia,ja] = 0.0

                    if(nroot>6)
                        #****************
                        #   p_par_u ghost terms
                        #****************
                        ja = 6*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(phi_A + psi_AN)
                        bmat[ia,ja] = -1.0*(phi_B + psi_BN)

                        ja = 7*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(psi_A + phi_AU)
                        bmat[ia,ja] = -1.0*(psi_B + phi_BU)

                        ja = 8*nbasis+jb + ja0
                        amat[ia,ja] = -0.5*(-1.0*sig_A)
                        bmat[ia,ja] = -0.5*(-1.0*sig_B)

                        ja = 9*nbasis+jb + ja0
                        amat[ia,ja] = 1.5*(-1.0*sig_A)
                        bmat[ia,ja] = 1.5*(-1.0*sig_B)

                        ja = 10*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 11*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
                        #****************
                        #  p_par_u trapped particle terms
                        #****************
                        ja = 12*nbasis+jb + ja0
                        amat[ia,ja] = phi_A + psi_AN
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = 13*nbasis+jb + ja0
                        amat[ia,ja] = -0.5*sig_A
                        bmat[ia,ja] = -0.5*sig_B

                        ja = 14*nbasis+jb + ja0
                        amat[ia,ja] = 1.5*sig_A
                        bmat[ia,ja] = 1.5*sig_B

                    end   # nroot >6
#*************************************************************
# p_tot_u equ #4
#*************************************************************
                    ia = 3*nbasis+ib + ia0

                    phi_A = N_j*im*w_s*(rlnsIS*hp3+rltsIS*1.5*(hr33-hp3))
                    if(vpar_model_in==0)
                        phi_A = phi_A + N_j*E_i*kpar_hp3p0*vpar[is]
                    end
                    phi_B = -hp3*E_i*N_j
                    sig_A = zero(K)
                    sig_B = zero(K)
                    psi_A = zero(K)
                    psi_B = zero(K)
                    phi_AU = zero(K)
                    phi_BU = zero(K)
                    psi_AN = zero(K)
                    psi_BN = zero(K)
                    if(use_bpar_in)
                        sig_A = -betae_sig*(asJS*tausJS*zsIS/massIS) * (im*w_s*(rlnsIS*h10p3 + rltsIS*1.5*(h10r33-h10p3)))
                        sig_B = betae_sig*h10p3*asJS*tausJS*zsIS*zsIS /(tausIS*massIS)
                        sig_A = sig_A - damp_sig*sig_B
                    end
                    if(use_bper_in)
                        psi_A = -betae_psi*J_j*im*w_s*vpar_shear[is]*hr13b0
                        if(vpar_model_in==0)
                            psi_A = psi_A -betae_psi*J_j*im*w_cd*vpar[is]*wdhr13b0
                            psi_B = betae_psi*M_i*J_j*vpar[is]*hr13b0/vsIS
                            phi_AU = betae_psi*U0*J_j*(im*w_s*(rlnsIS*hp3bp+rltsIS*1.5*(hr33bp-hp3bp))
                                    + E_i*kpar_hp3bp*vpar[is])
                            phi_BU = -betae_psi*U0*hp3bp*E_i*J_j
                            psi_AN = betae_psi*U0*N_j*im*(w_cd*vpar[is]*wdhr13bp +w_s*vpar_shear[is]*hr13bp)
                            psi_BN = -betae_psi*U0*M_i*N_j*vpar[is]*hr13bp/vsIS
                        end
                    end
                    #****************
                    #   p_tot_u untrapped terms
                    #****************
                    ja = jb + ja0
                    amat[ia,ja] = (phi_A + psi_AN
                            +2.0*tausIS*(modw_d1*hv3rht1/abs(zsIS) +w_d1*im*hv3iht1/zsIS)
                            +2.0*tausIS*(modw_d1*hv4rht3/abs(zsIS) +w_d1*im*hv4iht3/zsIS))
                    bmat[ia,ja] = phi_B + psi_BN

                    ja = nbasis+jb + ja0
                    amat[ia,ja] = k_par1*grad_hu3*vsIS + psi_A + phi_AU
                    bmat[ia,ja] = psi_B + phi_BU

                    ja = 2*nbasis+jb + ja0
                    amat[ia,ja] = (-0.5*sig_A
                            -im*w_d1*(tausIS/zsIS)*0.5*wdhu3
                            -2.0*tausIS*(modw_d1*hv3r/abs(zsIS) +w_d1*im*hv3i/zsIS))
                    bmat[ia,ja] = -0.5*sig_B

                    ja = 3*nbasis+jb + ja0
                    amat[ia,ja] =  (1.5*sig_A
                            -im*w_d1*(tausIS/zsIS)*1.5*wdhu33
                            -2.0*tausIS*(modw_d1*hv4r/abs(zsIS) +w_d1*im*hv4i/zsIS))
                    bmat[ia,ja] = d_ab + 1.5*sig_B

                    ja = 4*nbasis+jb + ja0
                    amat[ia,ja] = 0.0
                    bmat[ia,ja] = 0.0

                    ja = 5*nbasis+jb + ja0
                    amat[ia,ja] = -k_par*vsIS + am*gradB*vsIS
                    bmat[ia,ja] = 0.0

                    if(nroot>6)
                        #****************
                        #   p_tot_u ghost terms
                        #****************
                        ja = 6*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(phi_A + psi_AN)
                        bmat[ia,ja] = -1.0*(phi_B + psi_BN)

                        ja = 7*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(psi_A + phi_AU)
                        bmat[ia,ja] = -1.0*(psi_B + phi_BU)

                        ja = 8*nbasis+jb + ja0
                        amat[ia,ja] = -0.5*(-1.0*sig_A)
                        bmat[ia,ja] = -0.5*(-1.0*sig_B)

                        ja = 9*nbasis+jb + ja0
                        amat[ia,ja] = 1.5*(-1.0*sig_A)
                        bmat[ia,ja] = 1.5*(-1.0*sig_B)

                        ja = 10*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 11*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        #****************
                        #   p_tot_u trapped terms
                        #****************
                        ja = 12*nbasis+jb + ja0
                        amat[ia,ja] = phi_A + psi_AN
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = 13*nbasis+jb + ja0
                        amat[ia,ja] = -0.5*sig_A
                        bmat[ia,ja] = -0.5*sig_B

                        ja = 14*nbasis+jb + ja0
                        amat[ia,ja] = 1.5*sig_A
                        bmat[ia,ja] = 1.5*sig_B

                    end # nroot>6
#*************************************************************
# q_par_u equ #5
#*************************************************************
                    ia = 4*nbasis+ib + ia0

                    phi_A = N_j*im*w_s*hr11*vpar_shear[is]/vsIS
                    phi_B = zero(K)
                    if(vpar_model_in==0)
                        phi_A = (phi_A +N_j*im*w_cd*wdhr11p0*vpar[is]/vsIS
                                + d_1*(nuei_q1_u_1 + (5/3)*nuei_q1_q3_1+3*nuei_q1_q1_1)
                                *hr11*E_i*N_j*vpar[is]/vsIS)
                        phi_B = -E_i*N_j*hr11*vpar[is]/vsIS
                    end
                    sig_A = zero(K)
                    sig_B = zero(K)
                    psi_A = zero(K)
                    psi_B = zero(K)
                    phi_AU = zero(K)
                    phi_BU = zero(K)
                    psi_AN = zero(K)
                    psi_BN = zero(K)
                    if(use_bper_in)
                        psi_A = -betae_psi*J_j*vsIS*im*w_s*(rlnsIS*hr11b0+1.5*rltsIS*(hw113b0-hr11b0))
                        psi_B =betae_psi*M_i*J_j*hr11b0
                        psi_A = psi_A - damp_psi*psi_B
                        if(vpar_model_in==0)
                            psi_A = psi_A  -betae_psi*J_j*M_i*kpar_hr11b0*vpar[is]
                            phi_AU = (betae_psi*U0*J_j*im*(w_s*hr11bp*vpar_shear[is] +w_cd*wdhr11bp*vpar[is])/vsIS
                                    + d_1*(nuei_q1_u_1 + (5.0/3.0)*nuei_q1_q3_1+3.0*nuei_q1_q1_1)
                                    *hr11bp*E_i*betae_psi*U0*J_j*vpar[is]/vsIS)
                            phi_BU = -betae_psi*U0*E_i*J_j*hr11bp*vpar[is]/vsIS
                            psi_AN = (betae_psi*U0*N_j*(vsIS*im*w_s*(rlnsIS*hr11bp+1.5*rltsIS*(hw113bp-hr11bp))
                                    +M_i*kpar_hr11bp*vpar[is]))
                            psi_BN =-betae_psi*U0*M_i*N_j*hr11bp
                        end
                    end
                    #****************
                    #  q_par_u untrapped terms
                    #****************
                    ja =  jb + ja0
                    amat[ia,ja] = k_par1*kpar_hb1ht1*vsIS + phi_A + psi_AN
                    bmat[ia,ja] = phi_B + psi_BN

                    ja =   nbasis+jb + ja0
                    amat[ia,ja] = (psi_A + phi_AU +modk_par1*modkpar_hd1hu1*vsIS
                            - tausIS*(modw_d1*hv5r/abs(zsIS) + w_d1*im*hv5i/zsIS)
                            -d_ee*nuei_q1_u_1
                            +xnuei*(d_ab*xnu_q1_u_1 - d_ij*hu1*xnu_q1_q1_1
                            - d_ij*hu3*xnu_q1_q3_1))
                    bmat[ia,ja] = psi_B + phi_BU

                    ja = 2*nbasis+jb + ja0
                    amat[ia,ja] = (-k_par1*(kpar_hu1 -grad_hu1 - gradhr11p1)*vsIS
                            - k_par1*kpar_hb1*vsIS + (am+bm*1.5)*gradBhu1*vsIS
                            - d_11*k_par*vsIS*c06)
                    bmat[ia,ja] = 0.0

                    ja = 3*nbasis+jb + ja0
                    amat[ia,ja] = (-bm*4.5*gradBhu3*vsIS
                            + d_11*k_par*vsIS*c06 - d_11*k_par*vsIS*c08)
                    bmat[ia,ja] = 0.0

                    ja = 4*nbasis+jb + ja0
                    amat[ia,ja] = (-modk_par1*modkpar_hd1*vsIS
                            - tausIS*(modw_d1*hv6r/abs(zsIS) + w_d1*im*hv6i/zsIS)
                            -d_ee*nuei_q1_q1_1
                            +xnuei*d_ab*xnu_q1_q1_1)
                    bmat[ia,ja] = d_ab

                    ja = 5*nbasis+jb + ja0
                    amat[ia,ja] = (-tausIS*(modw_d1*hv7r/abs(zsIS) + w_d1*im*hv7i/zsIS)
                            -d_ee*nuei_q1_q3_1
                            +xnuei*d_ab*xnu_q1_q3_1)
                    bmat[ia,ja] = 0.0

                    if(nroot>6)
                        #****************
                        #  q_par_u ghost terms
                        #****************
                        ja = 6*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(phi_A + psi_AN)
                        bmat[ia,ja] = -1.0*(phi_B + psi_BN)

                        ja = 7*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(psi_A + phi_AU)
                        bmat[ia,ja] = -1.0*(psi_B + phi_BU)

                        ja = 8*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 9*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 10*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 11*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
                        #****************
                        #  q_par_u trapped terms
                        #****************
                        ja = 12*nbasis+jb + ja0
                        amat[ia,ja] = phi_A + psi_AN
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = 13*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 14*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                    end   # nroot >6
#*************************************************************
# q_tot_u equ #6
#*************************************************************
                    ia = 5*nbasis+ib + ia0

                    phi_A = N_j*im*w_s*vpar_shear[is]*hr13/vsIS
                    phi_B = zero(K)
                    if(vpar_model_in==0)
                        phi_A = (phi_A +N_j*im*w_cd*wdhr13p0*vpar[is]/vsIS
                                + d_1*(nuei_q3_u_1 + (5.0/3.0)*nuei_q3_q3_1)
                                *hr13*E_i*N_j*vpar[is]/vsIS)
                        phi_B = -E_i*N_j*hr13*vpar[is]/vsIS
                    end
                    sig_A = zero(K)
                    sig_B = zero(K)
                    psi_A = zero(K)
                    psi_B = zero(K)
                    phi_AU = zero(K)
                    phi_BU = zero(K)
                    psi_AN = zero(K)
                    psi_BN = zero(K)
                    if(use_bper_in)
                        psi_A = -betae_psi*J_j*vsIS*im*w_s*(rlnsIS*hr13b0+1.5*rltsIS*(hw133b0-hr13b0))
                        psi_B = hr13b0*betae_psi*M_i*J_j
                        psi_A = psi_A - damp_psi*psi_B
                        if(vpar_model_in==0)
                            psi_A = psi_A -betae_psi*J_j*M_i*kpar_hr13b0*vpar[is]
                            phi_AU = (betae_psi*U0*J_j*im*(w_s*vpar_shear[is]*hr13bp +w_cd*wdhr13bp*vpar[is])/vsIS
                                    + d_1*(nuei_q3_u_1 + (5.0/3.0)*nuei_q3_q3_1)
                                    *hr13bp*E_i*betae_psi*U0*J_j*vpar[is]/vsIS)
                            phi_BU = -betae_psi*U0*E_i*J_j*hr13bp*vpar[is]/vsIS
                            psi_AN = betae_psi*U0*N_j*(vsIS*im*w_s*(rlnsIS*hr13bp+1.5*rltsIS*(hw133bp-hr13bp))
                                    + M_i*kpar_hr13bp*vpar[is])
                            psi_BN = -betae_psi*U0*hr13bp*M_i*N_j
                        end
                    end
                    #****************
                    #  q_tot_u untrapped terms
                    #****************
                    ja = jb + ja0
                    amat[ia,ja] = (phi_A + psi_AN
                            +k_par1*(kpar_hb3ht3 -dhr13+ kpar_hb33ht1)*vsIS)
                    bmat[ia,ja] = phi_B + psi_BN

                    ja = nbasis+jb + ja0
                    amat[ia,ja] =  (psi_A  + phi_AU +
                            modk_par1*(modkpar_hd3hu3 + modkpar_hd33hu1)*vsIS
                            - tausIS*(modw_d1*hv8r/abs(zsIS) + w_d1*im*hv8i/zsIS)
                            -d_ee*nuei_q3_u_1
                            +xnuei*(d_ab*xnu_q3_u_1 - d_ij*hu3*xnu_q3_q3_1))
                    bmat[ia,ja] = psi_B + phi_BU

                    ja = 2*nbasis+jb + ja0
                    amat[ia,ja] = (- k_par1*(kpar_hu3 -grad_hu3 - gradhr13p1)*vsIS
                            - k_par1*kpar_hb33*vsIS + (am+bm*0.5)*gradBhu3*vsIS
                            -d_11*k_par*vsIS*c07)
                    bmat[ia,ja] = 0.0

                    ja = 3*nbasis+jb + ja0
                    amat[ia,ja] = (-k_par1*kpar_hb3*vsIS - bm*1.5*gradBhu33*vsIS
                            +d_11*k_par*vsIS*c07 )
                    bmat[ia,ja] = 0.0

                    ja = 4*nbasis+jb + ja0
                    amat[ia,ja] = (-modk_par1*modkpar_hd33*vsIS
                            - tausIS*(modw_d1*hv9r/abs(zsIS) + w_d1*im*hv9i/zsIS))
                    bmat[ia,ja] = 0.0

                    ja = 5*nbasis+jb + ja0
                    amat[ia,ja] = (-modk_par1*modkpar_hd3*vsIS
                            - tausIS*(modw_d1*hv10r/abs(zsIS) + w_d1*im*hv10i/zsIS)
                            -d_ee*nuei_q3_q3_1
                            +xnuei*d_ab*xnu_q3_q3_1)
                    bmat[ia,ja] = d_ab

                    if(nroot>6)
                        #****************
                        #  q_tot_u ghost terms
                        #****************
                        ja = 6*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(phi_A + psi_AN)
                        bmat[ia,ja] = -1.0*(phi_B + psi_BN)

                        ja = 7*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(psi_A + phi_AU)
                        bmat[ia,ja] = -1.0*(psi_B + phi_BU)

                        ja = 8*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 9*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 10*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 11*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
                        #****************
                        #  q_tot_u trapped terms
                        #****************
                        ja = 12*nbasis+jb + ja0
                        amat[ia,ja] = phi_A + psi_AN
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = 13*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 14*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                    end    # nroot>6

                    if(nroot>6)
#*************************************************************
# n_g equ #7
#*************************************************************
                        ia = 6*nbasis+ib + ia0

                        phi_A = N_j*im*w_s*(rlnsIS*gn + rltsIS*1.5*(gp3-gn))
                        if(vpar_model_in==0)
                            phi_A = phi_A + N_j*E_i*kpar_gnp0*vpar[is]
                        end
                        phi_B = -E_i*N_j*gn
                        phi_A = phi_A +xnu_phi_b*xnuei*xnu_n_b*phi_B
                        sig_A = zero(K)
                        sig_B = zero(K)
                        psi_A = zero(K)
                        psi_B = zero(K)
                        phi_AU = zero(K)
                        phi_BU = zero(K)
                        psi_AN = zero(K)
                        psi_BN = zero(K)
                        if(use_bpar_in)
                            sig_A = (-betae_sig*(asJS*tausJS*zsIS/massIS)*
                                    (im*w_s*(rlnsIS*g10n + rltsIS*1.5*(g10p3-g10n))))
                            sig_B = (betae_sig*g10n*asJS*tausJS*zsIS*zsIS
                                    /(tausIS*massIS))
                            sig_A = sig_A - damp_sig*sig_B
                        end
                        if(use_bper_in)
                            psi_A = -betae_psi*J_j*im*w_s*vpar_shear[is]*gp1b0
                            if(vpar_model_in==0)
                                psi_A = psi_A -betae_psi*J_j*im*w_cd*vpar[is]*wdgp1b0
                                psi_B = betae_psi*M_i*J_j*vpar[is]*gp1b0/vsIS
                                phi_AU = betae_psi*U0*J_j*(im*w_s*(rlnsIS*gnbp + rltsIS*1.5*(gp3bp-gnbp))
                                        + E_i*kpar_gnbp*vpar[is])
                                phi_BU = -betae_psi*U0*E_i*J_j*gnbp
                                psi_AN = betae_psi*U0*N_j*im*(w_cd*vpar[is]*wdgp1bp+w_s*vpar_shear[is]*gp1bp)
                                psi_BN = -betae_psi*U0*M_i*N_j*vpar[is]*gp1bp/vsIS
                            end
                        end
                        #****************
                        #  n_g untrapped terms
                        #****************
                        ja = jb + ja0
                        amat[ia,ja] = phi_A + psi_AN +d_ee*nuei_n_n*bn
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = nbasis+jb + ja0
                        amat[ia,ja] = psi_A + phi_AU
                        bmat[ia,ja] = psi_B + phi_BU

                        ja = 2*nbasis+jb + ja0
                        amat[ia,ja] = (-0.5*sig_A +d_ee*nuei_n_p1*bp1
                                -d_ee*(1.0 - ft2)*(bn*an*nuei_n_n + bp3*ap3*nuei_n_p3 + bp1*ap1*nuei_n_p1))
                        bmat[ia,ja] = -0.5*sig_B

                        ja = 3*nbasis+jb + ja0
                        amat[ia,ja] = (1.5*sig_A +d_ee*nuei_n_p3*bp3
                                +d_ee*(1.0 - ft2)*(bn*an*nuei_n_n + bp3*ap3*nuei_n_p3 + bp1*ap1*nuei_n_p1))
                        bmat[ia,ja] = 1.5*sig_B

                        ja = 4*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 5*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
                        #****************
                        #  n_g ghost terms
                        #****************
                        ja = 6*nbasis+jb + ja0
                        amat[ia,ja] = (-1.0*(phi_A +psi_AN)
                                -d_ee*nuei_n_n +xnuei*d_ab*xnu_n_b )
                        bmat[ia,ja] = d_ab - (phi_B + psi_BN)

                        ja = 7*nbasis+jb + ja0
                        amat[ia,ja] = -k_par*vsIS + am*gradB*vsIS - (psi_A + phi_AU)
                        bmat[ia,ja] = -(psi_B + phi_BU)

                        ja = 8*nbasis+jb + ja0
                        amat[ia,ja] = (-0.5*(-1.0*sig_A)
                                -0.5*im*w_dg*tausIS/zsIS-d_ee*nuei_n_p1)
                        bmat[ia,ja] = -0.5*(-1*sig_B)

                        ja = 9*nbasis+jb + ja0
                        amat[ia,ja] = (1.5*(-1.0*sig_A)
                                -1.5*im*w_dg*tausIS/zsIS -d_ee*nuei_n_p3)
                        bmat[ia,ja] = 1.5*(-1.0*sig_B)

                        ja = 10*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 11*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
                        #****************
                        # n_g trapped terms
                        #****************
                        ja = 12*nbasis+jb + ja0
                        amat[ia,ja] = phi_A + psi_AN
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = 13*nbasis+jb + ja0
                        amat[ia,ja] = -0.5*sig_A
                        bmat[ia,ja] = -0.5*sig_B

                        ja = 14*nbasis+jb + ja0
                        amat[ia,ja] = 1.5*sig_A
                        bmat[ia,ja] = 1.5*sig_B
#*************************************************************
# u_par_g equ #8
#*************************************************************
                        ia = 7*nbasis+ib + ia0

                        phi_A = N_j*im*w_s*vpar_shear[is]*gp1/vsIS
                        if(vpar_model_in==0)
                            phi_A = phi_A  + N_j*im*w_cd*wdgp1p0*vpar[is]/vsIS
                            phi_B = -E_i*N_j*gp1*vpar[is]/vsIS
                        end
                        sig_A = zero(K)
                        sig_B = zero(K)
                        psi_A = zero(K)
                        psi_B = zero(K)
                        phi_AU = zero(K)
                        phi_BU = zero(K)
                        psi_AN = zero(K)
                        psi_BN = zero(K)
                        if(use_bper_in)
                            psi_A = -betae_psi*J_j*vsIS*im*w_s*(rlnsIS*gp1b0+1.5*rltsIS*(gr13b0-gp1b0))
                            psi_B =betae_psi*M_i*J_j*gp1b0
                            psi_A = psi_A - damp_psi*psi_B
                            if(vpar_model_in==0)
                                psi_A = psi_A  -betae_psi*J_j*M_i*kpar_gp1b0*vpar[is]
                                phi_AU = betae_psi*U0*J_j*im*(w_s*vpar_shear[is]*gp1bp +w_cd*wdgp1bp*vpar[is])/vsIS
                                phi_BU = -betae_psi*U0*E_i*J_j*gp1bp*vpar[is]/vsIS
                                psi_AN = betae_psi*U0*N_j*(vsIS*im*w_s*(rlnsIS*gp1bp+1.5*rltsIS*(gr13bp-gp1bp))
                                        + M_i*kpar_gp1bp*vpar[is])
                                psi_BN = -betae_psi*U0*M_i*N_j*gp1bp
                            end
                        end
                        #****************
                        # u_par_g untrapped terms
                        #****************
                        ja = jb + ja0
                        amat[ia,ja] = phi_A + psi_AN
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = nbasis+jb + ja0
                        amat[ia,ja] = psi_A + phi_AU + d_ee*ft3*nuei_u_u
                        bmat[ia,ja] = psi_B + phi_BU

                        ja = 2*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 3*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 4*nbasis+jb + ja0
                        amat[ia,ja] = (d_ee*ft5*nuei_u_q1
                                -d_ee*(1.0 -ft2)*(ft3*nuei_u_u*1.25 +ft3*nuei_u_q3*35.0/12.0 +ft5*nuei_u_q1*6.25) )

                        bmat[ia,ja] = 0.0

                        ja = 5*nbasis+jb + ja0
                        amat[ia,ja] = (d_ee*ft3*nuei_u_q3
                                +d_ee*(1.0-ft2)*(ft3*nuei_u_u*2.25 +ft3*nuei_u_q3*5.25 +ft5*nuei_u_q1*11.25))
                        bmat[ia,ja] = 0.0
                        #****************
                        # u_par_g ghost terms
                        #****************
                        ja = 6*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(phi_A + psi_AN)
                        bmat[ia,ja] = -1.0*(phi_B + psi_BN)

                        ja = 7*nbasis+jb + ja0
                        amat[ia,ja] = (-1.0*(psi_A + phi_AU)
                                -d_ee*nuei_u_u_t -d_ee*nuei_u_u
                                +xnuei*(d_ab*xnu_u_u_1 - d_ij*gu3*xnu_u_q3_1)
                                +xnuei*d_ab*xnu_u_b)
                        bmat[ia,ja] = d_ab -1.0*(psi_B + phi_BU)

                        ja = 8*nbasis+jb + ja0
                        amat[ia,ja] = (-(k_par - k_par1*gradgp1p1)*vsIS
                                +(am + bm*0.5)*gradB*vsIS )
                        bmat[ia,ja] = 0.0

                        ja = 9*nbasis+jb + ja0
                        amat[ia,ja] = -bm*1.5*gradB*vsIS
                        bmat[ia,ja] = 0.0

                        ja = 10*nbasis+jb + ja0
                        amat[ia,ja] = (-0.5*im*w_dg*tausIS/zsIS
                            -d_ee*nuei_u_q1)
                        bmat[ia,ja] = 0.0

                        ja = 11*nbasis+jb + ja0
                        amat[ia,ja] = (-1.5*im*w_dg*tausIS/zsIS
                                -d_ee*nuei_u_q3_t -d_ee*nuei_u_q3
                                +xnuei*d_ab*xnu_u_q3_1)
                        bmat[ia,ja] = 0.0
                        #****************
                        # u_par_g trapped terms
                        #****************
                        ja = 12*nbasis+jb + ja0
                        amat[ia,ja] = phi_A + psi_AN
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = 13*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 14*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
#*************************************************************
# p_par_g equ #9
#*************************************************************
                        ia = 8*nbasis+ib + ia0

                        phi_A = N_j*im*w_s*(rlnsIS*gp1 + rltsIS*1.5*(gr13-gp1))
                        if(vpar_model_in==0)
                            phi_A = phi_A  +N_j*E_i*kpar_gp1p0*vpar[is]
                        end
                        phi_B = -E_i*N_j*gp1
                        phi_A = phi_A +xnu_phi_b*xnuei*xnu_p1_b*phi_B
                        sig_A = zero(K)
                        sig_B = zero(K)
                        psi_A = zero(K)
                        psi_B = zero(K)
                        phi_AU = zero(K)
                        phi_BU = zero(K)
                        psi_AN = zero(K)
                        psi_BN = zero(K)
                        if(use_bpar_in)
                            sig_A =( -betae_sig*(asJS*tausJS*zsIS/massIS)*
                                    (im*w_s*(rlnsIS*g10p1 + rltsIS*1.5*(g10r13-g10p1))))
                            sig_B = (betae_sig*g10p1*asJS*tausJS*zsIS*zsIS
                                    /(tausIS*massIS))
                            sig_A = sig_A - damp_sig*sig_B
                        end
                        if(use_bper_in)
                            psi_A = -betae_psi*J_j*im*w_s*vpar_shear[is]*gr11b0
                            if(vpar_model_in==0)
                                psi_A = psi_A  -betae_psi*J_j*im*w_cd*vpar[is]*wdgr11b0
                                psi_B = betae_psi*M_i*J_j*vpar[is]*gr11b0/vsIS
                                phi_AU = betae_psi*U0*J_j*(im*w_s*(rlnsIS*gp1bp + rltsIS*1.5*(gr13bp-gp1bp))
                                        + E_i*kpar_gp1bp*vpar[is])
                                phi_BU = -betae_psi*U0*E_i*J_j*gp1bp
                                psi_AN = betae_psi*U0*N_j*im*(w_cd*vpar[is]*wdgr11bp+w_s*vpar_shear[is]*gr11bp)
                                psi_BN = -betae_psi*U0*M_i*N_j*vpar[is]*gr11bp/vsIS
                            end
                        end
                        #****************
                        # p_par_g untrapped terms
                        #****************
                        ja = jb + ja0
                        amat[ia,ja] = phi_A + psi_AN +d_ee*nuei_p1_n*bn
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = nbasis+jb + ja0
                        amat[ia,ja] = psi_A + phi_AU
                        bmat[ia,ja] = psi_B + phi_BU

                        ja = 2*nbasis+jb + ja0
                        amat[ia,ja] = (-0.5*sig_A +d_ee*nuei_p1_p1*bp1
                                -d_ee*(1.0 - ft2)*(bn*an*nuei_p1_n + bp3*ap3*nuei_p1_p3 + bp1*ap1*nuei_p1_p1))
                        bmat[ia,ja] = -0.5*sig_B

                        ja = 3*nbasis+jb + ja0
                        amat[ia,ja] = (1.5*sig_A +d_ee*nuei_p1_p3*bp3
                                +d_ee*(1.0 - ft2)*(bn*an*nuei_p1_n + bp3*ap3*nuei_p1_p3 + bp1*ap1*nuei_p1_p1))
                        bmat[ia,ja] = 1.5*sig_B

                        ja = 4*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 5*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
                        #****************
                        # p_par_g ghost terms
                        #****************

                        ja = 6*nbasis+jb + ja0
                        amat[ia,ja] = (-1.0*(phi_A + psi_AN)
                                +2.0*tausIS*(modw_d1*gu1rgt1/abs(zsIS) +w_d1*im*gu1igt1/zsIS)
                                +2.0*tausIS*(modw_d1*gu2rgt3/abs(zsIS) +w_d1*im*gu2igt3/zsIS)
                                -d_ee*nuei_p1_n
                                -xnuei*d_ij*xnu_p1_1*(gt1-ft2*gt3))
                        bmat[ia,ja] = -1.0*(phi_B + psi_BN)

                        ja = 7*nbasis+jb + ja0
                        amat[ia,ja] = k_par1*grad_gu1*vsIS -1.0*(psi_A + phi_AU)
                        bmat[ia,ja] = -1.0*(psi_B + phi_BU)

                        ja = 8*nbasis+jb + ja0
                        amat[ia,ja] = (-0.5*(-1.0*sig_A)
                                -im*w_d1*(tausIS/zsIS)*(0.5*wdgu1+1.5*wdgu3)
                                -2.0*tausIS*(modw_d1*gu1r/abs(zsIS) +w_d1*im*gu1i/zsIS)
                                -d_ee*nuei_p1_p1_t -d_ee*nuei_p1_p1
                                +xnuei*d_ab*(xnu_p1_1 + xnu_p1_b))
                        bmat[ia,ja] = d_ab -0.5*(-1.0*sig_B)

                        ja = 9*nbasis+jb + ja0
                        amat[ia,ja] = (1.5*(-1.0*sig_A)
                                -2.0*tausIS*(modw_d1*gu2r/abs(zsIS) +w_d1*im*gu2i/zsIS)
                                -d_ee*nuei_p1_p3_t -d_ee*nuei_p1_p3
                                -xnuei*d_ab*xnu_p1_1*ft2)
                        bmat[ia,ja] = 1.5*(-1.0*sig_B)

                        ja = 10*nbasis+jb + ja0
                        amat[ia,ja] = -k_par*vsIS +(am +bm)*gradB*vsIS
                        bmat[ia,ja] = 0.0

                        ja = 11*nbasis+jb + ja0
                        amat[ia,ja] = -bm*3.0*gradB*vsIS
                        bmat[ia,ja] = 0.0
                        #****************
                        # p_par_g trapped terms
                        #****************
                        ja = 12*nbasis+jb + ja0
                        amat[ia,ja] = phi_A + psi_AN
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = 13*nbasis+jb + ja0
                        amat[ia,ja] = -0.5*sig_A
                        bmat[ia,ja] = -0.5*sig_B

                        ja = 14*nbasis+jb + ja0
                        amat[ia,ja] = 1.5*sig_A
                        bmat[ia,ja] = 1.5*sig_B
#*************************************************************
# p_tot_g equ #10
#*************************************************************
                        ia = 9*nbasis+ib + ia0

                        phi_A = N_j*im*w_s*(rlnsIS*gp3+rltsIS*1.5*(gr33-gp3))
                        if(vpar_model_in==0)
                            phi_A = phi_A + N_j*E_i*kpar_gp3p0*vpar[is]
                        end
                        phi_B = -gp3*E_i*N_j
                        phi_A = phi_A +xnu_phi_b*xnuei*xnu_p3_b*phi_B
                        sig_A = zero(K)
                        sig_B = zero(K)
                        psi_A = zero(K)
                        psi_B = zero(K)
                        phi_AU = zero(K)
                        phi_BU = zero(K)
                        psi_AN = zero(K)
                        psi_BN = zero(K)
                        if(use_bpar_in)
                            sig_A = (-betae_sig*(asJS*tausJS*zsIS/massIS)*
                                    (im*w_s*(rlnsIS*g10p3 + rltsIS*1.5*(g10r33-g10p3))))
                            sig_B = (betae_sig*g10p3*asJS*tausJS*zsIS*zsIS
                                    /(tausIS*massIS))
                            sig_A = sig_A - damp_sig*sig_B
                        end
                        if(use_bper_in)
                            psi_A = -betae_psi*J_j*im*w_s*vpar_shear[is]*gr13b0
                            if(vpar_model_in==0)
                                psi_A = psi_A -betae_psi*J_j*im*w_cd*vpar[is]*wdgr13b0
                                psi_B = betae_psi*M_i*J_j*vpar[is]*gr13b0/vsIS
                                phi_AU = (betae_psi*U0*J_j*(im*w_s*(rlnsIS*gp3bp+rltsIS*1.5*(gr33bp-gp3bp))
                                        + E_i*kpar_gp3bp*vpar[is]))
                                phi_BU = -betae_psi*U0*gp3bp*E_i*J_j
                                psi_AN = betae_psi*U0*N_j*im*(w_cd*vpar[is]*wdgr13bp+w_s*vpar_shear[is]*gr13bp)
                                psi_BN = -betae_psi*U0*M_i*N_j*vpar[is]*gr13bp/vsIS
                            end
                        end
                        #****************
                        #  p_tot_g untrapped terms
                        #****************
                        ja = jb + ja0
                        amat[ia,ja] = phi_A + psi_AN +d_ee*nuei_p3_n*bn
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = nbasis+jb + ja0
                        amat[ia,ja] = psi_A + phi_AU
                        bmat[ia,ja] = psi_B + phi_BU

                        ja = 2*nbasis+jb + ja0
                        amat[ia,ja] = (-0.5*sig_A +d_ee*nuei_p3_p1*bp1
                                -d_ee*(1.0 - ft2)*(bn*an*nuei_p3_n + bp3*ap3*nuei_p3_p3 + bp1*ap1*nuei_p3_p1))
                        bmat[ia,ja] = -0.5*sig_B

                        ja = 3*nbasis+jb + ja0
                        amat[ia,ja] = (1.5*sig_A +d_ee*nuei_p3_p3*bp3
                                +d_ee*(1.0 - ft2)*(bn*an*nuei_p3_n + bp3*ap3*nuei_p3_p3 + bp1*ap1*nuei_p3_p1))
                        bmat[ia,ja] = 1.5*sig_B

                        ja = 4*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 5*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
                        #****************
                        #  p_tot_g ghost terms
                        #****************
                        ja = 6*nbasis+jb + ja0
                        amat[ia,ja] = (-1.0*(phi_A + psi_AN)
                                +2.0*tausIS*(modw_d1*gu3rgt1/abs(zsIS) +w_d1*im*gu3igt1/zsIS)
                                +2.0*tausIS*(modw_d1*gu4rgt3/abs(zsIS) +w_d1*im*gu4igt3/zsIS)
                                -d_ee*nuei_p3_n)
                        bmat[ia,ja] = -1.0*(phi_B + psi_BN)

                        ja = 7*nbasis+jb + ja0
                        amat[ia,ja] = k_par1*grad_gu3*vsIS -1.0*(psi_A + phi_AU)
                        bmat[ia,ja] = -1.0*(psi_B + phi_BU)

                        ja = 8*nbasis+jb + ja0
                        amat[ia,ja] = (-0.5*(-1.0*sig_A)
                            -im*w_d1*(tausIS/zsIS)*0.5*wdgu3
                            -2.0*tausIS*(modw_d1*gu3r/abs(zsIS) +w_d1*im*gu3i/zsIS)
                            -d_ee*nuei_p3_p1)
                        bmat[ia,ja] = -0.5*(-1.0*sig_B)

                        ja = 9*nbasis+jb + ja0
                        amat[ia,ja] = (1.5*(-1.0*sig_A)
                                -im*w_d1*(tausIS/zsIS)*1.5*wdgu33
                                -2.0*tausIS*(modw_d1*gu4r/abs(zsIS) +w_d1*im*gu4i/zsIS)
                                -d_ee*nuei_p3_p3
                                +xnuei*d_ab*xnu_p3_b)
                        bmat[ia,ja] = d_ab + 1.5*(-1.0*sig_B)

                        ja = 10*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 11*nbasis+jb + ja0
                        amat[ia,ja] = -k_par*vsIS + am*gradB*vsIS
                        bmat[ia,ja] = 0.0
                        #****************
                        #  p_tot_g trapped terms
                        #****************
                        ja = 12*nbasis+jb + ja0
                        amat[ia,ja] = phi_A + psi_AN
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = 13*nbasis+jb + ja0
                        amat[ia,ja] = -0.5*sig_A
                        bmat[ia,ja] = -0.5*sig_B

                        ja = 14*nbasis+jb + ja0
                        amat[ia,ja] = 1.5*sig_A
                        bmat[ia,ja] = 1.5*sig_B
#*************************************************************
#  q_par_g equ #11
#*************************************************************
                        ia = 10*nbasis+ib + ia0
                    #
                        phi_A = N_j*im*w_s*vpar_shear[is]*gr11/vsIS
                        if(vpar_model_in==0)
                            phi_A = phi_A + N_j*im*w_cd*wdgr11p0*vpar[is]/vsIS
                            phi_B = -E_i*N_j*gr11*vpar[is]/vsIS
                        end
                        sig_A = zero(K)
                        sig_B = zero(K)
                        psi_A = zero(K)
                        psi_B = zero(K)
                        phi_AU = zero(K)
                        phi_BU = zero(K)
                        psi_AN = zero(K)
                        psi_BN = zero(K)
                        if(use_bper_in)
                            psi_A = -betae_psi*J_j*vsIS*im*w_s*(rlnsIS*gr11b0+1.5*rltsIS*(gw113b0-gr11b0))
                            psi_B = gr11b0*betae_psi*M_i*J_j
                            psi_A = psi_A - damp_psi*psi_B
                            if(vpar_model_in==0)
                                psi_A = psi_A  -betae_psi*J_j*M_i*kpar_gr11b0*vpar[is]
                                phi_AU = betae_psi*U0*J_j*im*(w_s*vpar_shear[is]*gr11bp +w_cd*wdgr11bp*vpar[is])/vsIS
                                phi_BU = -betae_psi*U0*E_i*J_j*gr11bp*vpar[is]/vsIS
                                psi_AN = -betae_psi*U0*N_j*(vsIS*im*w_s*(rlnsIS*gr11bp+1.5*rltsIS*(gw113bp-gr11bp))
                                        +M_i*kpar_gr11bp*vpar[is])
                                psi_BN = -gr11bp*betae_psi*U0*M_i*N_j
                            end
                        end
                        #****************
                        # q_par_g untrapped terms
                        #****************
                        ja = jb + ja0
                        amat[ia,ja] = phi_A + psi_AN
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = nbasis+jb + ja0
                        amat[ia,ja] = psi_A + phi_AU +d_ee*ft3*nuei_q1_u
                        bmat[ia,ja] = psi_B + phi_BU

                        ja = 2*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 3*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 4*nbasis+jb + ja0
                        amat[ia,ja] = (d_ee*ft5*nuei_q1_q1
                                -d_ee*(1.0 -ft2)*(ft3*nuei_q1_u*1.25 +ft3*nuei_q1_q3*35.0/12.0+ft5*nuei_q1_q1*6.25))
                        bmat[ia,ja] = 0.0

                        ja = 5*nbasis+jb + ja0
                        amat[ia,ja] = (d_ee*ft3*nuei_q1_q3
                                +d_ee*(1.0 -ft2)*(ft3*nuei_q1_u*2.25 +ft3*nuei_q1_q3*5.25 +ft5*nuei_q1_q1*11.25))
                        bmat[ia,ja] = 0.0
                        #****************
                        # q_par_g ghost terms
                        #****************
                        ja = 6*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(phi_A + psi_AN) + k_par1*kpar_gb1gt1*vsIS
                        bmat[ia,ja] = -1.0*(phi_B + psi_BN)

                        ja = 7*nbasis+jb + ja0
                        amat[ia,ja] = (-1.0*(psi_A + phi_AU) +modk_par1*modkpar_gd1gu1*vsIS
                                - tausIS*(modw_d1*gu5r/abs(zsIS) + w_d1*im*gu5i/zsIS)
                                -d_ee*nuei_q1_u_t -d_ee*nuei_q1_u
                                +xnuei*(d_ab*xnu_q1_u_1*ft2 - d_ij*gu1*xnu_q1_q1_1
                                - d_ij*gu3*xnu_q1_q3_1*ft2))
                        bmat[ia,ja] = -1.0*(psi_B + phi_BU)

                        ja = 8*nbasis+jb + ja0
                        amat[ia,ja] =  (-k_par1*(kpar_gu1 -grad_gu1 - gradgr11p1)*vsIS
                                - k_par1*kpar_gb1*vsIS + (am+bm*1.5)*gradBgu1*vsIS)
                        bmat[ia,ja] = 0.0

                        ja = 9*nbasis+jb + ja0
                        amat[ia,ja] = -bm*4.5*gradBgu3*vsIS
                        bmat[ia,ja] = 0.0

                        ja = 10*nbasis+jb + ja0
                        amat[ia,ja] = (-modk_par1*modkpar_gd1*vsIS
                                - tausIS*(modw_d1*gu6r/abs(zsIS) + w_d1*im*gu6i/zsIS)
                                -d_ee*nuei_q1_q1_t -d_ee*nuei_q1_q1
                                +xnuei*d_ab*(xnu_q1_q1_1 + xnu_q1_b))
                        bmat[ia,ja] = d_ab

                        ja = 11*nbasis+jb + ja0
                        amat[ia,ja] = (- tausIS*(modw_d1*gu7r/abs(zsIS) + w_d1*im*gu7i/zsIS)
                                -d_ee*nuei_q1_q3_t -d_ee*nuei_q1_q3
                                +xnuei*d_ab*xnu_q1_q3_1*ft2)
                        bmat[ia,ja] = 0.0
                        #****************
                        # q_par_g trapped terms
                        #****************
                        ja = 12*nbasis+jb + ja0
                        amat[ia,ja] = phi_A + psi_AN
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = 13*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 14*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
#*************************************************************
# q_tot_g equ #12
#*************************************************************
                        ia = 11*nbasis+ib + ia0

                        phi_A = N_j*im*w_s*vpar_shear[is]*gr13/vsIS
                        if(vpar_model_in==0)
                            phi_A = phi_A + N_j*im*w_cd*wdgr13p0*vpar[is]/vsIS
                            phi_B = -E_i*N_j*gr13*vpar[is]/vsIS
                        end
                        sig_A = zero(K)
                        sig_B = zero(K)
                        psi_A = zero(K)
                        psi_B = zero(K)
                        phi_AU = zero(K)
                        phi_BU = zero(K)
                        psi_AN = zero(K)
                        psi_BN = zero(K)
                        if(use_bper_in)
                            psi_A = -betae_psi*J_j*vsIS*im*w_s*(rlnsIS*gr13b0+1.5*rltsIS*(gw133b0-gr13b0))
                            psi_B = gr13b0*betae_psi*M_i*J_j
                            psi_A = psi_A - damp_psi*psi_B
                            if(vpar_model_in==0)
                                psi_A = psi_A  -betae_psi*J_j*M_i*kpar_gr13b0*vpar[is]
                                phi_AU = betae_psi*U0*J_j*im*(w_s*vpar_shear[is]*gr13bp +w_cd*wdgr13bp*vpar[is])/vsIS
                                phi_BU = -betae_psi*U0*E_i*J_j*gr13bp*vpar[is]/vsIS
                                psi_AN = betae_psi*U0*N_j*(vsIS*im*w_s*(rlnsIS*gr13bp+1.5*rltsIS*(gw133bp-gr13bp))
                                        +M_i*kpar_gr13bp*vpar[is])
                                psi_BN = -gr13bp*betae_psi*U0*M_i*N_j
                            end
                        end
                        #****************
                        # q_tot_g untrapped terms
                        #****************
                        ja = jb + ja0
                        amat[ia,ja] = phi_A + psi_AN
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = nbasis+jb + ja0
                        amat[ia,ja] = psi_A + phi_AU + d_ee*ft3*nuei_q3_u
                        bmat[ia,ja] = psi_B + phi_BU

                        ja = 2*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 3*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 4*nbasis+jb + ja0
                        amat[ia,ja] = (d_ee*ft5*nuei_q3_q1
                                -d_ee*(1.0 -ft2)*(ft3*nuei_q3_u*1.25 +ft3*nuei_q3_q3*35.0/12.0 +ft5*nuei_q3_q1*6.25))
                        bmat[ia,ja] = 0.0

                        ja = 5*nbasis+jb + ja0
                        amat[ia,ja] = (d_ee*ft3*nuei_q3_q3
                                +d_ee*(1.0 -ft2)*(ft3*nuei_q3_u*2.25 + ft3*nuei_q3_q3*5.25 + ft5*nuei_q3_q1*11.25))
                        bmat[ia,ja] = 0.0
                        #****************
                        # q_tot_g ghost terms
                        #****************
                        ja = 6*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(phi_A + psi_AN) + k_par1*(kpar_gb3gt3 -dgr13 + kpar_gb33gt1)*vsIS
                        bmat[ia,ja] = -1.0*(phi_B + psi_BN)

                        ja = 7*nbasis+jb + ja0
                        amat[ia,ja] = (-1.0*(psi_A + phi_AU) +
                                modk_par1*(modkpar_gd3gu3 + modkpar_gd33gu1)*vsIS
                                - tausIS*(modw_d1*gu8r/abs(zsIS) + w_d1*im*gu8i/zsIS)
                                -d_ee*nuei_q3_u_t -d_ee*nuei_q3_u
                                +xnuei*(d_ab*xnu_q3_u_1 - d_ij*gu3*xnu_q3_q3_1))
                        bmat[ia,ja] = -1.0*(psi_B + phi_BU)

                        ja = 8*nbasis+jb + ja0
                        amat[ia,ja] = (- k_par1*(kpar_gu3 -grad_gu3 - gradgr13p1)*vsIS
                                - k_par1*kpar_gb33*vsIS + (am+bm*0.5)*gradBgu3*vsIS)
                        bmat[ia,ja] = 0.0

                        ja = 9*nbasis+jb + ja0
                        amat[ia,ja] = - k_par1*kpar_gb3*vsIS - bm*1.5*gradBgu33*vsIS
                        bmat[ia,ja] = 0.0

                        ja = 10*nbasis+jb + ja0
                        amat[ia,ja] = (-modk_par1*modkpar_gd33*vsIS
                                - tausIS*(modw_d1*gu9r/abs(zsIS) + w_d1*im*gu9i/zsIS)
                                -d_ee*nuei_q3_q1)
                        bmat[ia,ja] = 0.0

                        ja = 11*nbasis+jb + ja0
                        amat[ia,ja] = (-modk_par1*modkpar_gd3*vsIS
                                - tausIS*(modw_d1*gu10r/abs(zsIS) +w_d1* im*gu10i/zsIS)
                                -d_ee*nuei_q3_q3_t -d_ee*nuei_q3_q3
                                +xnuei*d_ab*(xnu_q3_q3_1 + xnu_q3_b))
                        bmat[ia,ja] = d_ab
                        #****************
                        # q_tot_g trapped terms
                        #****************
                        ja = 12*nbasis+jb + ja0
                        amat[ia,ja] = phi_A + psi_AN
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = 13*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 14*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
#*************************************************************
# n_t equ #13
#*************************************************************
                        ia = 12*nbasis+ib + ia0

                        phi_A = N_j*im*w_s*(rlnsIS*gn + rltsIS*1.5*(gp3-gn))
                        phi_B = -gn*E_i*N_j
                        phi_A = phi_A +xnu_phi_b*xnuei*xnu_n_b*phi_B
                        sig_A = zero(K)
                        sig_B = zero(K)
                        psi_A = zero(K)
                        psi_B = zero(K)
                        phi_AU = zero(K)
                        phi_BU = zero(K)
                        psi_AN = zero(K)
                        psi_BN = zero(K)
                        if(use_bpar_in)
                            sig_A = (-betae_sig*(asJS*tausJS*zsIS/massIS)*
                                    (im*w_s*(rlnsIS*g10n + rltsIS*1.5*(g10p3-g10n))))
                            sig_B = (betae_sig*g10n*asJS*tausJS*zsIS*zsIS
                                    /(tausIS*massIS))
                            sig_A = sig_A - damp_sig*sig_B
                        end
                        if(use_bper_in)
                            psi_A = -betae_psi*J_j*im*w_s*vpar_shear[is]*gp1b0
                            if(vpar_model_in==0)
                                psi_A = psi_A -betae_psi*J_j*im*w_cd*vpar[is]*wdgp1b0
                                psi_B = betae_psi*M_i*J_j*vpar[is]*gp1b0/vsIS
                                phi_AU = betae_psi*U0*J_j*im*w_s*(rlnsIS*gnbp + rltsIS*1.5*(gp3bp-gnbp))
                                phi_BU = -betae_psi*U0*gnbp*E_i*J_j
                                psi_AN = betae_psi*U0*N_j*im*(w_cd*vpar[is]*wdgp1bp+w_s*vpar_shear[is]*gp1bp)
                                psi_BN = -betae_psi*U0*M_i*N_j*vpar[is]*gp1bp/vsIS
                            end
                        end
                        #****************
                        # n_t untrapped terms
                        #****************
                        ja = jb + ja0
                        amat[ia,ja] = phi_A + psi_AN +d_ee*nuei_n_n*bn
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 2*nbasis+jb + ja0
                        amat[ia,ja] = (-0.5*sig_A +d_ee*nuei_n_p1*bp1
                                -d_ee*(1.0 - ft2)*(bn*an*nuei_n_n + bp3*ap3*nuei_n_p3 + bp1*ap1*nuei_n_p1))
                        bmat[ia,ja] = -0.5*sig_B

                        ja = 3*nbasis+jb + ja0
                        amat[ia,ja] = (1.5*sig_A +d_ee*nuei_n_p3*bp3
                                +d_ee*(1.0 - ft2)*(bn*an*nuei_n_n + bp3*ap3*nuei_n_p3 + bp1*ap1*nuei_n_p1))
                        bmat[ia,ja] = 1.5*sig_B

                        ja = 4*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 5*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
                    #****************
                    # n_t ghost terms
                    #****************
                        ja = 6*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(phi_A + psi_AN)
                        bmat[ia,ja] = -1.0*(phi_B + psi_BN)

                        ja = 7*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 8*nbasis+jb + ja0
                        amat[ia,ja] = -0.5*(-1.0*sig_A)
                        bmat[ia,ja] = -0.5*(-1.0*sig_B)

                        ja = 9*nbasis+jb + ja0
                        amat[ia,ja] = 1.5*(-1.0*sig_A)
                        bmat[ia,ja] = 1.5*(-1.0*sig_B)

                        ja = 10*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 11*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
                    #****************
                    # n_t trapped terms
                    #****************
                        ja = 12*nbasis+jb + ja0
                        amat[ia,ja] = phi_A + psi_AN -d_ee*nuei_n_n +xnuei*d_ab*xnu_n_b
                        bmat[ia,ja] = d_ab + phi_B + psi_BN

                        ja = 13*nbasis+jb + ja0
                        amat[ia,ja] = (-0.5*sig_A -0.5*im*w_dg*tausIS/zsIS
                                -d_ee*nuei_n_p1)
                        bmat[ia,ja] = -0.5*sig_B

                        ja = 14*nbasis+jb + ja0
                        amat[ia,ja] = (1.5*sig_A -1.5*im*w_dg*tausIS/zsIS
                                -d_ee*nuei_n_p3)
                        bmat[ia,ja] = 1.5*sig_B
#*************************************************************
# p_par_t equ #14
#*************************************************************
                        ia = 13*nbasis+ib + ia0

                        phi_A = N_j*im*w_s*(rlnsIS*gp1 + rltsIS*1.5*(gr13-gp1))
                        phi_B = -gp1*E_i*N_j
                        phi_A = phi_A +xnu_phi_b*xnuei*xnu_p1_b*phi_B
                        sig_A = zero(K)
                        sig_B = zero(K)
                        psi_A = zero(K)
                        psi_B = zero(K)
                        phi_AU = zero(K)
                        phi_BU = zero(K)
                        psi_AN = zero(K)
                        psi_BN = zero(K)
                        if(use_bpar_in)
                            sig_A = (-betae_sig*(asJS*tausJS*zsIS/massIS)*
                                    (im*w_s*(rlnsIS*g10p1 + rltsIS*1.5*(g10r13-g10p1))))
                            sig_B = (betae_sig*g10p1*asJS*tausJS*zsIS*zsIS
                                    /(tausIS*massIS))
                            sig_A = sig_A - damp_sig*sig_B
                        end
                        if(use_bper_in)
                            psi_A = -betae_psi*J_j*im*w_s*vpar_shear[is]*gr11b0
                            if(vpar_model_in==0)
                                psi_A = psi_A  -betae_psi*J_j*im*w_cd*vpar[is]*wdgr11b0
                                psi_B = betae_psi*M_i*J_j*vpar[is]*gr11b0/vsIS
                                phi_AU = betae_psi*U0*J_j*im*w_s*(rlnsIS*gp1bp + rltsIS*1.5*(gr13bp-gp1bp))
                                phi_BU = -betae_psi*U0*gp1bp*E_i*J_j
                                psi_AN = betae_psi*U0*N_j*im*(w_cd*vpar[is]*wdgr11bp+ w_s*vpar_shear[is]*gr11bp)
                                psi_BN = -betae_psi*U0*M_i*N_j*vpar[is]*gr11bp/vsIS
                            end
                        end
                    #****************
                    # p_par_t  untrapped terms
                    #****************
                        ja = jb + ja0
                        amat[ia,ja] = phi_A + psi_AN +d_ee*nuei_p1_n*bn
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 2*nbasis+jb + ja0
                        amat[ia,ja] = (-0.5*sig_A +d_ee*nuei_p1_p1*bp1
                                -d_ee*(1.0 - ft2)*(bn*an*nuei_p1_n + bp3*ap3*nuei_p1_p3 + bp1*ap1*nuei_p1_p1))
                        bmat[ia,ja] = -0.5*sig_B

                        ja = 3*nbasis+jb + ja0
                        amat[ia,ja] = 1.5*sig_A
                        bmat[ia,ja] = (1.5*sig_B +d_ee*nuei_p1_p3*bp3
                                +d_ee*(1.0 - ft2)*(bn*an*nuei_p1_n + bp3*ap3*nuei_p1_p3 + bp1*ap1*nuei_p1_p1))

                        ja = 4*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 5*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
                    #****************
                    # p_par_t  ghost terms
                    #****************
                        ja = 6*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(phi_A + psi_AN)
                        bmat[ia,ja] = -1.0*(phi_B + psi_BN)

                        ja = 7*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 8*nbasis+jb + ja0
                        amat[ia,ja] = -0.5*(-1.0*sig_A)
                        bmat[ia,ja] = -0.5*(-1.0*sig_B)

                        ja = 9*nbasis+jb + ja0
                        amat[ia,ja] = 1.5*(-1.0*sig_A)
                        bmat[ia,ja] = 1.5*(-1.0*sig_B)

                        ja = 10*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 11*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
                    #****************
                    # p_par_t  trapped terms
                    #****************

                        ja = 12*nbasis+jb + ja0
                        amat[ia,ja] = (phi_A +psi_AN
                                +2.0*tausIS*(modw_d1*gu1rgt1/abs(zsIS) +w_d1*im*gu1igt1/zsIS)
                                +2.0*tausIS*(modw_d1*gu2rgt3/abs(zsIS) +w_d1*im*gu2igt3/zsIS)
                                -d_ee*nuei_p1_n -xnuei*d_ij*xnu_p1_1*(gt1-ft2*gt3))
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = 13*nbasis+jb + ja0
                        amat[ia,ja] =  (-0.5*sig_A
                                -im*w_d1*(tausIS/zsIS)*(0.5*wdgu1+1.5*wdgu3)
                                -2.0*tausIS*(modw_d1*gu1r/abs(zsIS) +w_d1*im*gu1i/zsIS)
                                -d_ee*nuei_p1_p1_t -d_ee*nuei_p1_p1
                                +xnuei*d_ab*(xnu_p1_1 + xnu_p1_b))
                        bmat[ia,ja] = d_ab -0.5*sig_B

                        ja = 14*nbasis+jb + ja0
                        amat[ia,ja] = (1.5*sig_A
                                -2.0*tausIS*(modw_d1*gu2r/abs(zsIS) +w_d1*im*gu2i/zsIS)
                                -d_ee*nuei_p1_p3_t -d_ee*nuei_p1_p3
                                -xnuei*d_ab*xnu_p1_1*ft2)
                        bmat[ia,ja] = 1.5*sig_B
#*************************************************************
# p_tot_t equ #15
#*************************************************************
                        ia = 14*nbasis+ib + ia0

                        phi_A = N_j*im*w_s*(rlnsIS*gp3+rltsIS*1.5*(gr33-gp3))
                        phi_B = -gp3*E_i*N_j
                        phi_A = phi_A +xnu_phi_b*xnuei*xnu_p3_b*phi_B
                        sig_A = zero(K)
                        sig_B = zero(K)
                        psi_A = zero(K)
                        psi_B = zero(K)
                        phi_AU = zero(K)
                        phi_BU = zero(K)
                        psi_AN = zero(K)
                        psi_BN = zero(K)
                        if(use_bpar_in)
                            sig_A =( -betae_sig*(asJS*tausJS*zsIS/massIS)*
                                    (im*w_s*(rlnsIS*g10p3 + rltsIS*1.5*(g10r33-g10p3))))
                            sig_B = (betae_sig*g10p3*asJS*tausJS*zsIS*zsIS
                                    /(tausIS*massIS))
                            sig_A = sig_A - damp_sig*sig_B
                        end
                        if(use_bper_in)
                            psi_A = -betae_psi*J_j*im*w_s*vpar_shear[is]*gr13b0
                            if(vpar_model_in==0)
                                psi_A = psi_A -betae_psi*J_j*im*w_cd*vpar[is]*wdgr13b0
                                psi_B = betae_psi*M_i*J_j*vpar[is]*gr13b0/vsIS
                                phi_AU = betae_psi*U0*J_j*im*w_s*(rlnsIS*gp3bp+rltsIS*1.5*(gr33bp-gp3bp))
                                phi_BU = -betae_psi*U0*gp3bp*E_i*J_j
                                psi_AN = betae_psi*U0*N_j*im*(w_cd*vpar[is]*wdgr13bp+w_s*vpar_shear[is]*gr13bp)
                                psi_BN = -betae_psi*U0*M_i*N_j*vpar[is]*gr13bp/vsIS
                            end
                        end
                    #****************
                    # p_tot_t untrapped terms
                    #****************
                        ja = jb + ja0
                        amat[ia,ja] = phi_A + psi_AN +d_ee*nuei_p3_n*bn
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 2*nbasis+jb + ja0
                        amat[ia,ja] = (-0.5*sig_A +d_ee*nuei_p3_p1*bp1
                                -d_ee*(1.0 - ft2)*(bn*an*nuei_p3_n + bp3*ap3*nuei_p3_p3 + bp1*ap1*nuei_p3_p1))
                        bmat[ia,ja] = -0.5*sig_B

                        ja = 3*nbasis+jb + ja0
                        amat[ia,ja] = (1.5*sig_A +d_ee*nuei_p3_p3*bp3
                                +d_ee*(1.0 - ft2)*(bn*an*nuei_p3_n + bp3*ap3*nuei_p3_p3 + bp1*ap1*nuei_p3_p1))
                        bmat[ia,ja] = 1.5*sig_B

                        ja = 4*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 5*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
                    #****************
                    # p_tot_t ghost terms
                    #****************
                        ja = 6*nbasis+jb + ja0
                        amat[ia,ja] = -1.0*(phi_A + psi_AN)
                        bmat[ia,ja] = -1.0*(phi_B + psi_BN)

                        ja = 7*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 8*nbasis+jb + ja0
                        amat[ia,ja] = -0.5*(-1.0*sig_A)
                        bmat[ia,ja] = -0.5*(-1.0*sig_B)

                        ja = 9*nbasis+jb + ja0
                        amat[ia,ja] = 1.5*(-1.0*sig_A)
                        bmat[ia,ja] = 1.5*(-1.0*sig_B)

                        ja = 10*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0

                        ja = 11*nbasis+jb + ja0
                        amat[ia,ja] = 0.0
                        bmat[ia,ja] = 0.0
                    #****************
                    # p_tot_t trapped terms
                    #****************
                        ja = 12*nbasis+jb + ja0
                        amat[ia,ja] = (phi_A + psi_AN
                                +2.0*tausIS*(modw_d1*gu3rgt1/abs(zsIS) +w_d1*im*gu3igt1/zsIS)
                                +2.0*tausIS*(modw_d1*gu4rgt3/abs(zsIS) +w_d1*im*gu4igt3/zsIS)
                                -d_ee*nuei_p3_n )
                        bmat[ia,ja] = phi_B + psi_BN

                        ja = 13*nbasis+jb + ja0
                        amat[ia,ja] = (-0.5*sig_A -im*w_d1*(tausIS/zsIS)*0.5*wdgu3
                                -2.0*tausIS*(modw_d1*gu3r/abs(zsIS) +w_d1*im*gu3i/zsIS)
                                -d_ee*nuei_p3_p1)
                        bmat[ia,ja] = -0.5*sig_B

                        ja = 14*nbasis+jb + ja0
                        amat[ia,ja] = (1.5*sig_A -im*w_d1*(tausIS/zsIS)*1.5*wdgu33
                                -2.0*tausIS*(modw_d1*gu4r/abs(zsIS) +w_d1*im*gu4i/zsIS)
                                -d_ee*nuei_p3_p3
                                +xnuei*d_ab*xnu_p3_b)
                        bmat[ia,ja] = d_ab + 1.5*sig_B

                    end # nroot > 6
                end # end of js loop
            end # end of ib loop
        end # end of jb loop
    end # end of is loop
    #*************************************************************
    # find the eigenvalues and eigenvectors
    #*************************************************************

    if inputs.SMALL != 0.0 && !inputs.FIND_EIGEN
        @warn "FIND_EIGEN is being ignored since inputs.SMALL!=0.0"
    end

    # calculate eigenvalues/eigenvectors
    if inputs.SMALL == 0.0 && !inputs.FIND_EIGEN && !inputs.IFLUX
        if inputs.FIND_WIDTH
            error("If FIND_EIGEN false, FIND_WIDTH should also be false")
        end

        sigma = 0.0
        if !isnan(inputs.EIGEN_SPECTRUM[ky_index])
            sigma = inputs.EIGEN_SPECTRUM[ky_index]  #array of eigenvalues from input file
        end

        if sigma != 0.0
            try              
                nev1=inputs.NMODES            
                L=construct_linear_map(sparse(amat), sparse(bmat), sigma)
                λ, v, _ = KrylovKit.eigsolve(L, size(amat)[1], nev1, :LM) 
               # printl("Success!!!!----------------------------")  - this solver almost never works
                return λ, v[1], NaN, NaN
            catch e
                @warn "KrylovKit.eigsolve() can't find eigen for ky = $(inputs.KY_SPECTRUM[ky_index]), using gesv!+geev! to find all eigenvalues: $(e)"
            end
        else
            @warn "no growth rate initial guess given for ky = $(inputs.KY_SPECTRUM[ky_index]), using gesv!+geev! to find all eigenvalues"
        end
    end

    use_tm = inputs.USE_TRANSPORT_MODEL
    if !use_gpu
        if !use_tm
            # Eigenvalues via the standard solve B⁻¹A (gesv! + geev!) instead of
            # the generalized ggev!. For these well-conditioned pencils the two
            # agree to ~1e-11 relative, but gesv!+geev! is ~4x faster on the small
            # dense matrices here (e.g. 270×270: 179 ms -> 45 ms). This also makes
            # the CPU non-transport path match the GPU path, which already solves
            # this way (see the use_gpu branch below).
            if inputs.IFLUX || find_eigenvector
                if workspace === nothing
                    amat_copy = copy(amat)
                    bmat_copy = copy(bmat)
                else
                    # same values, caller-owned storage (the solve destroys both copies)
                    amat_copy = copyto!(workspace.amat_copy, amat)
                    bmat_copy = copyto!(workspace.bmat_copy, bmat)
                end
                alpha = _standard_eigenvalues_via_solve(amat_copy, bmat_copy; use_gpu=use_gpu)
            else
                alpha = _standard_eigenvalues_via_solve(amat, bmat; use_gpu=use_gpu)
            end
            return alpha, fill(NaN*im,(1,1))
        end

        if inputs.IFLUX || find_eigenvector
            amat_copy = copy(amat)
            bmat_copy = copy(bmat)
            alpha = _standard_eigenvalues_via_solve(amat_copy, bmat_copy; use_gpu=use_gpu)
        else
            alpha = _standard_eigenvalues_via_solve(amat, bmat; use_gpu=use_gpu)
        end

    else
        # use_gpu=true: skip the CPU _generalized_eigenvalues call that the
        # non-transport-model path would otherwise run (its result is unused).
        amat_copy = copy(amat)
        bmat_copy = copy(bmat)
        alpha = _standard_eigenvalues_via_solve(amat_copy, bmat_copy; use_gpu=use_gpu)
    end

    return alpha, fill(NaN*im,(1,1))

end

# Fast path: use LAPACK for ComplexF64, generic Julia otherwise
function _generalized_eigenvalues(A::Matrix{ComplexF64}, B::Matrix{ComplexF64})
    (alpha, beta, _, _) = ggev!('N','N', A, B)
    return alpha ./ beta
end

function _generalized_eigenvalues(A::Matrix{K}, B::Matrix{K}) where {K<:Complex}
    F = eigen(A, B)
    return F.values
end

# Debug/benchmark hook: when TJLF_DUMP_PENCILS=<dir> is set, save up to
# TJLF_DUMP_PENCILS_MAX (default 64) raw (A,B) pencils that reach the dense eigensolve, as
# Serialization blobs pencil_<k>.jls containing (A=Matrix{ComplexF64}, B=Matrix{ComplexF64}).
# Used by TJLFEP/build/ad/benchmark_krylov_eigs.jl to test alternative eigensolvers on REAL
# TJLF matrices instead of random ones. Off (zero overhead beyond one ENV lookup) by default.
const _PENCIL_DUMP_COUNT = Threads.Atomic{Int}(0)
function _maybe_dump_pencil(A::Matrix{ComplexF64}, B::Matrix{ComplexF64})
    dir = get(ENV, "TJLF_DUMP_PENCILS", "")
    isempty(dir) && return
    # TJLF_DUMP_PENCILS_MINSIZE filters to full-nbasis pencils: a grid scan also solves many
    # cheap coarse-nbasis width-search passes (small n), which we usually don't want to profile.
    minsize = parse(Int, get(ENV, "TJLF_DUMP_PENCILS_MINSIZE", "0"))
    size(A, 1) < minsize && return
    maxn = parse(Int, get(ENV, "TJLF_DUMP_PENCILS_MAX", "64"))
    k = Threads.atomic_add!(_PENCIL_DUMP_COUNT, 1) + 1
    k > maxn && return
    mkpath(dir)
    path = joinpath(dir, "pencil_$(lpad(k, 4, '0')).jls")
    Serialization.serialize(path, (A = copy(A), B = copy(B)))
    return
end

# Opt-in eigensolve hook: when set (Ref holds a callable), _standard_eigenvalues_via_solve routes
# the (A,B) pencil to it instead of the dense per-pencil solve, returning Vector{ComplexF64}
# eigenvalues. This is the seam TJLFEP's batched shift-invert driver uses to collect pencils
# across grid combos and solve them in one batched GPU sweep (see TJLFEP kwscale inner=:batched_si).
# Only valid for the eigenvalue-only path (IFLUX=false, no eigenvectors); nothing by default so the
# standard behavior is bit-for-bit unchanged.
const _EIGSOLVE_HOOK = Ref{Any}(nothing)

function _standard_eigenvalues_via_solve(A::Matrix{ComplexF64}, B::Matrix{ComplexF64}; use_gpu::Bool=false)
    _maybe_dump_pencil(A, B)
    let h = _EIGSOLVE_HOOK[]
        h === nothing || return h(A, B)::Vector{ComplexF64}
    end
    if use_gpu
        # _gpu_solve! does getrf(B) + getrs(B,A) + Xgeev!(B⁻¹A) entirely on GPU
        # and returns Vector{ComplexF64} eigenvalues directly.
        return _gpu_solve!(A, B)
    else
        (A, B, _) = gesv!(B, A)
        return geev!('N','N', A)[1]
    end
end

function _standard_eigenvalues_via_solve(A::Matrix{K}, B::Matrix{K}; use_gpu::Bool=false) where {K<:Complex}
    alpha = _herm_eigen(B \ A).values
    return alpha
end

# ── TJLF-owned eigen entrypoints ─────────────────────────────────────────────
# These wrap LinearAlgebra.eigen so that AD (ForwardDiff.Dual element types) can
# be specialized in TJLFForwardDiffExt WITHOUT committing type piracy on
# LinearAlgebra.eigen. Default implementations just forward to Base/LAPACK; the
# extension adds Dual-typed methods to these owned functions instead.
#   _sym_eigen   : real symmetric matrices (modwd!)
#   _herm_eigen  : complex matrices (modkpar!, B\A standard solve)
# Both return an object supporting positional destructuring `(values, vectors)`
# and `.values` / `.vectors` field access (Eigen or NamedTuple).
_sym_eigen(A::Symmetric) = eigen(A)
_herm_eigen(A::AbstractMatrix) = eigen(A)
