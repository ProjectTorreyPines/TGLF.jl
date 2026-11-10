using Test
using TJLF

include("runtests_input.jl")

include("runtests_regressions.jl")

include("runtests_sat.jl")

include("runtests_core.jl")

include("runtests_EM.jl")

include("runtests_tm.jl")

include("runtests_kygrid.jl")

include("runtests_eigen.jl")

include("runtests_batched_si.jl")

include("runtests_workspace.jl")

println("\n*** Starting ForwardDiff AD tests (this may take ~20 min due to Dual-number compilation) ***\n")
include("runtests_ad.jl")
