using Test
using TJLF

# Single-ky linear solves (USE_TRANSPORT_MODEL=false) through a reusable LSWorkspace must be
# bit-identical to solves that allocate their own buffers, whatever the buffers held before.
@testset "LSWorkspace reuse is bit-identical" begin
    function single_ky_input(ky, width, nbasis)
        input = readInput(joinpath(@__DIR__, "test_EM", "input.tglf"))
        input.USE_TRANSPORT_MODEL = false
        input.KYGRID_MODEL = 0
        # KYGRID_MODEL=0 uses NKY spectrum entries; resize the per-ky arrays to match
        input.KY_SPECTRUM = fill(ky, input.NKY)
        input.WIDTH_SPECTRUM = fill(width, input.NKY)
        input.EIGEN_SPECTRUM = fill(complex(0.0, 0.0), input.NKY)
        input.KY = ky
        input.NMODES = 4
        input.NBASIS_MIN = nbasis
        input.NBASIS_MAX = nbasis
        input.WIDTH = width
        input.FIND_WIDTH = false
        input.IFLUX = true
        input.VEXB_SHEAR = 0.0   # single-ky path has no reference spectrum for the ExB pass
        return input
    end
    same(a, b) = isequal(a.eigenvalue, b.eigenvalue) && isequal(a.QL_weights, b.QL_weights) &&
                 isequal(a.field_weight_out, b.field_weight_out)

    ws = TJLF.LSWorkspace()
    n_unstable = 0
    # nbasis changes mid-sequence to exercise the resize path
    for (ky, width, nbasis) in ((0.3, 1.65, 4), (0.5, 1.2, 4), (0.3, 1.65, 6), (1.0, 0.8, 4))
        ref = TJLF.run(single_ky_input(ky, width, nbasis))
        got = TJLF.run(single_ky_input(ky, width, nbasis); workspace=ws)
        @test same(ref, got)
        # poison every reusable buffer: the next solve must not read stale contents
        for m in (ws.amat, ws.bmat, ws.amat_copy, ws.bmat_copy, ws.zmat)
            fill!(m, complex(NaN, NaN))
        end
        fill!(ws.eigenvector, complex(NaN, NaN))
        got2 = TJLF.run(single_ky_input(ky, width, nbasis); workspace=ws)
        @test same(ref, got2)
        n_unstable += count(>(0), ref.eigenvalue[:, 1, 1])
    end
    @test n_unstable > 0   # the comparison covered real eigenvector / QL-weight work
end
