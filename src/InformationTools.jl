module InformationTools
using Distributions
export omega2, omega2!, fval, fval!, qrfactchk, calculate_f, calculate_omega2, computessr,
       computessr!, poisson_deviance_grid, poisson_deviance_grid_pval, poisson_deviance_pval,
       makeperms

function check_ss_effect(x)
    if x < -3e-9
        error("ss_effect = $x < -3e-9")
    end
end

function omega2!(out, fullrss::AbstractArray, fulldfe, reducedrss::AbstractArray, reduceddfe, n)
    broadcast!(omega2, out, fullrss, fulldfe, reducedrss, reduceddfe, n)
end

function omega2(fullrss, fulldfe, reducedrss, reduceddfe, n)
    ss_effect = reducedrss - fullrss
    check_ss_effect(ss_effect)
    df_effect = reduceddfe - fulldfe
    ms_effect = ss_effect/df_effect
    ms_error = fullrss/fulldfe
    (df_effect*(ms_effect - ms_error))/(df_effect*ms_effect + (n - df_effect)*ms_error)
end

function fval!(out, fullrss::AbstractArray, fulldfe, reducedrss::AbstractArray, reduceddfe)
    broadcast!(fval, out, fullrss, fulldfe, reducedrss, reduceddfe)
end

function fval(fullrss, fulldfe, reducedrss, reduceddfe)
    ss_effect = reducedrss - fullrss
    check_ss_effect(ss_effect)
    (ss_effect/(reduceddfe - fulldfe))./(fullrss/fulldfe)
end

qrfactchk(X) = ((r = rank(X)) == size(X, 2) || error("rank is $(r) but size is $(size(X, 2))"); qrfact(X))

function calculate_f(rast, dmat::LinAlg.QRCompactWY, reduced_dmat::LinAlg.QRCompactWY; stat::Symbol=:omega2)
    reducedrss = computessr(reduced_dmat, rast)
    fullrss = computessr(dmat, rast)
    fulldfe = size(rast, 1) - size(dmat, 2)
    reduceddfe = size(rast, 1) - size(reduced_dmat, 2)
    return fval.(fullrss, fulldfe, reducedrss, reduceddfe)
end

function calculate_omega2(rast, dmat::LinAlg.QRCompactWY, reduced_dmat::LinAlg.QRCompactWY; stat::Symbol=:omega2)
    reducedrss = computessr(reduced_dmat, rast)
    fullrss = computessr(dmat, rast)
    fulldfe = size(rast, 1) - size(dmat, 2)
    reduceddfe = size(rast, 1) - size(reduced_dmat, 2)
    return omega2.(fullrss, fulldfe, reducedrss, reduceddfe, size(rast, 1))
end

computessr(qrf::LinAlg.QRCompactWY, y::StridedMatrix) = vec(sumabs2(view(qrf[:Q]'y, size(qrf, 2)+1:size(qrf, 1), :), 1))
computessr!(qrf::LinAlg.QRCompactWY, y::StridedMatrix) =
    vec(sumabs2(view(LAPACK.gemqrt!('L', 'T', qrf.factors, qrf.T, y), size(qrf, 2)+1:size(qrf, 1), :), 1))
function computessr!(out::StridedVector, qrf::LinAlg.QRCompactWY, y::StridedMatrix)
    sumabs2!(reshape(out, 1, length(out)), view(LAPACK.gemqrt!('L', 'T', qrf.factors, qrf.T, y), size(qrf, 2)+1:size(qrf, 1), :))
    out
end

function cumsum2!(X)
    for j = 2:size(X, 2), i = 1:size(X, 1)
        @inbounds X[i, j] = X[i, j-1] + X[i, j]
    end
    X
end

@inline xlogy{T<:Real}(x::T, y::T) = x > zero(T) ? x * log(y) : zero(x)
function poisson_deviance_grid(Y, ningroup, minlength=1)
    @assert sum(ningroup) == size(Y, 1)
    cs = sum(Y, 1)
    cumsum2!(cs)
    groupcs = fill(-1, length(ningroup), size(Y, 2))
    lls = fill(NaN, size(Y, 2), size(Y, 2))
    for itp = 1:size(Y, 2)
        iidx = 0
        for igrp = 1:length(ningroup)
            grpend = iidx+ningroup[igrp]
            s = 0
            @simd for i = iidx+1:grpend
                @inbounds s += Y[i, itp]
            end
            groupcs[igrp, itp] = s
            iidx = grpend
        end
    end
    cumsum2!(groupcs)

    for istart = 1:size(Y, 2), iend = istart+minlength:size(Y, 2)
        groupll = 0.0
        iidx = 0
        for igrp = 1:length(ningroup)
            k = groupcs[igrp, iend] - (istart == 1 ? 0 : groupcs[igrp, istart-1])
            groupll += -k + xlogy(float(k), k/(ningroup[igrp]*(iend-istart+1)))
        end
        k = cs[iend] - (istart == 1 ? 0 : cs[istart-1])
        ll = -k + xlogy(float(k), k/(size(Y, 1)*(iend-istart+1)))
        lls[iend, istart] = groupll - ll
    end

    return lls
end

function poisson_deviance_grid_pval(Y, ningroup, perm, minlength=1)
    @assert sum(ningroup) == size(Y, 1)
    cs = sum(Y, 1)
    cumsum2!(cs)
    groupcs = fill(-1, length(ningroup), size(Y, 2))

    truevals = poisson_deviance_grid(Y, ningroup, minlength)
    trueval = maximum(truevals)
    isnan(trueval) && return NaN

    ngt = 0
    for ibs = 1:size(perm, 2)
        for itp = 1:size(Y, 2)
            iidx = 0
            for igrp = 1:length(ningroup)
                grpend = iidx+ningroup[igrp]
                s = 0
                @simd for i = iidx+1:grpend
                    @inbounds s += Y[perm[i, ibs], itp]
                end
                groupcs[igrp, itp] = s
                iidx = grpend
            end
        end
        cumsum2!(groupcs)

        bestll = -Inf
        for istart = 1:size(Y, 2), iend = istart+minlength:size(Y, 2)
            groupll = 0.0
            for igrp = 1:length(ningroup)
                k = groupcs[igrp, iend] - (istart == 1 ? 0 : groupcs[igrp, istart-1])
                groupll += -k + xlogy(float(k), k/(ningroup[igrp]*(iend-istart+1)))
            end
            k = cs[iend] - (istart == 1 ? 0 : cs[istart-1])
            ll = -k + xlogy(float(k), k/(size(Y, 1)*(iend-istart+1)))
            newll = groupll - ll
            bestll = max(newll, bestll)
        end

        ngt += bestll >= trueval
    end

    return ngt/size(perm, 2)
end

function poissonll(y, ningroup, perm, ibs)
    iidx = 0
    tk = 0
    groupll = 0.0
    for igrp = 1:length(ningroup)
        grpend = iidx+ningroup[igrp]
        k = 0
        @simd for i = iidx+1:grpend
            @inbounds k += y[perm[i, ibs]]
        end
        groupll += -k + xlogy(float(k), k/ningroup[igrp])
        tk += k
        iidx = grpend
    end
    groupll - (-tk + xlogy(float(tk), tk/length(y)))
end

function poisson_deviance_pval(y, ningroup, perm)
    @assert sum(ningroup) == size(y, 1)
    trueval = poissonll(y, ningroup, 1:length(y), 1)
    isnan(trueval) && return NaN

    ngt = 0
    for ibs = 1:size(perm, 2)
        ll = poissonll(y, ningroup, perm, ibs)
        ngt += ll >= trueval
    end

    return ngt/size(perm, 2)
end

function makeperms(n, nperm)
    perms = zeros(Int, n, nperm)
    inds = collect(1:n)
    for iperm = 1:nperm
        perms[:, iperm] = shuffle!(inds)
    end
    perms
end
end
