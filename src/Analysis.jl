module Analysis
using Sessions, TaskData, SpikeTools, InformationTools, JLD, Distributions,
      DataFrames, DataStructures, StatsBase, HypothesisTests

"""
    mapcells(cellfunc, sessfunc, name, [extradict, [reduce]])

The main analysis driver routine. Does the following things:

1. Calls `sessresult = sessfunc(session::Session)` on each session.
2. Calls
   `cellfunc(session::Session, sessresult, spikes::Vector{Float64}, unit::ClusterMetadata)`
   for each unit. If this method returns `nothing`, the unit is skipped. Otherwise it
   must return a Dict of data corresponding to the given unit. The keys of the Dict must
   be the same for all units.
3. For each key of the Dict returned by `cellfunc`, concatenates values across units into
   vectors.
4. Saves vectors in file `name.jld` named by key, along with `chs::Vector{ClusterMetadata}`
   (the metadata objects for analyzed units), and contents of `extradict` if provided.

If a `reduce` Dict is provided, for keys present in `reduce`, values are reduced across
units by calling the value (a function) for that key, rather than concatenated.
"""
function mapcells{F1,F2}(cellfunc::F1, sessfunc::F2, name::Symbol,
                         extradict::Union{Associative,Void}=nothing,
                         reduce::Dict{Symbol}=Dict{Symbol,Any}())
    sess = sessionnames()
    res = pmap(sess) do s
        sdata = readsession(s)
        sessres = sessfunc(sdata)
        spikes = sdata.spikes
        chs = ClusterMetadata[]
        dicts = Dict{Symbol}[]
        local reduceddict = Dict{Symbol,Any}()

        for ich = 1:length(spikes)
            isdefined(spikes, ich) || continue
            clusters = spikes[ich]
            for icluster = 1:length(clusters)
                sdata.unit_type[ich][icluster] == :A && continue
                ch = ClusterMetadata(sdata, ich, icluster)
                println(ch)
                local res = cellfunc(sdata, sessres, clusters[icluster], ch)
                isa(res, Void) && continue
                push!(chs, ch)
                push!(dicts, res)
                if !isempty(reduce)
                    for (k, v) in reduce
                        reduceddict[k] = haskey(reduceddict, k) ? v(reduceddict[k], res[k]) : res[k]
                        delete!(res, k)
                    end
                end
            end
        end
        (chs, dicts, reduceddict)
    end
    filter!(x->!isempty(x[1]), res)

    dictout::Dict{Symbol} = Dict(k=>typeof(v)[] for (k, v) in res[1][2][1])
    reduceddict = Dict{Symbol,Any}()
    for x in res
        if !isempty(reduce)
            for (k, v) in reduce
                reduceddict[k] = haskey(reduceddict, k) ? v(reduceddict[k], x[3][k]) : x[3][k]
            end
        end
        for d::Dict{Symbol} in x[2]
            for k::Symbol in keys(d)
                push!(dictout[k], d[k])
            end
        end
    end
    !isempty(reduce) && (dictout = merge(dictout, reduceddict))

    j = jldopen(joinpath(PROCESSED_DATA_DIR, "$(name).jld"), "w")
    try
        write(j, "chs", vcat([x[1] for x in res]...))
        for (k, v) = dictout
            write(j, string(k), v)
        end
        if !isa(extradict, Void)
            for (k, v) = extradict
                write(j, string(k), v)
            end
        end
    finally
        close(j)
    end
end

"""
    compute_permpval()

Computes p-values for selectivity at image presentation, using Poisson likelihood ratio
permutation test.
"""
function compute_permpval()
    NPERM = 10000
    BINS = 100:10:1000
    function sessfunc(session::Session)
        trials = stimulusframe(session.trials)
        onsets = trials[:onset]
        df = DataFrame(y=pool(trials[:stimulus]), onset=trials[:stimulus])
        dmat = ModelMatrix(ModelFrame(onset ~ y, df)).m
        perms = makeperms(size(df, 1), NPERM)

        stimulus = convert(Vector{Int}, trials[:stimulus])
        inds = vcat([find(stimulus .== i) for i = 1:maximum(stimulus)]...)
        ningroup = [sum(stimulus .== i) for i = 1:maximum(stimulus)]
        (onsets[inds], perms, ningroup)
    end

    mapcells(sessfunc, :permpval, Dict(:nperm => NPERM, :bins => BINS)) do session, x, spikes, ch
        (onsets, perms, ningroup) = x
        rast = raster(spikes, onsets, BINS/1000)'
        @time pval = poisson_deviance_grid_pval(rast, ningroup, perms)
        Dict(:pval => pval)
    end
end

"""
    compute_permpval_eachisi()

Computes p-values for selectivity at image presentation separately for each ISI, using
Poisson likelihood ratio permutation test.
"""
function compute_permpval_eachisi()
    NPERM = 10000
    BINS = 100:10:1000
    function sessfunc(session::Session)
        trials = stimulusframe(session.trials)
        onsets = trials[:onset].data
        isi = convert(Vector{Float64}, trials[:isi])
        perms = makeperms(sum(isi .== ISI[1]), NPERM)
        stimulus = convert(Vector{Int}, trials[:stimulus])
        inds = [vcat([find((stimulus .== i) & (isi .== x)) for i = 1:maximum(stimulus)]...) for x in ISI]
        ningroup = [[sum((stimulus .== i) & (isi .== x)) for i = 1:maximum(stimulus)] for x in ISI]
        (onsets, perms, inds, ningroup)
    end

    mapcells(sessfunc, :permpval_eachisi, Dict(:nperm => NPERM, :bins => BINS)) do session, x, spikes, ch
        ch.unit_type != :A || return
        (onsets, perms, inds, ningroup) = x
        rast = raster(spikes, onsets, BINS/1000)
        pval = zeros(length(ISI))
        for iisi = 1:length(ISI)
            @time pval[iisi] = poisson_deviance_grid_pval(rast[:, inds[iisi]]', ningroup[iisi], perms)
        end
        Dict(:pval => pval)
    end
end

"""
    compute_latency()

Computes latencies and durations of spiking responses by finding the window yielding the
maximal Poisson likelihood ratio.
"""
function compute_latency()
    LASTPT = 1500
    LASTPT_OPTWIN = 1000
    BASELINE_EPOCH = -1000:1000:0
    function sessfunc(session::Session)
        trials = stimulusframe(session.trials)
        @assert sort(unique(trials[:isi])) == ISI
        isitrials = [find(trials[:isi] .== isi) for isi in ISI]
        onset = trials[:onset].data
        offset = trials[:offset].data
        pos = convert(Vector{Int}, trials[:pos])
        df = DataFrame(y=pool(trials[:stimulus]), onset=onset)
        dmat = ModelMatrix(ModelFrame(onset ~ y, df)).m
        stimulus = convert(Vector{Int}, trials[:stimulus])
        trial_start_time = [session.trials[:onsets][i][1] for i = 1:size(session.trials, 1)]
        inds = vcat([find(stimulus .== i) for i = 1:maximum(stimulus)]...)
        ningroup = [sum(stimulus .== i) for i = 1:maximum(stimulus)]
        (isitrials, onset, offset, stimulus, inds, ningroup, trial_start_time, pos)
    end

    mapcells(sessfunc, :latency, Dict(:lastpt_optwin => LASTPT_OPTWIN, :lastpt => LASTPT)) do session, x, spikes, ch
        (isitrials, onset, offset, stimulus, inds, ningroup, trial_start_time, pos) = x
        nstim = maximum(stimulus)

        rast = raster(spikes, onset, (1:LASTPT)/1000)
        eta2_grid = poisson_deviance_grid(rast[1:LASTPT_OPTWIN, inds]', ningroup)
        ends, starts = ind2sub(size(eta2_grid), find(eta2_grid .== maximum(eta2_grid)))
        fr_optimal_window = zeros(nstim)
        if isempty(ends)
            stim_window = (-1, -1)
        else
            best = indmin(ends - starts)
            stim_window = (starts[best], ends[best]+1)
        end

        isi_windows = fill((-1, -1), length(isitrials))
        isi_windows_notfirst = fill((-1, -1), length(isitrials))

        for iisi = 1:length(isitrials)
            sortedisitrials = Int[]
            sortedisitrials_notfirst = Int[]
            ningroupisi = Int[]
            ningroupisi_notfirst = Int[]
            for istim = 1:nstim
                stimsel = isitrials[iisi][stimulus[isitrials[iisi]] .== istim]
                append!(sortedisitrials, stimsel)
                push!(ningroupisi, length(stimsel))
                stimsel_notfirst = isitrials[iisi][(stimulus[isitrials[iisi]] .== istim) & (pos[isitrials[iisi]] .!= 1)]
                append!(sortedisitrials_notfirst, stimsel_notfirst)
                push!(ningroupisi_notfirst, length(stimsel_notfirst))
            end

            eta2_grid = poisson_deviance_grid(rast[:, sortedisitrials]', ningroupisi)
            ends, starts = ind2sub(size(eta2_grid), find(eta2_grid .== maximum(eta2_grid)))
            if !isempty(ends)
                best = indmin(ends - starts)
                isi_windows[iisi] = (starts[best], ends[best]+1)
            end

            eta2_grid = poisson_deviance_grid(rast[:, sortedisitrials_notfirst]', ningroupisi_notfirst)
            ends, starts = ind2sub(size(eta2_grid), find(eta2_grid .== maximum(eta2_grid)))
            if !isempty(ends)
                best = indmin(ends - starts)
                isi_windows_notfirst[iisi] = (starts[best], ends[best]+1)
            end
        end

        bsdata = raster(spikes, trial_start_time, BASELINE_EPOCH/1000)
        bsscale = 1000/(last(BASELINE_EPOCH)-first(BASELINE_EPOCH))

        Dict(:window => stim_window,
             :isi_windows => isi_windows,
             :isi_windows_notfirst => isi_windows_notfirst,
             :baseline_mean => mean(bsdata)*bsscale,
             :baseline_std => std(bsdata)*bsscale
            )
    end
end

"""
    compute_alignedinfoisi()

Computes ω^2 for stimulus and interaction in different windows, along with the sum of
permutations across units.
"""
function compute_alignedinfoisi(; pos1::Bool=false)
    @load joinpath(PROCESSED_DATA_DIR, "latency.jld") window chs
    SMOOTHBY = 200
    NPERM = 1000
    PLOT_EPOCH = (-200, 1001)
    selch = allisisigunits()
    function sessfunc(session::Session)
        trials = stimulusframe(session.trials)
        nstim = maximum(trials[:stimulus])::Int

        trial_number = convert(Vector{Int}, trials[:trial_number])
        stimulus = convert(Vector{Int}, trials[:stimulus])
        dmat = zeros(length(stimulus), nstim)
        prevstim_dmat = zeros(length(stimulus), nstim)
        nextstim_dmat = zeros(length(stimulus), nstim)
        prevstim = zeros(length(stimulus))
        nextstim = zeros(length(stimulus))
        for i = 1:length(stimulus)
            if i+1 <= length(stimulus) && trial_number[i+1] == trial_number[i]
                nextstim_dmat[i, stimulus[i+1]] = 1
                nextstim[i] = stimulus[i+1]
            end
            if i-1 > 0 && trial_number[i-1] == trial_number[i]
                prevstim[i] = stimulus[i-1]
                prevstim_dmat[i, stimulus[i-1]] = 1
            end
            dmat[i, stimulus[i]] = 1
        end
        if pos1
            keep = trials[:pos] .== 1
            prevstim_dmat = zeros(sum(keep), 0)
            stimulus = stimulus[keep]
            dmat = dmat[keep, :]
            nextstim_dmat = nextstim_dmat[keep, 2:end]
            trials = trials[keep, :]
        end

        isitrials = Vector{Vector{Int}}(length(ISI))
        isiningroup = Vector{Vector{Int}}(length(ISI))
        for (iisi, isi) = enumerate(ISI)
            isitrials[iisi] = curisitrials = Int[]
            isiningroup[iisi] = curningroup = zeros(Int, nstim)
            for igrp = 1:nstim
                idx = find((trials[:isi] .== isi) & (trials[:stimulus] .== igrp))
                append!(curisitrials, idx)
                curningroup[igrp] = length(idx)
            end
        end

        isis = Vector{Float64}(trials[:isi])
        isiperms = fill(makeperms(length(isitrials[1]), NPERM), length(ISI))

        isi2perms = [begin
            perms = zeros(Int, length(isitrials[iisi])*2, NPERM)
            isistim = stimulus[(isis .== ISI[iisi]) | (isis .== ISI[iisi+1])]
            for stim = 1:nstim
                idx = find(isistim .== stim)
                shuffidx = copy(idx)
                for ibs = 1:NPERM
                    perms[idx, ibs] = shuffle!(shuffidx)
                end
            end
            perms
        end for iisi in 1:length(ISI)-1]

        (trials[:onset].data, isitrials, isiningroup, isis, dmat, prevstim_dmat, nextstim_dmat, isiperms, isi2perms)
    end

    mapcells(sessfunc, Symbol("alignedinfoisi$(pos1 ? "_pos1" : "")"),
             Dict(:fwhm => SMOOTHBY, :plot_epoch => PLOT_EPOCH, :nperm => NPERM),
             Dict(:alignedinfoperms => (a, b)->broadcast!(+, a, a, b),
                  :alignedinfoisidiffperms => (a, b)->broadcast!(+, a, a, b))) do session, x, spikes, ch
        ch in selch || return nothing
        (onset, isitrials, isiningroup, isi, dmat, prevstim_dmat, nextstim_dmat, isiperms, isi2perms) = x
        win = window[findfirst(chs, ch)]

        extract_times, extract_range, times = paddedbins(PLOT_EPOCH, 1, 5*SMOOTHBY)
        rast = boxcarsmooth(raster(spikes, onset, extract_times+win[1]/1000), SMOOTHBY, Int)[extract_range, :]

        all_dmat = [dmat nextstim_dmat prevstim_dmat]
        null_dmat = [ones(size(dmat, 1)) nextstim_dmat prevstim_dmat]
        alignedinfo = zeros(size(rast, 1), length(isitrials))
        alignedinfoperms = zeros(size(rast, 1), NPERM, length(isitrials))
        alignedinfoisidiff = zeros(size(rast, 1), length(isitrials)-1)
        alignedinfoisidiffperms = zeros(size(rast, 1), NPERM, length(isitrials)-1)
        @time for iisi = 1:length(isitrials)
            sel = isi .== ISI[iisi]
            isinulldmat = null_dmat[sel, :]
            isinulldmatqr = qrfactchk(isinulldmat)
            isidmat = all_dmat[sel, :]
            isidmatqr = qrfactchk(isidmat)
            isirast = rast[:, sel]'
            fullrss = computessr(isidmatqr, isirast)
            reducedrss = computessr(isinulldmatqr, isirast)
            fulldfe = size(isirast, 1) - size(all_dmat, 2)
            reduceddfe = size(isirast, 1) - size(isinulldmat, 2)
            omega2!(view(alignedinfo, :, iisi), fullrss, fulldfe, reducedrss, reduceddfe, size(isirast, 1))

            permisirast = zeros(Float64, size(isirast))
            isirasttmp = zeros(Float64, size(isirast))
            @assert size(isiperms[iisi], 1) == size(isirast, 1)
            for iperm = 1:NPERM
                copy!(permisirast, view(isirast, view(isiperms[iisi], :, iperm), :))
                computessr!(fullrss, isidmatqr, copy!(isirasttmp, permisirast))
                computessr!(reducedrss, isinulldmatqr, permisirast)
                omega2!(view(alignedinfoperms, :, iperm, iisi), fullrss, fulldfe, reducedrss, reduceddfe, size(isirast, 1))
            end

            if iisi <= length(isitrials) - 1
                sel = find((isi .== ISI[iisi]) | (isi .== ISI[iisi+1]))
                isirast = rast[:, sel]'
                nulldmat = zeros(length(sel), (size(null_dmat, 2)-1)*2)
                dmat2 = zeros(length(sel), size(dmat, 2)*2)
                for ipt = 1:length(sel)
                    if isi[sel[ipt]] == ISI[iisi]
                        nulldmat[ipt, 1:size(null_dmat, 2)-1] = view(null_dmat, sel[ipt], 2:size(null_dmat, 2))
                        dmat2[ipt, 1:size(dmat, 2)] = view(dmat, sel[ipt], :)
                    else
                        nulldmat[ipt, size(null_dmat, 2):end] = view(null_dmat, sel[ipt], 2:size(null_dmat, 2))
                        dmat2[ipt, size(dmat, 2)+1:end] = view(dmat, sel[ipt], :)
                    end
                end
                nulldmatqr = qrfactchk([isi[sel] .== ISI[iisi] nulldmat dmat[sel, :]])
                dmatqr = qrfactchk([nulldmat dmat2])
                permisirast = zeros(Float64, size(isirast))
                isirasttmp = zeros(Float64, size(isirast))
                reduceddfe2 = size(isirast, 1) - size(nulldmatqr, 2)
                fulldfe2 = size(isirast, 1) - size(dmatqr, 2)
                computessr!(fullrss, dmatqr, copy!(isirasttmp, isirast))
                computessr!(reducedrss, nulldmatqr, copy!(isirasttmp, isirast))
                omega2!(view(alignedinfoisidiff, :, iisi), fullrss, fulldfe2, reducedrss, reduceddfe2, size(isirast, 1))
                @assert size(isi2perms[iisi], 1) == size(isirast, 1)
                for iperm = 1:NPERM
                    copy!(permisirast, view(isirast, view(isi2perms[iisi], :, iperm), :))
                    computessr!(fullrss, dmatqr, copy!(isirasttmp, permisirast))
                    computessr!(reducedrss, nulldmatqr, permisirast)
                    omega2!(view(alignedinfoisidiffperms, :, iperm, iisi), fullrss, fulldfe2, reducedrss, reduceddfe2, size(isirast, 1))
                end
            end
        end
        alignedinfo[isnan(alignedinfo)] = 0
        alignedinfoperms[isnan(alignedinfoperms)] = 0
        alignedinfoisidiff[isnan(alignedinfoisidiff)] = 0
        alignedinfoisidiffperms[isnan(alignedinfoisidiffperms)] = 0

        Dict(:window => win,
             :alignedinfo => alignedinfo,
             :alignedinfoperms => alignedinfoperms,
             :alignedinfoisidiff => alignedinfoisidiff,
             :alignedinfoisidiffperms => alignedinfoisidiffperms)
    end
end

"""
    shuffleby(x::AbstractVector, NPERM)::Matrix{Int}

Create a `length(x)` × `NPERM` matrix of indices for permutations of `x`, permuting within
levels of `x`.
"""
function shuffleby(x::AbstractVector, NPERM)
    perms = zeros(Int, length(x), NPERM)
    for val in unique(x)
        inds = find(x .== val)
        shuffled = copy(inds)
        for i = 1:NPERM
            perms[inds, i] = shuffle!(shuffled)
        end
    end
    perms
end

"""
    rearrange(x::AbstractVector)

Returns indices that rearrange `x` so that levels are adjacent, as well as counts for each
group. 
"""
function rearrange(x::AbstractVector)
    uniques = sort!(unique(x))
    inds = vcat([find(x .== i) for i = uniques]...)
    ningroup = [sum(x .== i) for i = uniques]
    (inds, ningroup)
end

"""
    compute_permpval_isiwin()

Performs permutation test for image selectivity in different bins, for different groups of
ISIs.
"""
function compute_permpval_isiwin()
    selch = permselunits()
    @load joinpath(PROCESSED_DATA_DIR, "latency.jld") chs window
    NPERM = 10000
    WIN = [(0, 200), (200, 400), (400, 600), (700, 900)]
    function sessfunc(session::Session)
        trials = stimulusframe(session.trials)
        isi = convert(Vector{Float64}, trials[:isi])
        stimulus = convert(Vector{Float64}, trials[:stimulus])
        tgrp = Vector{Vector{Int}}(length(ISI))
        ningroup = Vector{Vector{Int}}(length(ISI))
        for i = 1:length(ISI)
            idx = find(isi .>= ISI[i])
            inds, ningroup[i] = rearrange(stimulus[idx])
            tgrp[i] = idx[inds]
        end
        next_stimulus = convert(Vector{Int}, trials[:next_stimulus])
        perms = [shuffleby(next_stimulus[x], NPERM) for x in tgrp]
        (trials[:onset].data, tgrp, ningroup, perms)
    end

    mapcells(sessfunc, :permpval_isiwin, Dict(:nperm => NPERM, :win => WIN, :minisi => ISI)) do session, x, spikes, ch
        ch in selch || return nothing
        stim_window = window[findfirst(chs, ch)]

        (onset, tgrp, ningroup, perms) = x
        pval_cumulative = ones(length(ISI), length(WIN))
        @time for (i, win) = enumerate(WIN)
            rast = vec(raster(spikes, onset+stim_window[1]/1000, (win[1]:win[2]-win[1]:win[2])/1000))
            for j = 1:length(ningroup)
                y = rast[tgrp[j]]
                pval_cumulative[j, i] = poisson_deviance_pval(y, ningroup[j], perms[j])
            end
        end
        Dict(
            :pval_cumulative => pval_cumulative
        )
    end
end

"""
    compute_permpval_maintenance()

Performs permutation test for maintenance period stimulus selectivity.
"""
function compute_permpval_maintenance()
    NPERM = 100000
    WIN = (300, 2400)
    function sessfunc(session::Session)
        trials = session.trials
        nstim = maximum(map(maximum, trials[:stimuli]))
        X = zeros(size(trials, 1), nstim)
        for i = 1:size(trials, 1)
            X[i, trials[i, :stimuli]::Vector{Int}] = 1
        end
        perms = makeperms(size(X, 1), NPERM)
        (trials[:maint].data, qrfactchk(X), perms)
    end

    mapcells(sessfunc, :permpval_maintenance, Dict(:nperm => NPERM, :win => WIN)) do session, x, spikes, ch
        (maint, qrf, perms) = x
        null = qrfact(ones(size(qrf, 1), 1))
        rast = map(Float64, raster(spikes, maint, (WIN[1]:WIN[2]-WIN[1]:WIN[2])/1000)')
        trueval = calculate_f(rast, qrf, null)[1]
        if isnan(trueval)
            pval = NaN
        else
            permvals = calculate_f(rast[perms], qrf, null)
            pval = mean(x->x >= trueval, permvals)
        end
        Dict(
            :pval => pval
        )
    end
end

"""
    compute_modulation_maintenance()

Computes ω^2 at maintenance, as well as mean responses by stimulus.
"""
function compute_modulation_maintenance()
    WIN = (300, 2400)
    function sessfunc(session::Session)
        trials = session.trials
        nstim = maximum(map(maximum, trials[:stimuli]))
        X = zeros(size(trials, 1), nstim)
        for i = 1:size(trials, 1)
            X[i, trials[i, :stimuli]::Vector{Int}] = 1
        end
        stimuli = trials[:stimuli]
        (trials[:maint].data, stimuli.data, nstim, qrfactchk(X), qrfactchk(ones(size(X, 1), 1)))
    end

    mapcells(sessfunc, :modulation_maintenance, Dict(:win => WIN)) do session, x, spikes, ch
        (maint, stimuli, nstim, X, null) = x
        rast = map(Float64, raster(spikes, maint, (WIN[1]:WIN[2]-WIN[1]:WIN[2])/1000)')
        means = [mean(rast[[i for i = 1:length(rast) if stim in stimuli[i]]]) for stim in 1:nstim]
        Dict(:pve => calculate_omega2(rast, X, null)[1],
             :means => means)
    end
end

"""
    compute_modulation_presentation()

Computes mean responses by stimulus at presentation, as well as baseline firing rate and
standard deviation.
"""
function compute_modulation_presentation()
    @load joinpath(PROCESSED_DATA_DIR, "latency.jld") chs window
    chs = chs::Vector{ClusterMetadata}
    window = window::Vector{Tuple{Int,Int}}

    function sessfunc(session::Session)
        trials = stimulusframe(session.trials)
        @assert sort(unique(trials[:isi])) == ISI
        isi = convert(Vector{Float64}, trials[:isi])
        onset = trials[:onset].data
        stimulus = convert(Vector{Int}, trials[:stimulus])::Vector{Int}
        trial_start_time = [session.trials[:onsets][i][1] for i = 1:size(session.trials, 1)]
        (isi, onset, stimulus, trial_start_time)
    end

    mapcells(sessfunc, :modulation_presentation) do session, x, spikes, ch
        (isi, onset, stimulus, trial_start_time) = x
        chind = findfirst(chs, ch)
        win = window[chind]
        win_rg = (win[1]:win[2]-win[1]:win[2])/1000
        @assert length(win_rg) == 2
        nstim = maximum(stimulus)
        stim_fr_win = zeros(nstim)
        for stim = 1:nstim
            stim_onset = onset[stimulus .== stim]
            stim_fr_win[stim] = mean(raster(spikes, stim_onset, win_rg))
        end
        baseline_fr_all, baseline_std_all = mean_and_std(raster(spikes, trial_start_time, -1:0))
        baseline_fr_win, baseline_std_win = mean_and_std(raster(spikes, trial_start_time, win_rg-1))
        Dict(:stim_fr_win => stim_fr_win,
             :baseline_fr_all => baseline_fr_all,
             :baseline_std_all => baseline_std_all,
             :baseline_fr_win => baseline_fr_win,
             :baseline_std_win => baseline_std_win,
             :win => win)
    end
end

"""
    compute_maintenance_topstim_holdout()

Computes normalized modulation for best stimulus at presentation during maintenance,
determining best stimulus and maintenance period modulation on different subsets of trials.
"""
function compute_maintenance_topstim_holdout()
    WIN = (300, 2400)
    BSWIN = (300, 1000)
    LASTPT = 1000
    @load joinpath(PROCESSED_DATA_DIR, "modulation_presentation.jld") chs baseline_std_all

    function sessfunc(session::Session)
        s1 = 1:4:size(session.trials, 1)
        s2 = [!(x in s1) for x in 1:size(session.trials, 1)]
        strials = stimulusframe(session.trials[s1, :])
        onset = strials[:onset].data
        stimulus = convert(Vector{Int}, strials[:stimulus])::Vector{Int}
        trial_start_time = [x[1] for x in session.trials[s1, :onsets]]
        inds = vcat([find(stimulus .== i) for i = 1:maximum(stimulus)]...)
        ningroup = [sum(stimulus .== i) for i = 1:maximum(stimulus)]

        trials = session.trials[s2, :]
        (onset, stimulus, inds, ningroup, trial_start_time, getindex.(trials[:onsets].data, 1), trials[:maint].data, trials[:stimuli].data)
    end

    mapcells(sessfunc, :maintenance_topstim_holdout, Dict(:win => WIN)) do session, x, spikes, ch
        chind = findfirst(chs, ch)
        chind == 0 && return
        (ponsets, pstimulus, inds, ningroup, pbaseline, mbaseline, maint, stimuli) = x

        rast = raster(spikes, ponsets, (1:LASTPT)/1000)
        eta2_grid = poisson_deviance_grid(rast[:, inds]', ningroup)
        ends, starts = ind2sub(size(eta2_grid), find(eta2_grid .== maximum(eta2_grid)))
        isempty(ends) && return
        best = indmin(ends - starts)
        pwin = (starts[best], ends[best]+1)

        rast = raster(spikes, ponsets, (pwin[1]:pwin[2]-pwin[1]:pwin[2])/1000)
        stimresp = [mean(rast[pstimulus .== i]) for i = 1:maximum(pstimulus)]
        baseline = mean(raster(spikes, pbaseline, (pwin[1]:pwin[2]-pwin[1]:pwin[2])/1000-1))
        beststim = indmax(abs(stimresp - baseline))
        ispositive = indmax(stimresp) == beststim

        rast = raster(spikes, maint, (WIN[1]:WIN[2]-WIN[1]:WIN[2])/1000)'
        g1 = find(x->beststim in x, stimuli)
        g2 = find(x->!(beststim in x), stimuli)
        pval_mwu = pvalue(MannWhitneyUTest(rast[g1], rast[g2]))
        mbsrast = raster(spikes, maint, (BSWIN[1]:BSWIN[2]-BSWIN[1]:BSWIN[2])/1000)'
        bsrast = raster(spikes, [pbaseline; mbaseline], (BSWIN[1]:BSWIN[2]-BSWIN[1]:BSWIN[2])/1000-1)'
        sdest = baseline_std_all[chind]*(WIN[2]-WIN[1])/1000

        mg2 = mean(rast[g2])
        sep = (mean(rast[g1]) - mg2)/sdest*ifelse(ispositive, 1, -1)
        sep_by_pos = [mean(rast[find(x->x[pos] == beststim, stimuli)]) - mg2 for pos = 1:length(stimuli[1])]/sdest*ifelse(ispositive, 1, -1)
        Dict(
            :pwin => pwin,
            :beststim => beststim,
            :ng1 => length(g1),
            :ng2 => length(g2),
            :ispositive => ispositive,
            :pval_mwu => pval_mwu,
            :sep => sep,
            :sep_by_pos => sep_by_pos
        )
    end
end

"""
    compute_maintenance_topstim_correct()

Computes normalized modulation for best stimulus at presentation during maintenance on
correct trials relative to incorrect trials, when the stimulus was probed.
"""
function compute_maintenance_topstim_correct()
    individual_unit_test(rast, g1, g2) =
        isempty(g1) || isempty(g2) ? NaN : pvalue(MannWhitneyUTest(rast[g1], rast[g2]))

    modulation(rast, g1, g2, sdest, ispositive) =
        ifelse(ispositive, 1, -1)*(mean(rast[g1]) - mean(rast[g2]))/sdest

    WIN = (300, 2400)
    EARLY_WIN = (0, 200)
    PERSIST_WIN = (300, 700)

    sel = permselunits()
    @load joinpath(PROCESSED_DATA_DIR, "modulation_presentation.jld") chs stim_fr_win baseline_fr_all baseline_fr_win win baseline_std_all
    stim_chs = chs
    ispositive = [maxabs(fr-bs) == maximum(fr-bs) for (fr, bs) in zip(stim_fr_win, baseline_fr_win)]
    beststim = [indmax(ifelse(ip, 1, -1)*x) for (x, ip) in zip(stim_fr_win, ispositive)]

    function sessfunc(session::Session)
        trials = session.trials
        (trials[:onsets].data, trials[:maint].data, trials[:probes].data, trials[:stimuli].data, trials[:probe1].data, trials[:probe2].data,
         convert(Vector{Float64}, trials[:isi]), convert(Vector{Bool}, trials[:correct]))
    end

    mapcells(sessfunc, :maintenance_topstim_correct, Dict(:win => WIN)) do session, x, spikes, ch
        ch in sel || return
        chind = findfirst(chs, ch)
        chind == 0 && return
        (onsets, maint, probes, stimuli, probe1, probe2, isi, correct) = x
        bs = beststim[chind]
        intrial = [bs in x for x in stimuli]
        baselines = [x[1]-1 for x in onsets]

        rast = vec(raster(spikes, maint, (WIN[1]:WIN[2]-WIN[1]:WIN[2])/1000))
        probed = (probe1 .== bs) | (probe2 .== bs)
        g1 = find(intrial & probed & correct)
        ng1 = length(g1)
        g2 = find(intrial & probed & !correct)
        ng2 = length(g2)
        pval_mwu_maint = individual_unit_test(rast, g1, g2)
        sep_maint = modulation(rast, g1, g2, baseline_std_all[chind]*(WIN[2]-WIN[1])/1000, ispositive[chind])

        pres = [y[findfirst(x, bs)] for (y, x) in zip(onsets[intrial], stimuli[intrial])]
        cwin = win[chind]
        g1 = find(probed[intrial] & correct[intrial])
        g2 = find(probed[intrial] & !correct[intrial])

        rast = vec(raster(spikes, pres, (cwin[1]+EARLY_WIN[1]:EARLY_WIN[2]-EARLY_WIN[1]:cwin[1]+EARLY_WIN[2])/1000))
        pval_mwu_pres = individual_unit_test(rast, g1, g2)
        sep_pres = modulation(rast, g1, g2, baseline_std_all[chind]*(EARLY_WIN[2]-EARLY_WIN[1])/1000, ispositive[chind])

        use = isi[intrial] .>= .5
        rast = vec(raster(spikes, pres, (cwin[1]+PERSIST_WIN[1]:PERSIST_WIN[2]-PERSIST_WIN[1]:cwin[1]+PERSIST_WIN[2])/1000))
        sdest = baseline_std_all[chind]*(PERSIST_WIN[2]-PERSIST_WIN[1])/1000
        g1 = find(probed[intrial] & correct[intrial] & use)
        g2 = find(probed[intrial] & !correct[intrial] & use)
        ng1_persist = length(g1)
        ng2_persist = length(g2)
        pval_mwu_persist = individual_unit_test(rast, g1, g2)
        sep_persist = modulation(rast, g1, g2, baseline_std_all[chind]*(PERSIST_WIN[2]-PERSIST_WIN[1])/1000, ispositive[chind])

        Dict(
            :ng1 => ng1,
            :ng2 => ng2,
            :ng1_persist => ng1_persist,
            :ng2_persist => ng2_persist,
            :ispositive => ispositive[chind],
            :pval_mwu_maint => pval_mwu_maint,
            :pval_mwu_pres => pval_mwu_pres,
            :pval_mwu_persist => pval_mwu_persist,
            :sep_maint => sep_maint,
            :sep_persist => sep_persist,
            :sep_pres => sep_pres
        )
    end
end

"""
    bs95ci(X)

Bootstrap 95% confidence interval; trials along first dimension of X
"""
function bs95ci(X::AbstractArray, bs=rand(1:size(X, 1), size(X, 1), 1000))
    @assert size(X, 1) == size(bs, 1)
    tmp = zeros(size(bs, 2))
    ci = zeros(2, size(X, 2))
    for ip = 1:size(X, 2)
        for i = 1:size(bs, 2)
            mu = 0.0
            @inbounds @simd for j = 1:size(bs, 1)
                mu += X[bs[j, i], ip]
            end
            mu /= size(X, 1)
            tmp[i] = mu
        end
        ci[:, ip] = quantile!(tmp, [0.025, 0.975])
    end
    ci
end

"""
    bsvs0(X)

Single-sample bootstrap test that X is > 0
"""
function bsvs0(X::AbstractVector, nbs=10000)
    n = 0
    gen = Base.Random.RangeGenerator(1:size(X, 1))
    generator = Base.Random.GLOBAL_RNG
    for i = 1:nbs
        mu = 0.0
        @inbounds for j = 1:size(X, 1)
            mu += X[rand(generator, gen)]
        end
        n += mu <= 0
    end
    n/nbs
end

"""
    permselunits(X)::Vector{ClusterMetadata}

ClusterMetadata for units significant at presentation.
"""
function permselunits(; alpha=0.001)
    @load joinpath(PROCESSED_DATA_DIR, "permpval.jld") pval chs
    chsel = find(x->(x.unit_type == :MU || x.unit_type == :SU), chs)
    selchs = chs[chsel][pval[chsel] .< alpha]
end

"""
    allisisigunits(X)::Vector{ClusterMetadata}

ClusterMetadata for units significant at presentation separately at each ISI.
"""
function allisisigunits(; alpha=0.05)
    @load joinpath(PROCESSED_DATA_DIR, "permpval_eachisi.jld") chs pval
    chsel = find(x->(x.unit_type == :MU || x.unit_type == :SU), chs)
    chs[chsel[[all(x->x < alpha, x) for x in pval[chsel]]]]
end

"""
    maintselunits(X)::Vector{ClusterMetadata}

ClusterMetadata for units significant at maintenance.
"""
function maintselunits(; alpha=0.01)
    @load joinpath(PROCESSED_DATA_DIR, "permpval_maintenance.jld") pval chs
    chsel = find(x->(x.unit_type == :MU || x.unit_type == :SU), chs)
    selchs = chs[chsel][pval[chsel] .< alpha]
end
end