module Summary
using Sessions, TaskData, Analysis, SpikeTools, PlotTools, PyPlot, PyCall, JLD,
      HypothesisTests, Distributions, DataFrames, DataStructures, StatsBase, RCall
@pyimport matplotlib.gridspec as gridspec
@pyimport matplotlib.patches as patches
const slice = pybuiltin(:slice)

## Results

unit_mask(chs, region=:MTL, su_only::Bool=false) =
    [(ch.unit_type == :SU || (!su_only && ch.unit_type == :MU)) && ismatch(REGIONS[region], ch.site)
     for ch in chs]

function to_pct(x::AbstractArray{Bool})
    @sprintf "%.2g%% (%d/%d)" mean(x)*100 sum(x) length(x)
end

function region_breakdown(chs, sig::AbstractArray{Bool})
    println("    $(to_pct(sig[[ismatch(REGIONS[:PHG], ch.site) for ch in chs]])) Parahippocampal Units")
    println("    $(to_pct(sig[[ismatch(REGIONS[:EC], ch.site) for ch in chs]])) Entorhinal Units")
    println("    $(to_pct(sig[[ismatch(REGIONS[:H], ch.site) for ch in chs]])) Hippocampal Units")
    println("    $(to_pct(sig[[ismatch(REGIONS[:A], ch.site) for ch in chs]])) Amygdala Units")
end

function stimsel_info()
    regions = [:PHG, :EC, :H, :A]
    @load joinpath(PROCESSED_DATA_DIR, "permpval.jld") chs pval
    selch = unit_mask(chs, :MTL)
    chs = chs[selch]
    pval = pval[selch]
    sig = pval .< 0.001

    println("$(to_pct(sig)) units responsive")
    region_breakdown(chs, sig)
    proportions = [begin
        maskregion = unit_mask(chs, region)
        (sum(sig[maskregion]), sum(maskregion))
    end for region in regions]
    for ir1 = 1:length(regions), ir2 = ir1+1:length(regions)
        p1 = proportions[ir1]
        p2 = proportions[ir2]
        @printf "%3s vs. %3s p = %.3g\n" regions[ir1] regions[ir2] pvalue(FisherExactTest(p1[1], p1[2]-p1[1], p2[1], p2[2]-p2[1]))
    end
end

function showstimulus(f)
    text(0.5, 0.5, f, horizontalalignment="center", verticalalignment="center", fontsize=9)
    axis("off")
end

function plot_stimulus_all(epoch, session, channel_number, cluster_number; header::Bool=true, highlight=nothing, out="pdf")
    const LIMS = epoch == :maintenance ? (-1000, 2400) : (-200, 800)
    spikes = session.spikes[channel_number][cluster_number]
    sesstrials = session.trials
    trials = stimulusframe(sesstrials)
    nstim = length(session.stimuli)
    ticks = epoch == :maintenance ? (0:1000:2000) : 0:600:1000
    gcf()[:set_size_inches](4.488, 2.1)
    if header
        suptitle("$(session.folder) channel $(channel_number) cluster $(cluster_number) $(session.sites[channel_number]) $epoch")
        top = 0.93
    else
        top = 0.995
    end
    gs = gridspec.GridSpec(1, nstim, left=0.06, right=0.98, bottom=0.135, top=top, wspace=0.2, hspace=0.0)
    local ax
    nthimg = 4
    nthraster = 3
    rasterbin = epoch == :maintenance ? 200 : 50
    if isa(highlight, Void)
        highlighttrials = Int[]
    else
        stims = sesstrials[:stimuli].data::Vector{Vector{Int}}
        highlighttrials = find(Bool[highlight in x for x in stims])
    end

    for stim = 1:nstim
        gs_item = gs[:__getitem__]((0, stim-1))
        inner_gs = gridspec.GridSpecFromSubplotSpec(nthimg, 1, subplot_spec=gs_item, hspace=0.02)
        subplot(inner_gs[:__getitem__](0))
        showstimulus(stim)

        gs_item = inner_gs[:__getitem__](slice(1, nthimg))
        inner_gs = gridspec.GridSpecFromSubplotSpec(nthraster, 1, subplot_spec=gs_item, hspace=0.0)
        axa = subplot(inner_gs[:__getitem__]((slice(0, nthraster-1))), zorder=1000)
        local usetrials::Vector{Int}
        if epoch == :sample
            stims = sesstrials[:stimuli].data::Vector{Vector{Int}}
            usetrials = find([stim in x for x in stims])
            onsets = trials[:onset][trials[:stimulus] .== stim].data
        elseif epoch == :maintenance
            stims = sesstrials[:stimuli].data::Vector{Vector{Int}}
            usetrials = find([stim in x for x in stims])
            onsets = sesstrials[:maint][usetrials]
        end
        t1 = findin(usetrials, setdiff(usetrials, highlighttrials))
        t2 = findin(usetrials, highlighttrials)

        plotraster(spikes, onsets[[t1; t2]], LIMS)
        if !isa(highlight, Void)
            axhline(length(t1)+0.5, color=(1, 0, 0, 0.5), linewidth=0.5)
        end
        axis("off")
        axvline(0, ymin=-1/(nthraster-1), ymax=1, linestyle=":", color="darkcyan", dashes=(1, 1),
                   linewidth=0.5, clip_on=false)

        if stim == 1
            ylabel("Presentation Number")
            bax = ax = subplot(inner_gs[:__getitem__](nthraster-1))
        else
            gca()[:set_yticklabels](())
            bax = subplot(inner_gs[:__getitem__](nthraster-1), sharey=ax)
        end
        rast = raster(spikes, onsets, (LIMS[1]:rasterbin:LIMS[2])/1000)
        bar(LIMS[1]:rasterbin:LIMS[2]-rasterbin, mean(rast, 2)*1000/rasterbin, rasterbin, color="k", edgecolor="none", clip_on=false)
        if !isa(highlight, Void)
            bar(LIMS[1]:rasterbin:LIMS[2]-rasterbin, mean(rast[:, t1], 2)*1000/rasterbin,
                rasterbin, color=(1, 0, 0, 0.5), edgecolor="none", clip_on=false)
        end
        if epoch == :maintenance
            xticks(ticks, ["0"; [string(Int(x/1000)) for x in ticks[2:end]]], fontsize=7)
            stim == 1 && xlabel("Time After Mask Offset (Seconds)", fontsize=7, labelpad=2, horizontalalignment="left", x=-0.6)
        else
            xticks(ticks, fontsize=7)
            stim == 1 && xlabel("Time After Stimulus Onset (ms)", fontsize=7, labelpad=2, horizontalalignment="left", x=-0.6)
        end
        if stim == 1
            ylabel("Firing Rate (Hz)", fontsize=7, labelpad=3)
            removespines(["top", "right"])
            gca()[:yaxis][:set_ticks_position]("left")
            gca()[:xaxis][:set_ticks_position]("bottom")
        else
            bax[:axes][:get_yaxis]()[:set_visible](false)
            removespines(["left", "top", "right"])
            gca()[:xaxis][:set_ticks_position]("bottom")
        end
        tick_params(axis="y", pad=2)
        tick_params(axis="x", direction="out", length=3, pad=1)
        xlim(LIMS[1], LIMS[2])
    end
    sca(ax)
    yl = ax[:get_ylim]()
    yticks((0:max(div(yl[2], 3), 1):yl[2]), fontsize=7)
    savefig("$(session.folder)_ch$(channel_number)c$(cluster_number)_$(session.sites[channel_number])_$(session.unit_type[channel_number][cluster_number])_all_$(epoch).$(out)", dpi=300)
end

## Responses Persist Until Subsequent Image Presentations

function read_resptime(about; notfirst::Bool=false)
    @load joinpath(PROCESSED_DATA_DIR, "latency.jld") chs isi_windows isi_windows_notfirst
    if notfirst
        isi_windows = isi_windows_notfirst
    end
    sel = findin(chs, intersect(Analysis.allisisigunits()))
    if about == :latency
        x = [isi_windows[isel][iisi][1] for iisi = 1:4, isel = sel]
    elseif about == :duration
        x = [isi_windows[isel][iisi][2]-isi_windows[isel][iisi][1] for iisi = 1:4, isel = sel]
    else
        error("invalid about $about")
    end
    x, chs[sel]
end

function resptime_anova(about::Symbol=:duration; notfirst::Bool=false)
    x, chs = read_resptime(about, notfirst=notfirst)
    isi = vec(repmat([1, 2, 3, 4], 1, size(x, 2)))
    region = vec(repmat(reshape([string(siteregion(x)) for x in chs], 1, length(chs)), size(x, 1), 1))
    cell = vec(repmat(reshape(1:length(chs), 1, length(chs)), size(x, 1), 1))
    [@assert sum((isi .== y) & (cell .== x)) == 1 for x in cell, y in isi]
    println()
    println(R"""
    library(ez)
    df = data.frame(x=$(vec(x)), isi=factor($isi), region=factor($region), cell=factor($cell))
    ezANOVA(df, x, cell, within=isi, between=region, type=2, detailed=TRUE)
    """)
    for i = 1:length(ISI), j = i+1:length(ISI)
        test = OneSampleTTest(x[i, :], x[j, :])
        @printf("ISI %4d vs. %4d  t(%d) = %.1f p = %.3g\n", ISI[i]*1000, ISI[j]*1000, test.df, test.t, pvalue(test))
    end
end

function resptime_table(about::Symbol=:duration)
    x, chs = read_resptime(about)
    regions = [:PHG, :EC, :H, :A, :All]
    txt = Array(String, length(regions), length(ISI)+1)
    for i = 1:length(regions)
        masked = x[:, unit_mask(chs, regions[i])]
        for j = 1:length(ISI)
            txt[i, j] = @sprintf "%d (%d)" mean(masked[j, :]) std(masked[j, :])
        end
        txt[i, end] = @sprintf "%d (%d)" mean(mean(masked, 1)) std(mean(masked, 1))
    end
    [["" reshape(["ISI $(x)" for x in round(Int, ISI*1000)], 1, 4) "Average"]; regions txt]
end


function plot_isi_single(session, channel_number, cluster_number, stim)
    @load joinpath(PROCESSED_DATA_DIR, "latency.jld") chs isi_windows
    isiwins = isi_windows[findfirst(chs, ClusterMetadata(session, channel_number, cluster_number))]

    const LIMS = (-200, 1500)
    spikes = session.spikes[channel_number][cluster_number]
    trials = stimulusframe(session.trials)
    nstim = maximum(trials[:stimulus])::Int
    byval::Vector{Float64} = sort(unique(trials[:isi]))
    gcf()[:set_size_inches](3.34, 1.5)
    # ioff()
    gs = gridspec.GridSpec(2, 4, left=0.1, right=0.99, bottom=0.2, top=0.9, wspace=0.15, hspace=0.0)
    rasterbin = 100
    local ax
    for ibv = 1:length(byval)
        subplot(gs[:__getitem__]((0, ibv-1)), zorder=1000)
        onsets = trials[:onset][(trials[:stimulus] .== stim) & (trials[:isi] .== byval[ibv])].data
        plotraster(spikes, onsets, LIMS)
        xlim(LIMS[1], LIMS[2])
        axis("off")
        isi = Int(byval[ibv]*1000)
        title("ISI = $(isi) ms", y=0.98, fontsize=8)
        axvline(isiwins[ibv][1], ymin=-1, ymax=1, linestyle=":", color="darkcyan", dashes=(1, 1),
                linewidth=0.5, clip_on=false)
        axvline(isiwins[ibv][2], ymin=-1, ymax=1, linestyle=":", color="darkcyan", dashes=(1, 1),
                linewidth=0.5, clip_on=false)
        println("ISI = $(isi) ms onset = $(isiwins[ibv][1]) offset = $(isiwins[ibv][2]) duration = $(isiwins[ibv][2]-isiwins[ibv][1])")

        if ibv == 1
            setp(gca()[:get_yticklabels](), fontsize=9)
        else
            gca()[:set_yticklabels](())
        end

        if ibv == 1
            bax = ax = subplot(gs[:__getitem__]((1, ibv-1)))
        else
            bax = subplot(gs[:__getitem__]((1, ibv-1)), sharey=ax)
        end
        rast = raster(spikes, onsets, (LIMS[1]:rasterbin:LIMS[2])/1000)
        bar(LIMS[1]:rasterbin:LIMS[2]-rasterbin, mean(rast, 2)*1000/rasterbin, rasterbin, color="k", edgecolor="none", clip_on=false)
        xticks(0:500:1500, ("0", "", "1000", ""), fontsize=7)
        removespines(["top", "right"])
        bax[:yaxis][:set_ticks_position]("left")
        bax[:xaxis][:set_ticks_position]("bottom")
        if ibv == 1
            setp(bax[:get_yticklabels](), fontsize=7)
            setp(bax[:get_xticklabels](), fontsize=7)
            ylabel("Firing Rate (Hz)", fontsize=7, labelpad=3)
            removespines(["top", "right"])
            gca()[:yaxis][:set_ticks_position]("left")
            gca()[:xaxis][:set_ticks_position]("bottom")
            xlabel("Time After Stimulus Onset (ms)", fontsize=7, labelpad=2, horizontalalignment="left", x=-0.4)
        else
            bax[:set_xticklabels](())
            bax[:axes][:get_yaxis]()[:set_visible](false)
            removespines(["left", "top", "right"])
            gca()[:xaxis][:set_ticks_position]("bottom")
        end
        tick_params(axis="y", pad=3)
        tick_params(axis="x", direction="out", length=3, pad=1)
        xlim(LIMS[1], LIMS[2])
    end
    yl = ax[:get_ylim]()
    yticks((0:max(div(yl[2], 4), 1):yl[2]), fontsize=7)
    savefig("$(session.folder)_ch$(channel_number)c$(cluster_number)_isi_stim$stim.pdf", dpi=300)
end

function plot_resptime(about::Symbol=:duration)
    gcf()[:set_size_inches](3.34, 1.5)
    x, chs = read_resptime(about)
    scatter(repeat(ISI*1000, outer=size(x, 2))+(rand(length(x))-0.5)*100, vec(x), 20, marker=".", alpha=0.4, color="k", edgecolor="none")
    means = zeros(4)
    stds = zeros(4)
    for i = 1:4
        xi = x[i, :]
        means[i] = mean(xi)
        stds[i] = std(xi)/sqrt(length(xi))
    end
    errorbar(ISI*1000, means, stds, color="k")
    removespines()
    gca()[:yaxis][:set_ticks_position]("left")
    gca()[:xaxis][:set_ticks_position]("bottom")
    xticks(ISI*1000, fontsize=7)
    yticks(fontsize=7)
    xlim(-100, 900)
    ylim(0, 1085)
    xlabel("Inter-Stimulus Interval (ms)", fontsize=7, labelpad=3)
    ylabel("Response Duration (ms)", fontsize=7, labelpad=2, y=0.45)
    gap = 10
    ht1 = 1025
    len1 = 25
    makeline(ISI[1]*1000, ISI[2]*1000-gap, ht1, ht1, len1, pvalue(OneSampleTTest(x[1, :], x[2, :])))
    makeline(ISI[2]*1000+gap, ISI[3]*1000-gap, ht1, ht1, len1, pvalue(OneSampleTTest(x[2, :], x[3, :])))
    makeline(ISI[3]*1000+gap, ISI[4]*1000, ht1, ht1, len1, pvalue(OneSampleTTest(x[3, :], x[4, :])))
    tight_layout(rect=[-.03, -.07, 1.04, 1.0])
    savefig("$(about)_by_isi.pdf")
end

function resptime_region_differences(about::Symbol=:duration)
    x, chs = read_resptime(about)
    regions = [:PHG, :EC, :H, :A]
    region = [siteregion(ch) for ch in chs]
    regionnum = [findfirst(regions, x) for x in region]
    meandur = vec(mean(x, 1))
    means = zeros(4)
    stds = zeros(4)
    for i = 1:4
        xi = meandur[regionnum .== i]
        means[i] = mean(xi)
        stds[i] = std(xi)/sqrt(length(xi))
    end
    PHG = meandur[regionnum .== 1]
    EC = meandur[regionnum .== 2]
    H = meandur[regionnum .== 3]
    A = meandur[regionnum .== 4]
    @printf("PHG: %d ms\n", mean(PHG))
    @printf("EC: %d ms\n", mean(EC))
    @printf("H: %d ms\n", mean(H))
    @printf("A: %d ms\n", mean(A))
    for ir1 = 1:length(regions), ir2 = ir1+1:length(regions)
        md1 = meandur[regionnum .== ir1]
        md2 = meandur[regionnum .== ir2]
        uet = UnequalVarianceTTest(md1, md2)
        @printf("%3s vs. %3s t(%d) = %.2g p = %.2g\n",
            regions[ir1],
            regions[ir2],
            uet.df,
            uet.t,
            pvalue(uet))
    end
end

function plot_alignedisi(; ch=nothing, compare::Bool=false, pos1::Bool=false)
    @load(joinpath(PROCESSED_DATA_DIR, "alignedinfoisi$(pos1 ? "_pos1" : "").jld"), chs, plot_epoch, alignedinfo,
          alignedinfoperms, alignedinfoisidiff, alignedinfoisidiffperms)
    if compare
        alignedinfo = alignedinfoisidiff
        alignedinfoperms = alignedinfoisidiffperms
    end
    _, sel = read_resptime(:duration)
    @assert sel == chs
    info = cat(3, alignedinfo...)
    alignedinfoperms ./= length(chs)

    clf()
    gcf()[:set_size_inches](3.34, 2)
    win = plot_epoch
    xv = win[1]:win[2]-1
    xvsel = 1:length(xv)
    xv = xv[xvsel]
    info = info[xvsel, :, :]
    mu = mean(info, 3)
    sd = std(info, 3)
    lim = 0.0
    limmin = 0.0
    colors = ["b", "g", "r", "c"]
    compare && shift!(colors)
    for i = 1:size(info, 2)
        y = vec(mu[:, i])
        yerr = vec(sd[:, i])/sqrt(size(info, 3))
        lim = max(lim, maximum(y+yerr))
        limmin = min(limmin, minimum(y-yerr))
        errbarpatch(xv, y-yerr, y+yerr, color=colors[i])
        plot(xv, y, linewidth=0.5, color=colors[i])
    end
    xlabel("Time after Response Onset (ms)", fontsize=7, labelpad=3)
    if compare
        lstrings = ["ISI $(round(Int, ISI[i+1]*1000-200)) vs. $(round(Int, ISI[i]*1000-200)) ms" for i = 1:length(ISI)-1]
    else
        lstrings = ["ISI $(round(Int, x)) ms" for x in ISI*1000-200]
    end
    if pos1 && compare
        ylabel("Proportion Variance Explained (\$\\omega^2\$)", fontsize=7, labelpad=1)
    else
        ylabel("Proportion Variance Explained (\$\\omega^2\$)", fontsize=7)
    end
    removespines()
    gca()[:yaxis][:set_ticks_position]("left")
    gca()[:xaxis][:set_ticks_position]("bottom")
    tick_params(axis="x", length=2)
    xticks(fontsize=7)
    yticks(fontsize=7)
    if pos1 && compare
        tight_layout(rect=[-.06, -.06, 1.04, 1.047])
    else
        tight_layout(rect=[-.07, -.06, 1.04, 1.047])
    end
    xlim((-200, 1000))
    limmin = pos1 || compare ? floor(limmin, 4) : 0.0
    plot([0, 0], [limmin, lim*.99], color="k", linestyle=":")
    plot([200, 200], [limmin, lim*.99], color="b", linestyle=":")
    plot([400, 400], [limmin, lim*.99], color="g", linestyle=":")
    plot([700, 700], [limmin, lim*.99], color="r", linestyle=":")
    plot([1000, 1000], [limmin, lim*.99], color="c", linestyle=":", clip_on=false)
    legend(lstrings, fontsize=7, loc="upper right", bbox_to_anchor=(1, 0.9), borderaxespad=0)

    lim = lim+(lim - limmin)*0.09
    ylim(limmin, lim)
    spacing = lim/50

    maxes = maximum(alignedinfoperms, 1)
    thr = mapslices(x->quantile(x, 0.95), maxes, 2)
    for iisi = 1:size(info, 2)
        pv = vec(mean(maxes[1:1, :, iisi] .>= mu[:, iisi], 2))
        sig = pv .< 0.05
        sigs = find(sig)::Vector{Int}
        biggestsig = (-1, -1)
        for i = 1:length(sigs)
            j = i + 1
            while j <= length(sigs) && sigs[j] == sigs[j-1]+1
                sigs[j] - sigs[i] > biggestsig[2] - biggestsig[1] && (biggestsig = (sigs[i], sigs[j]))
                j += 1
            end
        end
        if compare
            println("ISI $(ISI[iisi+1]*1000) - $(ISI[iisi]*1000) min p = $(minimum(pv)) @ $(xv[indmin(pv)])")
            sum(sig) > 0 && println("ISI $(ISI[iisi+1]*1000) - $(ISI[iisi]*1000) significant in $((xv[biggestsig[1]], xv[biggestsig[2]]))")
        else
            println("ISI $(ISI[iisi]*1000) min p = $(minimum(pv)) @ $(xv[indmin(pv)])")
            sum(sig) > 0 && println("ISI $(ISI[iisi]*1000) significant in $((xv[biggestsig[1]], xv[biggestsig[2]]))")
        end
        scatter(xv[find(sig)], fill(lim-(4.5-iisi)*spacing, sum(sig)), 1, edgecolor="none", color=colors[iisi], zorder=10000, marker="s")
    end
    (pos1 || compare) && axhline(0, color="k", zorder=-1)
    savefig("alignedisi$(compare ? "_compare" : "")$(pos1 ? "_pos1" : "").pdf")
end

## Maintenance Period Activity

function wm_units(epoch=:all)
    @load joinpath(PROCESSED_DATA_DIR, "permpval.jld") chs pval
    chs_pres = chs
    pval_pres = pval
    if epoch == :all
        @load joinpath(PROCESSED_DATA_DIR, "permpval_maintenance.jld") pval chs
    elseif epoch == :last
        @load joinpath(PROCESSED_DATA_DIR, "permpval_maintenance_pos.jld") pval_last chs
        pval = pval_last
    elseif epoch == :notlast
        @load joinpath(PROCESSED_DATA_DIR, "permpval_maintenance_pos.jld") pval_notlast chs
        pval = pval_notlast
    end
    mask = unit_mask(chs)
    pval = pval[mask]
    chs = chs[mask]
    pval_pres = pval_pres[indexin(chs, chs_pres)]
    presalpha = 0.001

    for alpha2 in (0.01, 0.05, 0.001)
        selsel2 = pval[pval_pres .< presalpha] .< alpha2
        @printf("%s of stimulus-selective units significant at maintenance at α = %s (p = %.2g)\n",
        to_pct(selsel2), alpha2, pvalue(BinomialTest(selsel2, alpha2), tail=:right))
    end

    for alpha2 in (0.01, 0.05, 0.001)
        selsel3 = pval[pval_pres .>= presalpha] .< alpha2
        @printf("%s of non-stimulus-selective units significant at maintenance at α = %s (p = %.2g)\n",
        to_pct(selsel3), alpha2, pvalue(BinomialTest(selsel3, alpha2), tail=:right))
    end

    alpha = 0.01
    selmask = pval .< alpha
    rgd = Dict(region => selmask[(pval_pres .< presalpha) & [ismatch(Sessions.REGIONS[region], ch.site) for ch in chs]] for region in (:PHG, :EC, :H, :A))
    for region in (:PHG, :EC, :H, :A)
        selmask = pval .< alpha
        selrg = rgd[region]
        @printf("%s of units in %s show significant maintenance period activity at α = %s (p = %.2g)\n",
                to_pct(selrg), region, alpha, pvalue(BinomialTest(selrg, alpha), tail=:right))
    end
end

const FIGS2_UNITS = [
    "399e11sb-2007-7-21_16-34-42 LPHG ch46 clu1 MU",
    "405e23sb-2008-3-2_10-52-20 LPHG ch54 clu1 MU",
    "409e4sb-2008-7-23_16-8-10 LPHG ch1 clu1 SU",
    "409e19sb-2008-7-29_3-9-13 LPHG ch8 clu1 SU",
    "409e32sb-2008-8-1_3-2-50 LPHG ch6 clu1 MU",
    "409e32sb-2008-8-1_3-2-50 LPHG ch8 clu1 MU",
    "398e12sb-2007-5-16_15-55-53 REC ch2 clu2 SU",
    "396e33sb-2006-12-11_15-56-47 LMH ch43 clu2 SU",
    "399e11sb-2007-7-21_16-34-42 LAH ch36 clu2 MU",
    "399e11sb-2007-7-21_16-34-42 RA ch5 clu1 MU",
    "399e11sb-2007-7-21_16-34-42 RA ch7 clu2 SU",
    "399e11sb-2007-7-21_16-34-42 RA ch7 clu5 SU",
    "399e11sb-2007-7-21_16-34-42 RA ch8 clu3 SU",
    "407e13sb-2008-4-11_14-36-46 RA ch1 clu2 MU",
    "407e13sb-2008-4-11_14-36-46 RA ch5 clu2 SU"
]

function plot_stimulus_both(session, channel_number, cluster_number; highlight=nothing, out="pdf", small::Bool=false)
    small ? gcf()[:set_size_inches](2, 2.5) : gcf()[:set_size_inches](4.488, 3.5)
    fs = small ? 5 : 7
    nstim = length(session.stimuli)
    gs = gridspec.GridSpec(1, nstim, left=0.093, right=0.98, bottom=0.09, top=1.006, wspace=0.2, hspace=0.0)
    nthimg = 10
    nthraster = 3
    gss = Array{Any}(nstim)
    for stim = 1:nstim
        gs_item = gs[:__getitem__]((0, stim-1))
        inner_gs = gridspec.GridSpecFromSubplotSpec(nthimg, 1, subplot_spec=gs_item, hspace=0.02)
        subplot(inner_gs[:__getitem__](0))
        showstimulus(stim)
        gs_item = inner_gs[:__getitem__](slice(1, nthimg))
        gss[stim] = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_item, hspace=0.3)
    end
    for (iepoch, epoch) in enumerate((:sample, :maintenance))
        const LIMS = epoch == :maintenance ? (0, 2400) : (-200, 800)
        spikes = session.spikes[channel_number][cluster_number]
        trials = stimulusframe(session.trials)
        ticks = epoch == :maintenance ? (0:1000:2000) : 0:600:1000
        local ax
        rasterbin = epoch == :maintenance ? 200 : 50
        if isa(highlight, Void)
            highlighttrials = Int[]
        else
            stims = session.trials[:stimuli].data::Vector{Vector{Int}}
            highlighttrials = find(Bool[highlight in x for x in stims])
        end

        for stim = 1:nstim
            inner_gs = gridspec.GridSpecFromSubplotSpec(nthraster, 1, subplot_spec=gss[stim][:__getitem__]((iepoch-1, 0)), hspace=0.0)
            axa = subplot(inner_gs[:__getitem__]((slice(0, nthraster-1))), zorder=1000)
            local usetrials::Vector{Int}
            if epoch == :sample
                stims = session.trials[:stimuli].data::Vector{Vector{Int}}
                usetrials = find([stim in x for x in stims])
                onsets = trials[:onset][trials[:stimulus] .== stim].data
            elseif epoch == :maintenance
                stims = session.trials[:stimuli].data::Vector{Vector{Int}}
                usetrials = find([stim in x for x in stims])
                onsets = session.trials[:maint][usetrials]
            end
            t1 = findin(usetrials, setdiff(usetrials, highlighttrials))
            t2 = findin(usetrials, highlighttrials)

            plotraster(spikes, onsets[[t1; t2]], LIMS)
            if !isa(highlight, Void)
                axhline(length(t1)+0.5, color=(1, 0, 0, 0.5), linewidth=0.5)
            end
            axis("off")
            axvline(0, ymin=-1/(nthraster-1), ymax=1, linestyle=":", color="darkcyan", dashes=(1, 1),
                       linewidth=0.5, clip_on=false)

            if stim == 1
                ylabel("Presentation Number")
                bax = ax = subplot(inner_gs[:__getitem__](nthraster-1))
            else
                gca()[:set_yticklabels](())
                bax = subplot(inner_gs[:__getitem__](nthraster-1), sharey=ax)
            end
            rast = raster(spikes, onsets, (LIMS[1]:rasterbin:LIMS[2])/1000)
            bar(LIMS[1]:rasterbin:LIMS[2]-rasterbin, mean(rast, 2)*1000/rasterbin, rasterbin, color="k", edgecolor="none", clip_on=false)
            if !isa(highlight, Void)
                bar(LIMS[1]:rasterbin:LIMS[2]-rasterbin, mean(rast[:, t1], 2)*1000/rasterbin,
                    rasterbin, color=(1, 0, 0, 0.5), edgecolor="none", clip_on=false)
            end
            if epoch == :maintenance
                xticks(ticks, ["0"; [string(Int(x/1000)) for x in ticks[2:end]]], fontsize=fs)
                stim == 1 && xlabel("Time After Mask Offset (Seconds)", fontsize=fs, labelpad=small ? 1 : 2, horizontalalignment="left", x=-0.6)
            else
                xticks(ticks, fontsize=fs)
                stim == 1 && xlabel("Time After Stimulus Onset (ms)", fontsize=fs, labelpad=small ? 1 : 2, horizontalalignment="left", x=-0.6)
            end
            if stim == 1
                ylabel("Firing Rate (Hz)", fontsize=fs, labelpad=small ? 1 : 3)
                removespines(["top", "right"])
                gca()[:yaxis][:set_ticks_position]("left")
                gca()[:xaxis][:set_ticks_position]("bottom")
            else
                bax[:axes][:get_yaxis]()[:set_visible](false)
                removespines(["left", "top", "right"])
                gca()[:xaxis][:set_ticks_position]("bottom")
            end
            tick_params(axis="y", pad=2)
            tick_params(axis="x", direction="out", length=3, pad=1)
            xlim(LIMS[1], LIMS[2])
        end
        sca(ax)
        yl = ax[:get_ylim]()
        yticks((0:max(div(yl[2], 3), 1):yl[2]), fontsize=fs)
    end
    savefig("$(session.folder)_ch$(channel_number)c$(cluster_number)_$(session.sites[channel_number])_$(session.unit_type[channel_number][cluster_number])_samplemaintenance$(small ? "_small" : "").$(out)", dpi=300)
end

function selunits_table()
    patientid(ch) = parse(Int, ch.session[1:3])
    function nperpatient(patients, chs)
        n = zeros(Int, length(patients))
        for ch in chs
            n[searchsortedfirst(patients, patientid(ch))] += 1
        end
        n
    end

    @load joinpath("processed", "permpval_maintenance.jld") chs pval
    allunit = chs[Summary.unit_mask(chs)]
    chpatient = map(patientid, allunit)
    vissel = Analysis.permselunits()
    maintsel = intersect(vissel, chs[pval .< 0.01])
    patients = sort(unique(chpatient))
    nallunit = nperpatient(patients, allunit)
    nvissel = nperpatient(patients, vissel)
    nmaintsel = nperpatient(patients, maintsel)
    push!(nallunit, sum(nallunit))
    push!(nvissel, sum(nvissel))
    push!(nmaintsel, sum(nmaintsel))
    println("Patient\t# Units\t# Vis\t# Maint\t% Vis\t\t% Maint")
    for i = 1:length(patients)+1
        visci = rcopy(R"binom.test($(nvissel[i]), $(nallunit[i]))$conf.int")
        if nvissel[i] == 0
            maintci = [NaN, NaN]
        else
            maintci = rcopy(R"binom.test($(nmaintsel[i]), $(nvissel[i]))$conf.int")
        end
        patient = i > length(patients) ? "TOTAL" : patients[i]
        @printf("%s\t%d\t%d\t%d\t%.0f%% (%.0f%%-%.0f%%)\t%.0f%% (%.0f%%-%.0f%%)\n",
        patient, nallunit[i], nvissel[i], nmaintsel[i],
        nvissel[i]./nallunit[i]*100, visci[1]*100, visci[2]*100,
        nmaintsel[i]./nvissel[i]*100, maintci[1]*100, maintci[2]*100)
    end
end

function meancortest(x, y, nperm)
    cors = map(cor, x, y)
    meancor = tanh(mean(atanh(cors)))
    permcor = zeros(nperm)
    pvalunit = zeros(length(x))
    for ix = 1:length(x)
        permx = copy(x[ix])
        cury = y[ix]
        unitgt = 0
        for iperm = 1:nperm
            shuffle!(permx)
            v = cor(permx, cury)
            permcor[iperm] += atanh(v)
            unitgt += v >= cors[ix]
        end
        pvalunit[ix] = unitgt/nperm
    end
    map!(x->tanh(x/length(cors)), permcor)
    (meancor, cors, mean(permcor .>= meancor) + 1/nperm, pvalunit)
end

function corscatter()
    wd = 6.5
    ht = 4.55
    regions = [:PHG, :EC, :H, :A]
    region_names = ["PHC", "EC", "Hippocampus", "Amygdala"]
    gcf()[:set_size_inches](wd, ht)
    subplots_adjust(wspace=0.35, hspace=0.4, left=0.072, right=0.995, bottom=0.08, top=0.965)
    selselch = intersect(Analysis.maintselunits(), Analysis.permselunits())
    @load joinpath("processed", "modulation_maintenance.jld") means chs
    sel = findin(chs, selselch)
    chsel = chs
    n = map(length, means[sel])
    @load joinpath("processed", "modulation_presentation.jld") stim_fr_win win chs baseline_fr_all baseline_fr_win
    sel2 = findin(chs, selselch)
    sel2region = [findfirst(regions, siteregion(x)) for x in chs[sel2]]
    sp = sortperm(selselch, lt=(x,y)-> siteregion(x) == siteregion(y) ? string(x.site)[1] < string(y.site)[1] :
                                                                        findfirst(regions, siteregion(x)) < findfirst(regions, siteregion(y)))
    sel = sel[sp]
    sel2 = sel2[sp]
    sel2region = sel2region[sp]
    @assert chsel[sel] == chs[sel2] == selselch[sp]
    fig3ch = findfirst(x->x.session == "399e11sb-2007-7-21_16-34-42" && x.channel == 8 && x.cluster == 2, selselch[sp])
    for (i, isel) = enumerate([1:fig3ch-1; fig3ch+1:length(sel); fig3ch])
        subplot(4, 6, i)
        x = stim_fr_win[sel2[isel]]./(win[sel2[isel]][2] - win[sel2[isel]][1])*1000
        y = means[sel[isel]]./(2400 - 300)*1000
        plot(x, y, color="k", marker="o", markerfacecolor="none", markeredgewidth=0.5, linestyle="none", markersize=2)
        rgx = maximum(x) - minimum(x)
        rgy = maximum(y) - minimum(y)
        stepx = rgx > 20 ? 10 : rgx > 10 ? 5 : rgx > 6 ? 3 : rgx > 6 ? 2 : rgx > 2 ? 1 : rgx > 1 ? 0.5 : rgx > 0.6 ? 0.3 : 0.1
        stepy = rgy > 10 ? 5 : rgy > 6 ? 2 : rgy > 3 ? 1 : rgy > 1.5 ? 0.5 : rgy > 0.6 ? 0.2 : rgy > 0.3 ? 0.1 : 0.05
        xticks(-100*stepx:stepx:100*stepx, fontsize=7)
        yticks(-100*stepy:stepy:100*stepy, fontsize=7)
        minx = minimum(x)-rgx/10
        maxx = maximum(x)+rgx/10
        miny = minimum(y)-rgy/10
        maxy = maximum(y)+rgy/10
        xlim(minx, maxx)
        ylim(miny, maxy)
        tick_params(axis="y", direction="out", pad=2, length=2)
        tick_params(axis="x", direction="out", pad=2, length=2)
        gca()[:yaxis][:set_ticks_position]("left")
        gca()[:xaxis][:set_ticks_position]("bottom")
        if (i - 1) % 6 == 0
            ylabel("Maintenance\nFiring Rate (Hz)", fontsize=7, labelpad=2)
        end
        if i > 18
            xlabel("Sample\nFiring Rate (Hz)", fontsize=7, labelpad=1)
        end
        left, bottom, right, top = gca()[:get_position]()[:extents]
        cortxt = "\$r = $(round(cor(x, y), 2))\$"
        idx = findfirst(FIGS2_UNITS, string(chsel[sel[isel]]))
        if idx != 0
            cortxt = "Fig. S2$('A'+idx-1)\n$cortxt"
        elseif isel == fig3ch
            cortxt = "Fig. 3\n$cortxt"
        end
        text(maxx-rgx/30, miny+rgy/50, cortxt, horizontalalignment="right", verticalalignment="baseline", fontsize=8)
        side = string(chs[sel2[isel]].site)[1] == 'L' ? "Left" : "Right"
        title(side*" "*region_names[sel2region[isel]], y=0.97, fontsize=8)
    end
    savefig("corscatter.pdf")
end

function responsestrength()
    function stimhist(stim_fr, baseline_fr)
        resps = [begin
            stimbs = stim .- bs
            bin = maximum(abs(stimbs)) == minimum(abs(stimbs)) ? indmax(stimbs) : indmax(abs(stimbs))
            stimbs[bin]
            end for (stim, bs) in zip(stim_fr, baseline_fr)]
        resps[!isfinite(resps)] = 0
        edges = -4:0.2:9
        respsc = copy(resps)
        respsc[resps .<= edges[1]] = edges[1]+0.01
        respsc[resps .>= edges[end]] = edges[end]-0.01
        h = fit(Histogram, respsc, edges, closed=:right)
        w = h.weights
        (resps, h.edges[1], w)
    end

    gcf()[:set_size_inches](3.5, 2)
    chs_maintenance = Analysis.maintselunits()
    chs_resp = Analysis.permselunits()
    @load joinpath(PROCESSED_DATA_DIR, "modulation_presentation.jld") chs stim_fr_win baseline_fr_win
    stim_fr = stim_fr_win
    baseline_fr = baseline_fr_win
    chsel = indexin(chs_resp, chs)
    vs, edge, w = stimhist(stim_fr[chsel], baseline_fr[chsel])
    bar(edge[1:end-1], w, step(edge), color=(0.5, 0.5, 0.5), linewidth=0.1)
    maintenancesel = findin(chs[chsel], chs_maintenance)
    selmod = vs[maintenancesel]
    nselmod = vs[setdiff(1:length(vs), maintenancesel)]
    selmod = sort(selmod)
    test = MannWhitneyUTest(abs(selmod), abs(nselmod))
    @printf "Median maintenance: %.1f; median no maintenance: %.1f; U = %.3g, p = %.2g\n" median(abs(selmod)) median(abs(nselmod)) test.U pvalue(test)
    showch = chs[chsel][sortperm(vs)]
    ylim(0.0, maximum(w)+1)
    selmod[selmod .>= edge[end]] = edge[end] - step(edge)/2
    selmod = sort!(selmod)
    y = zeros(length(selmod))
    for i = 2:length(selmod)
        i1 = i
        i2 = i-1
        while i2 > 0 && selmod[i1] - selmod[i2] < 0.1
            y[i] += 1
            i2 -= 1
        end
    end
    @assert length(y) == length(selmod)
    scatter(selmod, y/(maximum(w)+1)*25+1.5, edgecolor="k", facecolor="none", linewidth=0.5, clip_on=false, zorder=1000)
    removespines()
    gca()[:yaxis][:set_ticks_position]("left")
    gca()[:xaxis][:set_ticks_position]("bottom")

    xticks(Int(edge[1]):Int(edge[end]))
    xlim(first(edge), last(edge))
    ylabel("Number of Units")
    xlabel("Change in Firing Rate from Baseline (Number of Spikes)")
    tight_layout(pad=0.1)
    savefig("responsestrength.pdf")
end

function wm_correct(epoch=:maint)
    selch = intersect(Analysis.permselunits(), Analysis.maintselunits())
    jld = jldopen(joinpath(PROCESSED_DATA_DIR, "maintenance_topstim_correct.jld"))
    chs = read(jld, "chs")
    pval_mwu = read(jld, "pval_mwu_$epoch")
    sep = read(jld, "sep_$epoch")
    ng1 = read(jld, epoch == :persist ? "ng1_persist" : "ng1")
    ng2 = read(jld, epoch == :persist ? "ng2_persist" : "ng2")
    close(jld)
    alpha = 0.01

    for (title, mych) = [
            ("all visually selective", Analysis.permselunits()),
            ("visually selective, but not maintenance-selective",
            setdiff(Analysis.permselunits(), selch)),
            ("maintenance-selective", selch)
        ]
        println(title)
        sel = findin(chs, mych)
        sel = sel[(ng1[sel] .> 0) & (ng2[sel] .> 0)]
        @printf "  median %s correct trials, %s incorrect trials\n" median(ng1[sel]) median(ng2[sel])
        @printf("  at α = %s, %s significant (p = %.2g)\n", alpha,
                to_pct(pval_mwu[sel] .< alpha),
                pvalue(BinomialTest(pval_mwu[sel] .< alpha, alpha), tail=:right))
        test = OneSampleTTest(sep[sel])
        @printf("  population effect %.2f, t(%d) = %.3f, p = %.4g\n",
                test.xbar, test.df, test.t, pvalue(test))
        @printf "  bootstrap vs. 0 p = %.2g\n" 2*(min(Analysis.bsvs0(sep[sel], 1000000), Analysis.bsvs0(sep[sel], 1000000)))
    end
end
end