module PlotTools
using PyPlot
export plotraster, removespines, errbarpatch, makeline, plotpval

function errbarpatch(x, lower, upper; kwargs...)
    fill([x; reverse(x)], [lower; reverse(upper)]; linewidth=0, alpha=0.3, zorder=0, kwargs...)
end

function plotraster(spikes, onsets, lims)
    for i = 1:length(onsets)
        vlines(scale!(spikes[searchsortedfirst(spikes, onsets[i]+lims[1]/1000):searchsortedlast(spikes, onsets[i]+lims[2]/1000)]-onsets[i], 1000), i - 0.4, i + 0.4, linewidth=0.3)
    end
    xlim(lims[1], lims[2])
    curylim = gca()[:get_ylim]()
    if curylim[2] != 0.5 || curylim[1] < length(onsets)+0.5
        ylim(length(onsets)+0.5, 0.5)
    end
end

function removespines(spines=["top", "right"])
	ax = gca()
	for spine in spines
		ax[:spines][spine][:set_color]("none")
	end
end

function makeline(xstart, xend, ystart, yend, linelength, pv)
    plot([xstart, xstart, xend, xend], [ystart, max(ystart, yend)+linelength, max(ystart, yend)+linelength, yend], color="k")
    plotpval(xstart+(xend-xstart)/2, max(ystart, yend)+linelength, pv)
end

function plotpval(x, y, pv)
    if pv < 0.0001
        pv = "\$p < 10^{-$(floor(Integer, -log10(pv)))}\$"
    elseif pv < 0.001
        pv = "\$p = $(round(pv, 4))\$"
    elseif pv < 0.1
        pv = "\$p = $(round(pv, 3))\$"
    else
        pv = "\$p = $(round(pv, 2))\$"
    end
    text(x, y, verticalalignment="bottom", horizontalalignment="center", pv, fontsize=8)
end
end