module SpikeTools
export raster, paddedbins, boxcarsmooth

function raster(spike_times::Vector{Float64}, start_times, r)
    firstr = first(r)
    lastr = last(r)
    n = length(r) - 1
    rasters = zeros(Int32, length(r) - 1, length(start_times))

    n = length(spike_times)
    j = 1
    for i = 1:length(start_times)
        start_time = start_times[i]
        a = start_time + firstr
        b = start_time + lastr
        j = searchsortedfirst(spike_times, a, Base.Order.Forward)
        if j <= n
            j += spike_times[j] == a
            while j <= n && spike_times[j] < b
                rasters[searchsortedlast(r, spike_times[j] - start_time, Base.Order.Forward), i] += 1
                j += 1
            end
        end
    end

    rasters
end

function paddedbins(lims::Tuple{Real, Real}, step::Real, expand_by::Real)
    expand_by = div(expand_by, step)
    extract_times = (lims[1]-expand_by*step:step:lims[2]+expand_by*step)/1000
    extract_range = expand_by+1:length(extract_times)-expand_by-1
    times = extract_times[extract_range]*1000+step/2
    (extract_times, extract_range, times)
end

function boxcarsmooth{T<:Real}(rast, fwhm, outtype::Type{T}=Float64)
    convrast = zeros(T, size(rast))
    for j = 1:size(rast, 2)
        halffwhm = fwhm÷2
        v = 0
        for i = 1:2halffwhm
            v += rast[i, j]
        end
        convrast[halffwhm+1, j] = v
        for i = halffwhm+2:size(rast, 1)-halffwhm
            v -= rast[i-halffwhm-1, j]
            v += rast[i+halffwhm-1, j]
            convrast[i, j] = v
        end
    end
    convrast
end
end