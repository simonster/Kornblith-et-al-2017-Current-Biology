module Sessions
using DataFrames, NaturalSort, DataStructures, HDF5

export PROCESSED_DATA_DIR, REGIONS, REGION_LABELS, Session, ClusterMetadata, sessionnames, readsession, patient, siteregion

const RAW_DATA_DIR = joinpath(dirname(@__FILE__), "..", "raw")
const PROCESSED_DATA_DIR = joinpath(dirname(@__FILE__), "..", "processed")
const REGIONS = OrderedDict(
    :H => r"^..H$",
    :EC => r"^.EC",
    :A => r"^.A$",
    :PHG => r"^.(?:PHG|PHC|PG.?)$",
    :MTL => r"^.(?:.H|EC|A|PHG|PHC|PG.?)$",
    :All => r""
)
const REGION_LABELS = OrderedDict(
    :H => "Hippocampus",
    :EC => "Entorhinal Cortex",
    :A => "Amygdala",
    :PHG => "Parahippocampal Cortex",
    :MTL => "Medial Temporal Lobe",
    :All => "All"
)

type Session
    folder::String
    trials::DataFrame
    spikes::Vector{Vector{Vector{Float64}}}
    unit_type::Vector{Vector{Symbol}}
    sites::Vector{String}
    stimuli::Vector{String}
end

sessionnames() = sort!([x[1:rsearch(x, '.')-1] for x in readdir(RAW_DATA_DIR) if endswith(x, ".h5")], lt=natural)
splitcols(x) = [x[:, i] for i = 1:size(x, 2)]
function readsession(sname)::Session
    h5open(joinpath(RAW_DATA_DIR, "$sname.h5")) do f
        ftrials = f["trials"]
        stim = read(ftrials, "stimuli")
        nstim = maximum(stim)
        trials = DataFrame(correct=read(ftrials, "correct") .== 1,
                          isi=read(ftrials, "isi"),
                          probe1=read(ftrials, "probe1"),
                          probe2=read(ftrials, "probe2"),
                          stimuli=splitcols(stim),
                          onsets=splitcols(read(ftrials, "onsets")),
                          offsets=splitcols(read(ftrials, "offsets")),
                          maint=read(ftrials, "maint"),
                          probes=read(ftrials, "probes"))
        close(ftrials)
        fchannels = f["channels"]

        sites = read(f, "sites")
        spikes = Vector{Vector{Vector{Float64}}}(length(sites))
        unit_type = Vector{Vector{Symbol}}(length(sites))
        for i = 1:length(sites)
            st = "ch$(i)_spike_times"
            if exists(fchannels, st)
                spikes[i] = read(fchannels, st)
                unit_type[i] = map(Symbol, read(fchannels, "ch$(i)_unit_types"))
            end
        end

        # Not included in data files for reasons of confidentiality
        stimuli = ["Stimulus $i" for i in 1:nstim]

        Session(sname, trials, spikes, unit_type, sites, stimuli)
    end
end

immutable ClusterMetadata
    session::String
    site::String
    channel::Int
    cluster::Int
    unit_type::Symbol
end

ClusterMetadata(s::Session, ch::Int, cluster::Int) =
    ClusterMetadata(s.folder, s.sites[ch], ch, cluster, s.unit_type[ch][cluster])

Base.:(==)(x::ClusterMetadata, y::ClusterMetadata) =
    x.cluster == y.cluster && x.channel == y.channel && x.session == y.session
Base.hash(x::ClusterMetadata) = hash(x.session, hash(x.site, hash(x.channel, hash(x.cluster))))
Base.isless(x::ClusterMetadata, y::ClusterMetadata) =
    x.session != y.session ? natural(x.session, y.session) : x.channel != y.channel ? isless(x.channel, y.channel) :
    isless(x.cluster, y.cluster)
Base.show(io::IO, m::ClusterMetadata) =
    print(io, m.session, ' ', m.site, " ch", m.channel, " clu", m.cluster, ' ', m.unit_type)

function siteregion(m::ClusterMetadata)
    for (k, v) in REGIONS
        (k == :All || k == :MTL) && continue
        ismatch(v, m.site) && return k
    end
    ismatch(REGIONS[:MTL], m.site) && return :MTL
    :All
end
end
