module TaskData
using Sessions, DataFrames, JLD
export stimulusframe, sessinfo, ISI

const ISI = [0.0, 0.2, 0.5, 0.8]

function stimulusframe(trials::DataFrame)
    trial_stimuli = trials[:stimuli].data::Vector{Vector{Int}}
    trial_onsets = trials[:onsets].data::Vector{Vector{Float64}}
    trial_offsets = trials[:offsets].data::Vector{Vector{Float64}}
    trial_correct::Vector{Bool} = trials[:correct].data
    trial_isi = trials[:isi].data::Vector{Float64}
    trial_probe1 = trials[:probe1].data::Vector{Int}
    trial_probe2 = trials[:probe2].data::Vector{Int}

    stimulus = zeros(Int, sum(length, trial_stimuli))
    next_stimulus = zeros(Int, sum(length, trial_stimuli))
    prev_stimulus = zeros(Int, sum(length, trial_stimuli))
    trial_number = zeros(Int, length(stimulus))
    pos = zeros(Int, length(stimulus))
    onset = fill(NaN, length(stimulus))
    offset = fill(NaN, length(stimulus))
    isi = fill(NaN, length(stimulus))
    probed = fill(false, length(stimulus))
    correct = DataArray(Bool, length(stimulus))
    n = 0
    for i = 1:length(trial_stimuli)
        rg = n+1:n+length(trial_stimuli[i])
        trial_number[rg] = i
        pos[rg] = 1:length(rg)
        stimulus[rg] = trial_stimuli[i]
        prev_stimulus[rg[2:end]] = trial_stimuli[i][1:end-1]
        next_stimulus[rg[1:end-1]] = trial_stimuli[i][2:end]
        onset[rg] = trial_onsets[i]
        offset[rg] = trial_offsets[i]
        isi[rg] = trial_isi[i]
        pb = probed[rg] = (trial_stimuli[i] .== trial_probe1[i]) | (trial_stimuli[i] .== trial_probe2[i])
        @assert sum(pb) == 1
        correct[rg[findfirst(pb)]] = trial_correct[i]
        n = rg[end]
    end
    @assert n == length(stimulus)

    DataFrame(trial_number=trial_number, pos=pool(pos), stimulus=pool(stimulus),
              prev_stimulus=pool(prev_stimulus), next_stimulus=pool(next_stimulus),
              onset=onset, offset=offset, isi=pool(isi), probed=probed, correct=correct)
end

function sessinfo()
    patients = Int[]
    correct = Float64[]
    trials = Int[]
    stimuli = Int[]
    npertrial = Int[]
    maintenancetime = Tuple{Float64,Float64,Float64}[]
    snames = sessionnames()
    for sname in snames
        s = readsession(sname)
        push!(patients, parse(Int, sname[1:3]))
        push!(correct, mean(s.trials[:correct]))
        push!(trials, size(s.trials, 1))
        push!(stimuli, length(s.stimuli))
        push!(npertrial, length(s.trials[:stimuli][1]))
        maintenancetimes = s.trials[:probes] .- s.trials[:maint]
        push!(maintenancetime, (median(maintenancetimes), extrema(maintenancetimes)...))
    end
    DataFrame(session=snames, correct=correct, trials=trials, stimuli=stimuli,
              npertrial=npertrial, maintenancetime=maintenancetime)
end
end