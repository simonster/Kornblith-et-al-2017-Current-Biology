This repository contains the raw data and analysis files for:

Kornblith, S., Quiroga, R. Q., Koch, C., Fried, I., and Mormann, F. (2017). Persistent single neuron activity during working memory in the human medial temporal lobe. _Current Biology_, in press.

## Raw data

This repository contains the spike times and task information for all sessions analyzed in the paper in the [raw/](raw/) directory. The data is represented in [HDF5 format](https://support.hdfgroup.org/HDF5/), with the following structure:

- `channels/`
    - `ch<N>_spike_times`: Spike times in seconds for each unit on channel `N`, encoded as an array of Float64 VLENs (a "ragged" array).
    - `ch<N>_unit_types`: Vector of strings reflecting the type of each unit on channel `N`. "SU" = single unit; "MU" = multi unit; "A" = artifact.
- `sites`: Vector of site label strings for each channel. See [here](src/Sessions.jl#L11-L28) for regular expressions to select specific MTL subregions.
- `trials/`
    - `correct`: Vector of length `ntrials`. 1 if response on trial was correct, or 0 if incorrect.
    - `isi`: Vector of length `ntrials`. Inter-stimulus interval in seconds.
    - `probe1`: Vector of length `ntrials`. Index of left probe stimulus.
    - `probe2`: Vector of length `ntrials`. Index of right probe stimulus.
    - `stimuli`: `nstimuli_per_trial` × `ntrials` matrix. Indices of stimuli presented on each trial.
    - `onsets`: `nstimuli_per_trial` × `ntrials` matrix. Onset times in seconds for each stimulus presentation.
    - `offsets`: `nstimuli_per_trial` × `ntrials` matrix. Offset times in seconds for each stimulus presentation.
    - `maint`: Vector of length `ntrials`. Start of maintenance period in seconds.
    - `probes`: Vector of length `ntrials`. Start of probe presentation in seconds.

## Analysis code

### Requirements

The included code requires:

- [Julia](https://github.com/JuliaLang/julia) v0.5 and the following Julia packages:
    - [DataFrames.jl](https://github.com/JuliaStats/DataFrames.jl) (tested with v0.8.5)
    - [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) (tested with v0.12.0)
    - [DataStructures.jl](https://github.com/JuliaCollections/DataStructures.jl) (tested with v0.5.2)
    - [HypothesisTests.jl](https://github.com/JuliaStats/HypothesisTests.jl) (tested with v0.4.0)
    - [JLD.jl](https://github.com/JuliaIO/JLD.jl) (tested with v0.6.9)
    - [PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl) (tested with v2.2.4)
    - [PyCall.jl](https://github.com/stevengj/PyCall.jl) (tested with v1.9.0)
    - [StatsBase](https://github.com/JuliaStats/StatsBase.jl) (tested with v0.13.0)
    - [RCall](https://github.com/JuliaInterop/RCall.jl) (tested with v0.6.4)

For the full list of Julia packages used, along with their dependencies, versions, and versions of dependencies, see [this file](REQUIRE).

- [R](https://www.r-project.org/) (tested with v3.2.3) and the following R packages:
    - [exact2x2](https://cran.r-project.org/web/packages/exact2x2/index.html) (tested with v1.4.1)
    - [ez](https://cran.r-project.org/web/packages/ez/index.html) (tested with v4.3)
    - [pwr](https://cran.r-project.org/web/packages/pwr/index.html) (tested with v1.2)

### Jupyter notebook

This repository contains a [Jupyter notebook](Paper.ipynb) that reproduces the plots and values provided in the paper.

### Rerunning analyses

To rerun all analyses, run the Julia script [src/analyze_all.jl](src/analyze_all.jl). To parallelize across multiple cores, run `julia -p n src/analyze_all.jl` where `n` is the number of cores to use.

The analyses may take many hours to run, even on a fast machine. Since many analyses are based on permutation tests, exact numbers may differ slightly from those in the paper, because they will run with a different random seed. However, resulting p-values should be very similar to those presented.

### Docker image

We plan to provide a Docker image containing the data, code, and all dependencies shortly.
