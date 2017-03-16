This repository contains the raw data and analysis files for:

Kornblith, S., Quiroga, R. Q., Koch, C., Fried, I., and Mormann, F. (2017). [Persistent single neuron activity during working memory in the human medial temporal lobe](http://www.sciencedirect.com/science/article/pii/S0960982217301495). _Current Biology_, in press.

## Raw data

This repository contains the spike times and task information for all sessions analyzed in the paper in the [`raw/`](raw/) directory. Each file contains raw data for an individual recording session. The data are represented in [HDF5 format](https://support.hdfgroup.org/HDF5/), with the following structure:

- `channels/`
    - `ch<N>_spike_times`: Spike times in seconds for each unit on channel `N`, encoded as an array of Float64 VLENs (a "ragged" array).
    - `ch<N>_unit_types`: Vector of strings reflecting the type of each unit on channel `N`. "SU" = single unit; "MU" = multi unit; "A" = artifact.
- `sites`: Vector of site label strings for each channel. See [here](src/Sessions.jl#L11-L28) for regular expressions to select specific MTL subregions.
- `trials/`
    - `correct`: Vector of length `ntrials`. 1 if response on trial was correct, or 0 if incorrect.
    - `isi`: Vector of length `ntrials`. Length of inter-stimulus interval in seconds.
    - `probe1`: Vector of length `ntrials`. Index of left probe stimulus.
    - `probe2`: Vector of length `ntrials`. Index of right probe stimulus.
    - `stimuli`: `nstimuli_per_trial` × `ntrials` matrix. Indices of stimuli presented on each trial.
    - `onsets`: `nstimuli_per_trial` × `ntrials` matrix. Onset times in seconds for each stimulus presentation.
    - `offsets`: `nstimuli_per_trial` × `ntrials` matrix. Offset times in seconds for each stimulus presentation.
    - `maint`: Vector of length `ntrials`. Time of start of maintenance period in seconds.
    - `probes`: Vector of length `ntrials`. Time of start of probe presentation in seconds.

## Analysis code

### Requirements

The [Docker image](#docker-image) includes all necessary dependencies. If you do not use the Docker image, you will need to install:

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
    - [IJulia](https://github.com/JuliaLang/IJulia.jl) (for the notebook; tested with v1.4.1)

For the full list of Julia packages used, along with their dependencies, versions, and versions of dependencies, see [this file](REQUIRE).

- [R](https://www.r-project.org/) (tested with v3.2.3) and the following R packages:
    - [exact2x2](https://cran.r-project.org/web/packages/exact2x2/index.html) (tested with v1.4.1)
    - [ez](https://cran.r-project.org/web/packages/ez/index.html) (tested with v4.3)
    - [pwr](https://cran.r-project.org/web/packages/pwr/index.html) (tested with v1.2)

### Jupyter notebook

This repository contains a [Jupyter notebook](Paper.ipynb) that reproduces the plots and values provided in the paper.

### Rerunning analyses

The [`processed/`](processed/) directory contains all intermediate analysis results. If you wish, you can also re-run these analyses using the Julia script [src/analyze_all.jl](src/analyze_all.jl). To parallelize across multiple cores, run `julia -p n src/analyze_all.jl` where `n` is the number of cores to use.

The analyses may take many hours to run, even on a fast machine. Since many analyses are based on permutation tests, exact numbers may differ slightly from those in the provided data files and the paper, because they will run with a different random seed. However, resulting p-values and confidence intervals should be very similar.

### Docker image

A [Docker](https://www.docker.com/what-docker) image containing the code, data, and all dependencies including Julia and Jupyter notebook is [available]((https://cloud.docker.com/swarm/simonster/repository/docker/simonster/kornblith-et-al-current-biology-2017/general)). To use it, [install Docker](https://www.docker.com/community-edition#/download) and run the following command from the shell:

```sh
docker run -it -p=127.0.0.1:8888:8888 simonster/kornblith-et-al-2017-current-biology
```

This will fetch the Docker image and launch an instance of Jupyter notebook. You can run the notebook by navigating to [http://127.0.0.1:8888/notebooks/Paper.ipynb](http://127.0.0.1:8888/notebooks/Paper.ipynb). When asked for a token, enter the 48 character string following "http://.../?token=" that is shown in the shell.

Within the Docker image, the contents of this repository are located in `/root/.julia/v0.5/Kornblith-et-al-2017-Current-Biology`.

## License

### Data

The raw and processed data are made available under the [Open Database License](http://opendatacommons.org/licenses/odbl/1.0/). We consider intermediate results of analyses of the raw data (such as the processed data included here) to constitute derivative works. Thus, works making use of the raw data must make any such intermediate results available under the ODbL. Additionally, we kindly request that works making use of our data make the code necessary to reproduce their analyses publicly available.

### Code

The analysis code is made available under the MIT license:

> Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
