FROM library/julia:0.5

# Install system packages
RUN echo "deb http://cran.us.r-project.org/bin/linux/debian jessie-cran3/" >> /etc/apt/sources.list
RUN apt-get -y update && apt-get install -y --force-yes build-essential libzmq3 libhdf5-8 r-base

# Install Julia dependencies and Jupyter
ADD . /root/Kornblith-et-al-2017-Current-Biology
RUN julia -e 'Pkg.clone("/root/Kornblith-et-al-2017-Current-Biology", "Kornblith-et-al-2017-Current-Biology"); Pkg.add("IJulia")'
RUN rm -rf /root/Kornblith-et-al-2017-Current-Biology
RUN julia -e 'push!(LOAD_PATH, "/root/.julia/v0.5/Kornblith-et-al-2017-Current-Biology/src"); using Summary, IJulia'

# Install R dependencies
RUN R -e "install.packages(c('exact2x2', 'ez', 'pwr'), repos='http://cran.us.r-project.org')"

# Run Jupyter on Docker-external interface
RUN mkdir ~/.jupyter && echo "c.NotebookApp.ip = '172.17.0.2'" > ~/.jupyter/jupyter_notebook_config.py
CMD cd ~/.julia/v0.5/Kornblith-et-al-2017-Current-Biology && ~/.julia/v0.5/Conda/deps/usr/bin/jupyter notebook
