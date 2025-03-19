FROM continuumio/miniconda3

#TODO: add a user, don't run as root

# Copy the environment file
COPY env.yml /tmp/env.yml

RUN conda env create -f /tmp/env.yml && conda clean --all -y && \
    echo "source activate variant_calling_env" > ~/.bashrc
    
ENV PATH /opt/conda/envs/variant_calling_env/bin:$PATH



