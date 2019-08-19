FROM jupyter/minimal-notebook

# Octave
USER root

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python-sympy \
        octave \
        octave-symbolic \
        octave-miscellaneous \
        octave-io \
        octave-control \
        gnuplot \
        ghostscript && \
        apt-get -qq clean && rm -rf /var/lib/apt/lists/*

USER $NB_UID

RUN conda install --quiet --yes \
    'octave_kernel' && \
    conda clean --all -f -y && \
    fix-permissions $CONDA_DIR && \
    
