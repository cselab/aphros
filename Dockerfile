FROM ubuntu:20.04
ARG MAKEFLAGS='-j 2'
RUN DEBIAN_FRONTEND=noninteractive apt update
RUN DEBIAN_FRONTEND=noninteractive apt install --yes --quiet --no-install-recommends \
cmake \
g++ \
git \
libgsl-dev \
hdf5-tools \
libhdf5-openmpi-dev \
make \
man-db \
python3 \
python3-numpy \
rsync
RUN echo root:g | chpasswd
RUN useradd -ms /bin/sh aphros
USER aphros
WORKDIR /home/aphros
RUN GIT_SSL_NO_VERIFY=1 \
    git clone --quiet --single-branch --depth 1 https://github.com/cselab/aphros.git aphros
RUN cd aphros/deploy && ./install_setenv $HOME/.local/aphros
RUN { \
    . $HOME/.profile && \
    . ap.setenv && \
    cmake -B .work/deploy aphros/deploy && \
    (cd .work/deploy && make && make install) && \
    (cd aphros/src && make && make test) && \
    (cd aphros/examples && make build) ; \
    } 2>&1 | tee .log
