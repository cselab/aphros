FROM ubuntu:20.04
RUN DEBIAN_FRONTEND=noninteractive apt update
RUN DEBIAN_FRONTEND=noninteractive apt install --yes --quiet --no-install-recommends \
cmake \
g++ \
git \
hdf5-tools \
libhdf5-openmpi-dev \
make \
man-db \
python3 \
python3-numpy \
rsync
RUN echo root:g | chpasswd
SHELL ["/bin/bash", "-l", "-c"]
ENV MAKEFLAGS='-j4 VERBOSE=1'
ENV GIT_SSL_NO_VERIFY=1
RUN git clone --quiet --single-branch --depth 1 https://github.com/cselab/aphros.git aphros
RUN cd aphros/deploy && ./install_setenv $HOME/.local/aphros
RUN echo '. $HOME/.local/bin/ap.setenv 2>/dev/null' >> $HOME/.profile
RUN cmake -B .work/deploy aphros/deploy
RUN cd .work/deploy && make && make install
RUN cd aphros/src && make
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
RUN cd aphros/src && make test || true
RUN cd aphros/examples && make build || true
ENTRYPOINT ["/bin/bash", "-l", "-c", "ap.mfer \"$@\"", "ap.mfer"]
