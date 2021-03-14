FROM ubuntu:20.04
RUN DEBIAN_FRONTEND=noninteractive apt-get -qq update
RUN DEBIAN_FRONTEND=noninteractive apt-get -qq install --yes --no-install-recommends \
cmake \
dumb-init \
g++ \
git \
hdf5-tools \
libhdf5-mpich-dev \
make \
man-db \
mpich \
python3 \
python3-numpy \
rsync
RUN echo root:g | chpasswd
SHELL ["/bin/bash", "-l", "-c"]
ARG MAKEFLAGS=-j4
ENV GIT_SSL_NO_VERIFY=1
RUN git clone --quiet --single-branch --depth 1 https://github.com/cselab/aphros.git aphros
RUN cd aphros/deploy && ./install_setenv $HOME/.local/aphros
RUN echo '. $HOME/.local/bin/ap.setenv 2>/dev/null' >> $HOME/.profile
RUN mkdir -p .work/deploy
RUN cd .work/deploy && cmake /aphros/deploy
RUN cd .work/deploy && make
RUN cd .work/deploy && make install
RUN cd aphros/src && make
RUN cd aphros/src && dumb-init make test || true
RUN cd aphros/examples && make build || true
ENTRYPOINT ["dumb-init", "/bin/bash", "-l", "-c", "ap.mfer \"$@\"", "ap.mfer"]
