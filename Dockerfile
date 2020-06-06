FROM ubuntu:20.04
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
RUN groupadd -r user && useradd -g user -m user
USER user
ENTRYPOINT ["/bin/sh", "-l"]
WORKDIR /home/user
ARG MAKEFLAGS='-j 1 VERBOSE=1'
ENV GIT_SSL_NO_VERIFY=1
RUN git clone --quiet --single-branch --depth 1 https://github.com/cselab/aphros.git aphros
RUN cd aphros/deploy && ./install_setenv $HOME/.local/aphros
RUN echo ". ap.setenv" >> $HOME/.profile
RUN { \
    . $HOME/.profile && \
    cmake -B .work/deploy aphros/deploy && \
    (cd .work/deploy && make && make install) && \
    (cd aphros/src && make && make test) && \
    (cd aphros/examples && make build) ; \
    } 2>&1 | tee .log
