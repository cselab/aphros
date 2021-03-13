# hash:sha256:07458f591ce56f411b034c9239278836dc72efef3551e9c146078da78c220107
FROM registry.codeocean.com/codeocean/ubuntu:18.04.5
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qq update 
RUN apt-get -qq install --yes --no-install-recommends \
cmake \
dumb-init \
g++ \
git \
hdf5-tools \
libhdf5-openmpi-dev \
make \
man-db \
openmpi-bin \
paraview \
paraview-python \
python3 \
python3-numpy \
rsync \
ssh \
xvfb
RUN echo root:g | chpasswd
SHELL ["/bin/bash", "-l", "-c"]
ENV GIT_SSL_NO_VERIFY=1