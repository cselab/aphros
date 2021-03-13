# hash:sha256:07458f591ce56f411b034c9239278836dc72efef3551e9c146078da78c220107
FROM registry.codeocean.com/codeocean/ubuntu:18.04.5
ARG DEBIAN_FRONTEND=noninteractive
RUN apt update
RUN apt install --yes --quiet --no-install-recommends \
g++ \
git \
hdf5-tools \
libhdf5-mpich-dev \
mpich \
make \
man-db \
python3 \
python3-numpy \
cmake \
rsync
RUN echo root:g | chpasswd
SHELL ["/bin/bash", "-l", "-c"]
ENV GIT_SSL_NO_VERIFY=1