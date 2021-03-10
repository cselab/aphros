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
RUN groupadd -r user && useradd -g user -m user
USER user
SHELL ["/bin/bash", "-l", "-c"]
WORKDIR /home/user
ENV MAKEFLAGS='-j4 VERBOSE=1'
ARG GIT_SSL_NO_VERIFY=1
RUN git clone --quiet --single-branch --depth 1 https://github.com/cselab/aphros.git aphros
RUN cd aphros/deploy && ./install_setenv $HOME/.local/aphros
<<<<<<< HEAD
RUN echo '. $HOME/.local/bin/ap.setenv 2>/dev/null' >> $HOME/.profile
RUN cmake -B .work/deploy aphros/deploy
RUN cd .work/deploy && make && make install
RUN cd aphros/src && make
RUN cd aphros/src && make test || true
RUN cd aphros/examples && make build || true
ENTRYPOINT ["/bin/bash", "-l", "-c", "ap.mfer \"$@\"", "ap.mfer"]
=======
RUN echo ". ap.setenv" >> $HOME/.profile
RUN { \
    . $HOME/.profile && \
    cmake -B .work/deploy aphros/deploy && \
    (cd .work/deploy && make && make install) && \
    (cd aphros/src && make && make test) && \
    (cd aphros/examples && make build) ; \
    } 2>&1 | tee .log
CMD exec ap.mfer
>>>>>>> 3a8eb71ad087f6785af460c259cca0aba5199bed
