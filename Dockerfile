FROM debian:stretch AS rdkit-build-env

RUN apt-get update 
RUN apt-get install -yq --no-install-recommends \
    ca-certificates \
    build-essential \
    cmake \
    wget \
    libboost-dev \
    libboost-system-dev \
    libboost-thread-dev \
    libboost-serialization-dev \
    libboost-python-dev \
    libboost-iostreams-dev \
    libboost-system-dev \
    libboost-regex-dev \
    libcairo2-dev \
    libeigen3-dev \
    python3-dev \
    python3-pip \
    python3-numpy \
    python3-pandas \
    git\
    vim\
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*


COPY modified_rdkit /rdkit


RUN mkdir /rdkit/build
WORKDIR /rdkit/build

#RDK_OPTIMIZE_NATIVE=ON assumes container will be run on the same architecture on which it is built
RUN cmake -Wno-dev \
  -D RDK_INSTALL_INTREE=OFF \
  -D RDK_INSTALL_STATIC_LIBS=OFF \
  -D RDK_BUILD_INCHI_SUPPORT=ON \
  -D RDK_BUILD_AVALON_SUPPORT=ON \
  -D RDK_BUILD_PYTHON_WRAPPERS=ON \
  -D RDK_BUILD_CAIRO_SUPPORT=ON \
  -D RDK_USE_FLEXBISON=OFF \
  -D RDK_BUILD_THREADSAFE_SSS=ON \
  -D RDK_OPTIMIZE_NATIVE=ON \
  -D PYTHON_EXECUTABLE=/usr/bin/python3 \
  -D PYTHON_INCLUDE_DIR=/usr/include/python3.5 \
  -D PYTHON_NUMPY_INCLUDE_PATH=/usr/lib/python3/dist-packages/numpy/core/include \
  -D CMAKE_INSTALL_PREFIX=/usr \
  -D CMAKE_BUILD_TYPE=Release \
  ..

RUN make -j $(nproc) \
 && make -j $(nproc) install



RUN pip3 install setuptools && pip3 install  scipy  && pip3 install scikit-learn Cython mdtraj
WORKDIR /
RUN git clone https://github.com/hjuinj/cpeptools.git
WORKDIR /cpeptools
RUN pip3 install -e .


#COPY scripts /scripts
#WORKDIR /scripts
