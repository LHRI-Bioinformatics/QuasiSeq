FROM lhri/django_async
ENV PYTHONUNBUFFERED 1


# RUN apt-get -o Acquire::Check-Valid-Until=false update && apt-get install -y \

RUN apt-get update && apt-get install -y \
    git \
    rsync \
    make \
    cmake \
    curl \
    bzip2 \
    g++ \
    libz-dev \
    gfortran \
    graphviz \
    libjpeg-dev \
    libfreetype6-dev \
    wget \
    default-jre \
    ccache \
    unzip \
    gzip \
    vim \
    bwa

RUN  mkdir /opt/pacbio && git clone https://github.com/PacificBiosciences/pitchfork.git /opt/pacbio/pitchfork
WORKDIR /opt/pacbio/pitchfork
RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.16/src/hdf5-1.8.16.tar.gz -O /opt/pacbio/pitchfork/ports/thirdparty/hdf5/hdf5-1.8.16.tar.gz
RUN git checkout legacy_blasr
RUN make init PREFIX=/opt/pacbio
RUN make blasr PREFIX=/opt/pacbio
ENV PATH="/opt/pacbio/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/pacbio/lib:${LD_LIBRARY_PATH}"

RUN mkdir /opt/MATLAB
WORKDIR /opt/MATLAB
RUN wget https://www.mathworks.com/supportfiles/downloads/R2016b/deployment_files/R2016b/installers/glnxa64/MCR_R2016b_glnxa64_installer.zip
RUN unzip MCR_R2016b_glnxa64_installer.zip
RUN ./install -mode silent -agreeToLicense yes

WORKDIR /opt
RUN wget https://downloads.sourceforge.net/project/bbmap/BBMap_37.24.tar.gz
RUN tar -zxvf BBMap_37.24.tar.gz

WORKDIR /opt
RUN git clone https://github.com/yechengxi/Sparc.git
WORKDIR /opt/Sparc
RUN g++ -O3 -o Sparc *.cpp
RUN chmod +x /opt/Sparc/compiled/Sparc

ENV ROOT_DIR /data
ENV APP_DIR cluster_interface
#RUN mkdir ${ROOT_DIR};
RUN mkdir ${ROOT_DIR}/quasiSeqOut;
WORKDIR ${ROOT_DIR}

#RUN mkdir /config;
ADD requirements.pip /config/requirements.pip
RUN pip install -r /config/requirements.pip

RUN echo deb http://archive.debian.org/debian jessie-backports main >> /etc/apt/sources.list
RUN apt-get -o Acquire::Check-Valid-Until=false update && apt install -t jessie-backports  openjdk-8-jre-headless ca-certificates-java
#RUN update-alternatives --config java
ENV PATH="/usr/lib/jvm/java-8-openjdk-amd64/jre/bin:${PATH}"

RUN pip install gevent

ENV PATH="/opt/Sparc/compiled:${PATH}"

ENV ROOT_DIR /opt/quasi-seq/
ENV APP_DIR cluster_interface
ENV DATA_DIR /data
WORKDIR ${ROOT_DIR}

ENTRYPOINT ["/bin/bash","--init-file","/opt/quasi-seq/init.sh"]