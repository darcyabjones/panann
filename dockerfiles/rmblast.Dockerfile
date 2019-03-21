ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       cpio \
       libblas-dev \
       libboost-dev \
       libboost-test-dev \
       libboost-thread-dev \
       libboost-regex-dev \
       libboost-iostreams-dev \
       libboost-program-options-dev \
       libboost-system-dev \
       libboost-test-dev \
       libboost-filesystem-dev \
       libbz2-dev \
       libcurl4-gnutls-dev \
       libgmp-dev \
       libgnutls28-dev \
       liblapack-dev \
       liblzma-dev \
       libperl-dev \
       libxml2-dev \
       libxslt1-dev \
       zlib1g-dev \
       wget \
  && rm -rf /var/lib/apt/lists/*


ENV RMBLAST_PREFIX="/opt/rmblast"

WORKDIR "/tmp/blast"
RUN  wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-src.tar.gz \
  && wget http://www.repeatmasker.org/isb-2.6.0+-changes-vers2.patch.gz \
  && tar zxf ncbi-blast-2.6.0+-src.tar.gz \
  && rm ncbi-blast-2.6.0+-src.tar.gz \
  && gunzip isb-2.6.0+-changes-vers2.patch.gz \
  && cd ncbi-blast-2.6.0+-src \
  && patch -p1 < ../isb-2.6.0+-changes-vers2.patch \
  && cd "c++" \
  && ./configure \
       --with-mt \
       --prefix="${RMBLAST_PREFIX}" \
       -j $(grep -c ^processor /proc/cpuinfo) \
       --without-debug \
  && make \
  && make install


FROM debian:${DEBIAN_VERSION}

ENV RMBLAST_PREFIX="/opt/rmblast"
COPY --from=builder "${RMBLAST_PREFIX}" "${RMBLAST_PREFIX}"

RUN  apt-get update \
  && apt-get install -y \
       libblas3 \
       libbz2-1.0 \
       libcurl3-gnutls \
       libgmp10 \
       libgnutlsxx28 \
       liblapack3 \
       libstdc++6 \
       libxml2 \
       libxslt1.1 \
       lzma \
       perl \
       zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV PATH="${RMBLAST_PREFIX}/bin:${PATH}"
