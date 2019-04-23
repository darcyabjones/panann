ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

RUN  apt-get update \
  && apt-get install -y \
       autoconf \
       build-essential \
       git \
       libbz2-dev \
       libcurl4-openssl-dev \
       libgsl-dev \
       liblzma-dev \
       libncurses5-dev \
       libssl-dev \
       libperl-dev \
       zlib1g-dev \
       wget \
  && rm -rf /var/lib/apt/lists/*


ENV HTSLIB_PREFIX="/opt/htslib"
ARG HTSLIB_VERSION="1.9"
ARG HTSLIB_REPO="https://github.com/samtools/htslib.git"

ENV BCFTOOLS_PREFIX="/opt/bcftools"
ARG BCFTOOLS_VERSION="1.9"
ARG BCFTOOLS_REPO="https://github.com/samtools/bcftools.git"

ENV SAMTOOLS_PREFIX="/opt/samtools"
ARG SAMTOOLS_VERSION="1.9"
ARG SAMTOOLS_REPO="https://github.com/samtools/samtools.git"


WORKDIR /tmp/htslib
RUN  git clone ${HTSLIB_REPO} . \
  && git fetch --tags \
  && git checkout ${HTSLIB_VERSION} \
  && autoheader \
  && autoconf \
  && ./configure --prefix="${HTSLIB_PREFIX}" --enable-libcurl \
  && make -j $(grep -c ^processor /proc/cpuinfo) \
  && make -j $(grep -c ^processor /proc/cpuinfo) install


WORKDIR /tmp/bcftools
RUN  git clone ${BCFTOOLS_REPO} . \
  && git fetch --tags \
  && git checkout ${BCFTOOLS_VERSION} \
  && autoheader \
  && autoconf \
  && ./configure \
       --prefix="${BCFTOOLS_PREFIX}" \
       --enable-perl-filters \
       --enable-lib-gsl \
       --with-htslib="${HTSLIB_PREFIX}" \
  && make -j $(grep -c ^processor /proc/cpuinfo) \
  && make -j $(grep -c ^processor /proc/cpuinfo) install


WORKDIR /tmp/samtools
RUN  git clone ${SAMTOOLS_REPO} . \
  && git fetch --tags \
  && git checkout ${SAMTOOLS_VERSION} \
  && autoheader \
  && autoconf \
  && ./configure \
       --prefix="${SAMTOOLS_PREFIX}" \
       --with-htslib="${HTSLIB_PREFIX}" \
  && make -j $(grep -c ^processor /proc/cpuinfo) \
  && make -j $(grep -c ^processor /proc/cpuinfo) install


FROM debian:${DEBIAN_VERSION}

# Install runtimes
# Suprisingly libcurl3 contains runtime libs for libcurl4-dev.
# libssl1.1 ?? only for aws support
# libperl5.24 ?? for bcftools with perl filters
# libgsl2 ?? for bcftools polysomy command
# libncurses5 ?? for samtools tview command.
RUN  apt-get update \
  && apt-get install -y \
       libbz2-1.0 \
       libcurl3 \
       libgsl2 \
       libncurses5 \
       libperl5.24 \
       libssl1.1 \
       lzma \
       zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV HTSLIB_PREFIX="/opt/htslib"
ENV BCFTOOLS_PREFIX="/opt/bcftools"
ENV SAMTOOLS_PREFIX="/opt/samtools"

COPY --from=builder "${HTSLIB_PREFIX}" "${HTSLIB_PREFIX}"
COPY --from=builder "${BCFTOOLS_PREFIX}" "${BCFTOOLS_PREFIX}"
COPY --from=builder "${SAMTOOLS_PREFIX}" "${SAMTOOLS_PREFIX}"

ENV PATH="${PATH}:${SAMTOOLS_PREFIX}/bin:${BCFTOOLS_PREFIX}/bin:${HTSLIB_PREFIX}/bin"
ENV CPATH="${HTSLIB_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LD_LIBRARY_PATH}"
