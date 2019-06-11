ARG IMAGE="darcyabjones/base"

FROM ${IMAGE} as builder

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

ARG HTSLIB_TAG
ARG BCFTOOLS_TAG
ARG SAMTOOLS_TAG

ARG HTSLIB_REPO="https://github.com/samtools/htslib.git"
ARG BCFTOOLS_REPO="https://github.com/samtools/bcftools.git"
ARG SAMTOOLS_REPO="https://github.com/samtools/samtools.git"

ARG HTSLIB_PREFIX_ARG="/opt/htslib/${HTSLIB_TAG}"
ARG BCFTOOLS_PREFIX_ARG="/opt/bcftools/${BCFTOOLS_TAG}"
ARG SAMTOOLS_PREFIX_ARG="/opt/samtools/${SAMTOOLS_TAG}"

ENV HTSLIB_PREFIX="${HTSLIB_PREFIX_ARG}"
ENV BCFTOOLS_PREFIX="${BCFTOOLS_PREFIX_ARG}"
ENV SAMTOOLS_PREFIX="${SAMTOOLS_PREFIX_ARG}"

WORKDIR /tmp/htslib
RUN  git clone "${HTSLIB_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${HTSLIB_TAG}" \
  && autoheader \
  && autoconf \
  && ./configure --prefix="${HTSLIB_PREFIX}" --enable-libcurl \
  && make -j $(grep -c ^processor /proc/cpuinfo) \
  && make -j $(grep -c ^processor /proc/cpuinfo) install


WORKDIR /tmp/bcftools
RUN  git clone "${BCFTOOLS_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${BCFTOOLS_TAG}" \
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
RUN  git clone "${SAMTOOLS_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${SAMTOOLS_TAG}" \
  && autoheader \
  && autoconf \
  && ./configure \
       --prefix="${SAMTOOLS_PREFIX}" \
       --with-htslib="${HTSLIB_PREFIX}" \
  && make -j $(grep -c ^processor /proc/cpuinfo) \
  && make -j $(grep -c ^processor /proc/cpuinfo) install


FROM ${IMAGE}

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

ARG HTSLIB_TAG
ARG BCFTOOLS_TAG
ARG SAMTOOLS_TAG

ARG HTSLIB_PREFIX_ARG="/opt/htslib/${HTSLIB_TAG}"
ARG BCFTOOLS_PREFIX_ARG="/opt/bcftools/${BCFTOOLS_TAG}"
ARG SAMTOOLS_PREFIX_ARG="/opt/samtools/${SAMTOOLS_TAG}"

ENV HTSLIB_PREFIX="${HTSLIB_PREFIX_ARG}"
ENV BCFTOOLS_PREFIX="${BCFTOOLS_PREFIX_ARG}"
ENV SAMTOOLS_PREFIX="${SAMTOOLS_PREFIX_ARG}"

COPY --from=builder "${HTSLIB_PREFIX}" "${HTSLIB_PREFIX}"
COPY --from=builder "${BCFTOOLS_PREFIX}" "${BCFTOOLS_PREFIX}"
COPY --from=builder "${SAMTOOLS_PREFIX}" "${SAMTOOLS_PREFIX}"

ENV PATH="${PATH}:${SAMTOOLS_PREFIX}/bin:${BCFTOOLS_PREFIX}/bin:${HTSLIB_PREFIX}/bin"
ENV CPATH="${HTSLIB_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LD_LIBRARY_PATH}"
