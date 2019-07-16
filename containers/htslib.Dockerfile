ARG IMAGE="darcyabjones/base"

FROM ${IMAGE} as builder

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

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       autoconf \
       build-essential \
       ca-certificates \
       git \
       libbz2-dev \
       libcurl4-gnutls-dev \
       libgsl-dev \
       liblzma-dev \
       libncurses5-dev \
       libssl-dev \
       libperl-dev \
       zlib1g-dev \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates

# CA cert stuff required for git clone https

WORKDIR /tmp/htslib
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && git clone "${HTSLIB_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${HTSLIB_TAG}" \
  && autoheader \
  && autoconf \
  && ./configure --prefix="${HTSLIB_PREFIX}" --enable-libcurl \
  && make -j $(grep -c ^processor /proc/cpuinfo) \
  && make -j $(grep -c ^processor /proc/cpuinfo) install \
  && add_runtime_dep libbz2-1.0 libcurl3-gnutls libssl1.1 lzma zlib1g
# htslib: libssl (providing libcrypto) only required if need Amazon S3 support.

WORKDIR /tmp/bcftools
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && git clone "${BCFTOOLS_REPO}" . \
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
  && make -j $(grep -c ^processor /proc/cpuinfo) install \
  && add_runtime_dep libperl5.28 libgsl23 zlib1g
# Samtools also depends on htslib and it's dependencies.

WORKDIR /tmp/samtools
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && git clone "${SAMTOOLS_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${SAMTOOLS_TAG}" \
  && autoheader \
  && autoconf \
  && ./configure \
       --prefix="${SAMTOOLS_PREFIX}" \
       --with-htslib="${HTSLIB_PREFIX}" \
  && make -j $(grep -c ^processor /proc/cpuinfo) \
  && make -j $(grep -c ^processor /proc/cpuinfo) install \
  && add_runtime_dep libncurses5 zlib1g
# Samtools also depends on htslib and it's dependencies.


FROM ${IMAGE}

ARG HTSLIB_TAG
ARG BCFTOOLS_TAG
ARG SAMTOOLS_TAG

ARG HTSLIB_PREFIX_ARG="/opt/htslib/${HTSLIB_TAG}"
ARG BCFTOOLS_PREFIX_ARG="/opt/bcftools/${BCFTOOLS_TAG}"
ARG SAMTOOLS_PREFIX_ARG="/opt/samtools/${SAMTOOLS_TAG}"

ENV HTSLIB_PREFIX="${HTSLIB_PREFIX_ARG}"
ENV BCFTOOLS_PREFIX="${BCFTOOLS_PREFIX_ARG}"
ENV SAMTOOLS_PREFIX="${SAMTOOLS_PREFIX_ARG}"

LABEL htslib.version="${HTSLIB_TAG}"
LABEL bcftools.version="${BCFTOOLS_TAG}"
LABEL samtools.version="${SAMTOOLS_TAG}"

ENV PATH="${SAMTOOLS_PREFIX}/bin:${BCFTOOLS_PREFIX}/bin:${HTSLIB_PREFIX}/bin:${PATH}"
ENV CPATH="${HTSLIB_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LD_LIBRARY_PATH}"

COPY --from=builder "${HTSLIB_PREFIX}" "${HTSLIB_PREFIX}"
COPY --from=builder "${BCFTOOLS_PREFIX}" "${BCFTOOLS_PREFIX}"
COPY --from=builder "${SAMTOOLS_PREFIX}" "${SAMTOOLS_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/htslib.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
