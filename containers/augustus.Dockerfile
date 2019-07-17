ARG IMAGE
ARG HTSLIB_IMAGE

FROM "${HTSLIB_IMAGE}" as htslib_builder

FROM "${IMAGE}" as builder

ARG AUGUSTUS_COMMIT="8b1b14a7489e4545e89c8725dc33268f6c2a9117"
ARG AUGUSTUS_REPO="https://github.com/Gaius-Augustus/Augustus.git"
ARG AUGUSTUS_PREFIX_ARG="/opt/augustus/${AUGUSTUS_COMMIT}"
ENV AUGUSTUS_PREFIX="${AUGUSTUS_PREFIX_ARG}"

ARG HTSLIB_TAG="1.9"
ARG HTSLIB_REPO="https://github.com/samtools/htslib.git"

ARG SAMTOOLS_TAG="1.9"
ARG SAMTOOLS_REPO="https://github.com/samtools/samtools.git"

# Install htslib and samtools.
# Bam2wig in augustus depends on some intermediate samtools/htslib compilation
# rather than the actual headers/shared libraries, so I have to compile it
# separately.

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && apt-get update \
  && apt-get install -y \
       autoconf \
       build-essential \
       ca-certificates \
       git \
       libbz2-dev \
       libcurl4-gnutls-dev \
       liblzma-dev \
       libncurses5-dev \
       libssl-dev \
       zlib1g-dev \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates

WORKDIR /tmp/htslib
RUN  git clone ${HTSLIB_REPO} . \
  && git fetch --tags \
  && git checkout "tags/${HTSLIB_TAG}" \
  && autoheader \
  && autoconf \
  && ./configure --enable-libcurl \
  && make -j $(grep -c ^processor /proc/cpuinfo)

WORKDIR /tmp/samtools
RUN  git clone ${SAMTOOLS_REPO} . \
  && git fetch --tags \
  && git checkout "tags/${SAMTOOLS_TAG}" \
  && autoheader \
  && autoconf \
  && ./configure \
  && make -j $(grep -c ^processor /proc/cpuinfo)

# Install augustus

# This is for bam2wig
ENV TOOLDIR="/tmp"

WORKDIR /tmp/augustus
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       libbamtools-dev \
       libboost-all-dev \
       libboost-iostreams-dev \
       libboost-graph-dev \
       libcurl4-gnutls-dev \
       libgsl-dev \
       liblpsolve55-dev \
       libssl-dev \
       libsuitesparse-dev \
       libsqlite3-dev \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${AUGUSTUS_REPO}" . \
  && git fetch --tags \
  && git checkout "${AUGUSTUS_COMMIT}" \
  && mkdir bin \
  && sed -i "s/# SQLITE = true/SQLITE = true/g" common.mk \
  && sed -i "s/# COMPGENEPRED = true/COMPGENEPRED = true/g" common.mk \
  && sed -i 's~INSTALLDIR = .*~INSTALLDIR="${AUGUSTUS_PREFIX}"~g' Makefile \
  && cd auxprogs/bam2wig \
  && make \
  && cd /tmp/augustus \
  && make \
  && make install \
  && make test \
  && add_runtime_dep \
       libbamtools2.5.1 \
       libcurl3-gnutls \
       libgsl23 \
       libparallel-forkmanager-perl \
       libssl1.1 \
       libsqlite3-0 \
       lp-solve \
       zlib1g


FROM "${IMAGE}"

ARG AUGUSTUS_COMMIT
ARG AUGUSTUS_PREFIX_ARG="/opt/augustus/${AUGUSTUS_COMMIT}"
ENV AUGUSTUS_PREFIX="${AUGUSTUS_PREFIX_ARG}"
LABEL augustus.version="${AUGUSTUS_COMMIT}"

ENV PATH="${AUGUSTUS_PREFIX}/bin:${AUGUSTUS_PREFIX}/scripts:${PATH}"
ENV AUGUSTUS_CONFIG_PATH="${AUGUSTUS_PREFIX}/config"

COPY --from=builder "${AUGUSTUS_PREFIX}" "${AUGUSTUS_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/augustus.txt

ARG HTSLIB_TAG
ARG SAMTOOLS_TAG
ARG HTSLIB_PREFIX_ARG="/opt/htslib/${HTSLIB_TAG}"
ARG SAMTOOLS_PREFIX_ARG="/opt/samtools/${SAMTOOLS_TAG}"
ENV HTSLIB_PREFIX="${HTSLIB_PREFIX_ARG}"
ENV SAMTOOLS_PREFIX="${SAMTOOLS_PREFIX_ARG}"
LABEL htslib.version="${HTSLIB_TAG}"
LABEL samtools.version="${SAMTOOLS_TAG}"

ENV PATH="${SAMTOOLS_PREFIX}/bin:${BCFTOOLS_PREFIX}/bin:${HTSLIB_PREFIX}/bin:${PATH}"
ENV CPATH="${HTSLIB_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LD_LIBRARY_PATH}"

COPY --from=htslib_builder "${HTSLIB_PREFIX}" "${HTSLIB_PREFIX}"
COPY --from=htslib_builder "${SAMTOOLS_PREFIX}" "${SAMTOOLS_PREFIX}"
COPY --from=htslib_builder "${APT_REQUIREMENTS_FILE}" /build/apt/htslib.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"

# This is useful for testing.
# COPY --from=builder "/tmp/augustus/examples" "${AUGUSTUS_PREFIX}/examples"
# RUN augustus --species=human --UTR=on ${AUGUSTUS_PREFIX}/examples/example.fa
