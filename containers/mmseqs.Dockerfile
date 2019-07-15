ARG IMAGE

FROM "${IMAGE}" as builder

## Config variables
ARG MMSEQS_TAG="7-4e23d"
ARG MMSEQS_CMAKE_OPTIONS=""
ARG MMSEQS_REPO="https://github.com/soedinglab/MMseqs2.git"
ARG MMSEQS_PREFIX_ARG="/opt/mmseqs/${MMSEQS_TAG}"
ENV MMSEQS_PREFIX="${MMSEQS_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       cmake \
       git \
       libbz2-dev \
       libmpich-dev \
       xxd \
       zlib1g-dev \
  && rm -rf -- /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${MMSEQS_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${MMSEQS_TAG}" \
  && git submodule update --init \
  && mkdir -p build \
  && cd build \
  && cmake \
       ${MMSEQS_CMAKE_OPTIONS} \
       -DHAVE_MPI=1 \
       -DCMAKE_BUILD_TYPE=Release \
       -DCMAKE_INSTALL_PREFIX="${MMSEQS_PREFIX}" .. \
  && make \
  && make install \
  && add_runtime_dep \
       gawk \
       bash \
       grep \
       libbz2-1.0 \
       libgomp1 \
       libstdc++6 \
       mpich \
       zlib1g


FROM "${IMAGE}"

LABEL maintainer="darcy.ab.jones@gmail.com"

ARG MMSEQS_TAG="7-4e23d"
ARG MMSEQS_PREFIX_ARG="/opt/mmseqs/${MMSEQS_TAG}"
ENV MMSEQS_PREFIX="${MMSEQS_PREFIX_ARG}"
LABEL mmseqs.version="${MMSEQS_TAG}"

ENV PATH="${MMSEQS_PREFIX}/bin:${PATH}"

COPY --from=builder "${MMSEQS_PREFIX}" "${MMSEQS_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/mmseqs.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
