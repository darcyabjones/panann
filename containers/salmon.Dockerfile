ARG IMAGE

FROM "${IMAGE}" as builder

ARG SALMON_TAG="v0.13.1"
ARG SALMON_REPO="https://github.com/COMBINE-lab/salmon.git"
ARG SALMON_PREFIX_ARG="/opt/salmon/${SALMON_TAG}"
ENV SALMON_PREFIX="${SALMON_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       autoconf \
       build-essential \
       ca-certificates \
       cmake \
       curl \
       git \
       libbz2-dev \
       libboost-all-dev \
       liblzma-dev \
       libtbb-dev \
       unzip \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${SALMON_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${SALMON_TAG}" \
  && mkdir -p build \
  && cd build \
  && cmake -DCMAKE_INSTALL_PREFIX="${SALMON_PREFIX}" .. \
  && make \
  && make install \
  && make test \
  && add_runtime_dep libtbb2 libbz2-1.0 libgomp1 zlib1g liblzma5

# CA cert stuff sometimes required for git clone https


FROM "${IMAGE}"

ARG SALMON_TAG
ARG SALMON_PREFIX_ARG="/opt/salmon/${SALMON_TAG}"
ENV SALMON_PREFIX="${SALMON_PREFIX_ARG}"
LABEL salmon.version="${SALMON_TAG}"

ENV PATH "${SALMON_PREFIX}/bin:${PATH}"
ENV LD_LIBRARY_PATH "${SALMON_PREFIX}/lib:${LD_LIBRARY_PATH}"

COPY --from=builder "${SALMON_PREFIX}" "${SALMON_PREFIX}"
COPY --from=builder "/build/apt-requirements.txt" /build/apt/salmon.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
