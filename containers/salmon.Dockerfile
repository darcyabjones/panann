ARG IMAGE

FROM "${IMAGE}" as builder

ARG SALMON_TAG="v0.13.1"
ARG SALMON_REPO="https://github.com/COMBINE-lab/salmon.git"
ARG SALMON_PREFIX_ARG="/opt/salmon/${SALMON_TAG}"
ENV SALMON_PREFIX="${SALMON_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && . /build/base.sh \
  && echo "deb http://deb.debian.org/debian stretch-backports main" >> /etc/apt/sources.list \
  && apt-get update \
  && apt-get install -y \
       autoconf \
       build-essential \
       curl \
       git \
       libbz2-dev \
       libboost-all-dev \
       liblzma-dev \
       libtbb-dev \
       unzip \
       zlib1g-dev \
  && apt-get -t stretch-backports install -y cmake \
  && rm -rf /var/lib/apt/lists/* \
  && git clone "${SALMON_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${SALMON_TAG}" \
  && mkdir -p build \
  && cd build \
  && cmake -DCMAKE_INSTALL_PREFIX="${SALMON_PREFIX}" .. \
  && make \
  && make install \
  && make test \
  && prepend_path PATH "\${SALMON_PREFIX}/bin" \
  && prepend_path LD_LIBRARY_PATH "\${SALMON_PREFIX}/lib" \
  && add_runtime_dep libtbb2 libbz2-1.0 zlib1g liblzma5


FROM "${IMAGE}"

ARG SALMON_TAG
ARG SALMON_PREFIX_ARG="/opt/salmon/${SALMON_TAG}"
ENV SALMON_PREFIX="${SALMON_PREFIX_ARG}"

COPY --from=builder "${SALMON_PREFIX}" "${SALMON_PREFIX}"
COPY --from=builder "/build/apt-requirements.txt" "/build/apt/salmon.txt"
COPY --from=builder "/build/env.sh" "/build/env/salmon.sh"

RUN  set -eu \
  && . /build/base.sh \
  && for f in /build/env/*.sh; do . ${f}; done \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/*
