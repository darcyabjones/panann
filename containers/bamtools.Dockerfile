ARG IMAGE

FROM "${IMAGE}" as bamtools_builder

ARG BAMTOOLS_TAG="v2.5.1"
ARG BAMTOOLS_REPO="https://github.com/pezmaster31/bamtools.git"
ARG BAMTOOLS_PREFIX_ARG="/opt/bamtools/${BAMTOOLS_TAG}"
ENV BAMTOOLS_PREFIX="${BAMTOOLS_PREFIX_ARG}"

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
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${BAMTOOLS_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${BAMTOOLS_TAG}" \
  && mkdir -p build \
  && cd build \
  && cmake -DCMAKE_INSTALL_PREFIX="${BAMTOOLS_PREFIX}" .. \
  && make \
  && make install \
  && add_runtime_dep zlib1g


FROM "${IMAGE}"

LABEL maintainer="darcy.ab.jones@gmail.com"

ARG BAMTOOLS_TAG="v2.5.1"
ARG BAMTOOLS_PREFIX_ARG="/opt/bamtools/${BAMTOOLS_TAG}"
ENV BAMTOOLS_PREFIX="${BAMTOOLS_PREFIX_ARG}"
LABEL bamtools.version="${BAMTOOLS_TAG}"

ENV PATH="${BAMTOOLS_PREFIX}/bin:${PATH}"
ENV LIBRARY_PATH="${LD_LIBRARY_PATH}:${BAMTOOLS_PREFIX}/lib"
ENV CPATH="${CPATH}:${BAMTOOLS_PREFIX}/include"

COPY --from=bamtools_builder "${BAMTOOLS_PREFIX}" "${BAMTOOLS_PREFIX}"
COPY --from=bamtools_builder "${APT_REQUIREMENTS_FILE}" /build/apt/bamtools.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"

WORKDIR /
