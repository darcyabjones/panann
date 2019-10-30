ARG IMAGE

FROM "${IMAGE}" as sibeliaz_builder

ARG SIBELIAZ_COMMIT
ARG SIBELIAZ_REPO="https://github.com/medvedevgroup/SibeliaZ.git"
ARG SIBELIAZ_PREFIX_ARG="/opt/sibeliaz/${SIBELIAZ_COMMIT}"
ENV SIBELIAZ_PREFIX="${SIBELIAZ_PREFIX_ARG}"


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
       libtbb-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${SIBELIAZ_REPO}" . \
  && git fetch --tags \
  && git checkout "${SIBELIAZ_COMMIT}" \
  && mkdir build \
  && git submodule update --init --recursive \
  && cd build \
  && cmake .. -DCMAKE_INSTALL_PREFIX="${SIBELIAZ_PREFIX}" \
  && make install \
  && add_runtime_dep libgomp1 libtbb2 time


FROM "${IMAGE}"

ARG SIBELIAZ_COMMIT
ARG SIBELIAZ_PREFIX_ARG
ENV SIBELIAZ_PREFIX="${SIBELIAZ_PREFIX_ARG}"
LABEL sibeliaz.version="${SIBELIAZ_COMMIT}"

ENV PATH="${PATH}:${SIBELIAZ_PREFIX}/bin"

COPY --from=sibeliaz_builder "${SIBELIAZ_PREFIX}" "${SIBELIAZ_PREFIX}"
COPY --from=sibeliaz_builder "${APT_REQUIREMENTS_FILE}" /build/apt/sibeliaz.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"

WORKDIR /
