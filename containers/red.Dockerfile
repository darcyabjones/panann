ARG BUILD_IMAGE
ARG IMAGE

# Need edge distro for build of GCC8 dependency.

FROM "${BUILD_IMAGE}" as builder

ARG RED_TAG
ARG RED_REPO="https://github.com/TulsaBioinformaticsToolsmith/Red.git"
ARG RED_PREFIX_ARG="/opt/red/${RED_TAG}"
ENV RED_PREFIX="${RED_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       git \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${RED_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${RED_TAG}" \
  && cd src_2.0 \
  && sed -i 's/CXX = g++-8/CXX = g++/g' Makefile \
  && make bin \
  && make \
  && mkdir -p "${RED_PREFIX}/bin" \
  && mv ../bin/Red "${RED_PREFIX}/bin" \
  && add_runtime_dep libgomp1

# CA cert stuff sometime required for git clone https

FROM "${IMAGE}"

ARG RED_TAG
ARG RED_PREFIX_ARG="/opt/red/${RED_TAG}"
ENV RED_PREFIX="${RED_PREFIX_ARG}"

ENV PATH="${RED_PREFIX}/bin:${PATH}"

COPY --from=builder "${RED_PREFIX}" "${RED_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/red.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
