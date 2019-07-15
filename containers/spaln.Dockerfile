ARG IMAGE

FROM "${IMAGE}" as builder

ARG SPALN_TAG
ARG SPALN_REPO="https://github.com/ogotoh/spaln.git"
ARG SPALN_PREFIX_ARG="/opt/spaln/${SPALN_TAG}"
ENV SPALN_PREFIX="${SPALN_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       build-essential \
       ca-certificates \
       git \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${SPALN_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${SPALN_TAG}" \
  && cd src \
  && ./configure \
       --use_zlib=1 \
       --exec_prefix="${SPALN_PREFIX}/bin" \
       --table_dir="${SPALN_PREFIX}/table" \
       --alndbs_dir="${SPALN_PREFIX}/seqdb" \
  && make \
  && make install \
  && mv /tmp/perl "${SPALN_PREFIX}/perl" \
  && add_runtime_dep perl zlib1g

# CA cert stuff required for git clone https


FROM "${IMAGE}"

ARG SPALN_TAG
ARG SPALN_PREFIX_ARG="/opt/spaln/${SPALN_TAG}"
ENV SPALN_PREFIX="${SPALN_PREFIX_ARG}"
ENV ALN_TAB="${SPALN_PREFIX}/table"
ENV ALN_DBS="${SPALN_PREFIX}/seqdb"
LABEL spaln.version="${SPALN_TAG}"

ENV PATH="${SPALN_PREFIX}/bin:${SPALN_PREFIX}/perl:${PATH}"

COPY --from=builder "${SPALN_PREFIX}" "${SPALN_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/spaln.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
