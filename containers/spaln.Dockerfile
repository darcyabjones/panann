ARG IMAGE

FROM ${IMAGE} as builder

ARG SPALN_TAG
ARG SPALN_REPO="https://github.com/ogotoh/spaln.git"
ARG SPALN_PREFIX_ARG="/opt/spaln/${SPALN_TAG}"
ENV SPALN_PREFIX="${SPALN_PREFIX_ARG}"

WORKDIR /tmp
RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       git \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

# Runtime requires
#   perl
#   zlib1g

WORKDIR /tmp
RUN  git clone "${SPALN_REPO}" . \
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
  && mv /tmp/perl "${SPALN_PREFIX}/perl"


FROM ${IMAGE}

ARG SPALN_TAG
ARG SPALN_PREFIX_ARG="/opt/spaln/${SPALN_TAG}"
ENV SPALN_PREFIX="${SPALN_PREFIX_ARG}"
ENV ALN_TAB="${SPALN_PREFIX}/table"
ENV ALN_DBS="${SPALN_PREFIX}/seqdb"

COPY --from=builder "${SPALN_PREFIX}" "${SPALN_PREFIX}"

RUN  apt-get update \
  && apt-get install -y --no-install-recommends perl zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV PATH="${SPALN_PREFIX}/bin:${SPALN_PREFIX}/perl:${PATH}"
