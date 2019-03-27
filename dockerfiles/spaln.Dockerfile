ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

ENV SPALN_PREFIX="/opt/spaln"
ARG SPALN_VERSION="Ver.2.3.3"
ENV SPALN_REPO="https://github.com/ogotoh/spaln.git"

WORKDIR /tmp
RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       git \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp
RUN  git clone ${SPALN_REPO} . \
  && git fetch --tags \
  && git checkout tags/${SPALN_VERSION} \
  && cd src \
  && ./configure \
       --use_zlib=1 \
       --exec_prefix="${SPALN_PREFIX}/bin" \
       --table_dir="${SPALN_PREFIX}/table" \
       --alndbs_dir="${SPALN_PREFIX}/seqdb" \
  && make \
  && make install \
  && mv /tmp/perl ${SPALN_PREFIX}/perl


FROM debian:${DEBIAN_VERSION}

ENV SPALN_PREFIX="/opt/spaln"

COPY --from=builder ${SPALN_PREFIX} ${SPALN_PREFIX}

RUN  apt-get update \
  && apt-get install -y perl zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV PATH="${SPALN_PREFIX}:${PATH}"
