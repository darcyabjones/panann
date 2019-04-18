ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

ENV BWA_PREFIX="/opt/bwa"
ARG BWA_VERSION="v0.7.17"
ENV BWA_URL="https://github.com/lh3/bwa"

WORKDIR /tmp
RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       git \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && git clone ${BWA_URL} . \
  && git fetch --tags \
  && git checkout tags/${BWA_VERSION} \
  && make \
  && mkdir ${BWA_PREFIX} \
  && mv bwa qualfa2fq.pl xa2multi.pl ${BWA_PREFIX}

FROM debian:${DEBIAN_VERSION}

ENV BWA_PREFIX="/opt/bwa"
COPY --from=builder ${BWA_PREFIX} ${BWA_PREFIX}

RUN  apt-get update \
  && apt-get install -y perl zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV PATH="${BWA_PREFIX}:${PATH}"
