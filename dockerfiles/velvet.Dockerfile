ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder
ENV VELVET_PREFIX="/opt/velvet"
ARG VELVET_VERSION="v1.2.10"
ENV VELVET_URL="https://github.com/dzerbino/velvet/archive/${VELVET_VERSION}.tar.gz"

WORKDIR /tmp
RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       wget \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && wget ${VELVET_URL} \
  && tar zxf ${VELVET_VERSION}.tar.gz \
  && cd velvet* \
  && make 'OPENMP=1' 'MAXKMERLENGTH=81' 'CATEGORIES=5' \
  && mkdir ${VELVET_PREFIX} \
  && cp velveth velvetg ${VELVET_PREFIX}


FROM debian:${DEBIAN_VERSION}

ENV VELVET_PREFIX="/opt/velvet"
COPY --from=builder ${VELVET_PREFIX} ${VELVET_PREFIX}

RUN  apt-get update \
  && apt-get install -y libgomp1 zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV PATH="${VELVET_PREFIX}:${PATH}"
