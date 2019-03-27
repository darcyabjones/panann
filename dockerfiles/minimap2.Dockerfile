ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

ENV MINIMAP_PREFIX="/opt/minimap"
ARG MINIMAP_VERSION="v2.16"
ENV MINIMAP_REPO="https://github.com/lh3/minimap2.git"

ENV K8_PREFIX="/opt/k8"
ARG K8_VERSION="0.2.4"
ENV K8_URL="https://github.com/attractivechaos/k8/releases/download/v${K8_VERSION}/k8-${K8_VERSION}.tar.bz2"

WORKDIR /tmp
RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       curl \
       git \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && git clone ${MINIMAP_REPO} . \
  && git fetch --tags \
  && git checkout tags/${MINIMAP_VERSION} \
  && make extra \
  && mkdir ${MINIMAP_PREFIX} \
  && cp minimap2 sdust minimap2-lite ${MINIMAP_PREFIX} \
  && cp misc/paftools.js ${MINIMAP_PREFIX} \
  && curl -L ${K8_URL} \
   | tar jxf - \
  && mkdir ${K8_PREFIX} \
  && mv "k8-${K8_VERSION}/k8-$(uname -s)" ${K8_PREFIX}/k8


FROM debian:${DEBIAN_VERSION}

ENV MINIMAP_PREFIX="/opt/minimap"
ENV K8_PREFIX="/opt/k8"

COPY --from=builder ${MINIMAP_PREFIX} ${MINIMAP_PREFIX}
COPY --from=builder ${K8_PREFIX} ${K8_PREFIX}

RUN  apt-get update \
  && apt-get install -y zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV PATH="${MINIMAP_PREFIX}:${K8_PREFIX}:${PATH}"
