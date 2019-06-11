ARG IMAGE

FROM ${IMAGE} as builder

ARG MINIMAP_TAG
ARG MINIMAP_REPO="https://github.com/lh3/minimap2.git"
ARG MINIMAP_PREFIX_ARG="/opt/minimap/${MINIMAP_TAG}"
ENV MINIMAP_PREFIX="${MINIMAP_PREFIX_ARG}"

ARG K8_VERSION
ARG K8_URL="https://github.com/attractivechaos/k8/releases/download/v${K8_VERSION}/k8-${K8_VERSION}.tar.bz2"
ARG K8_PREFIX_ARG="/opt/k8/${K8_VERSION}"
ENV K8_PREFIX="${K8_PREFIX_ARG}"

WORKDIR /tmp
RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       curl \
       git \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && git clone "${MINIMAP_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${MINIMAP_TAG}" \
  && make extra \
  && mkdir -p "${MINIMAP_PREFIX}" \
  && cp minimap2 sdust minimap2-lite "${MINIMAP_PREFIX}" \
  && cp misc/paftools.js "${MINIMAP_PREFIX}" \
  && curl -L "${K8_URL}" \
   | tar jxf - \
  && mkdir -p "${K8_PREFIX}" \
  && mv "k8-${K8_VERSION}/k8-$(uname -s)" "${K8_PREFIX}/k8"


FROM ${IMAGE}

ARG MINIMAP_TAG
ARG MINIMAP_PREFIX_ARG="/opt/minimap/${MINIMAP_TAG}"
ENV MINIMAP_PREFIX="${MINIMAP_PREFIX_ARG}"

ARG K8_VERSION
ARG K8_PREFIX_ARG="/opt/k8/${K8_VERSION}"
ENV K8_PREFIX="${K8_PREFIX_ARG}"

COPY --from=builder "${MINIMAP_PREFIX}" "${MINIMAP_PREFIX}"
COPY --from=builder "${K8_PREFIX}" "${K8_PREFIX}"

RUN  apt-get update \
  && apt-get install -y zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV PATH="${MINIMAP_PREFIX}:${K8_PREFIX}:${PATH}"
