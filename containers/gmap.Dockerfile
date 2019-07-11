ARG IMAGE

FROM "${IMAGE}" as builder

ARG GMAP_VERSION="2019-05-12"
ARG GMAP_URL="http://research-pub.gene.com/gmap/src/gmap-gsnap-${GMAP_VERSION}.tar.gz"
ARG GMAP_PREFIX_ARG="/opt/gmap/${GMAP_VERSION}"
ENV GMAP_PREFIX="${GMAP_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       libbz2-dev \
       perl \
       wget \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && wget -O gmap.tar.gz "${GMAP_URL}" \
  && tar -zxf gmap.tar.gz \
  && rm gmap.tar.gz \
  && cd gmap* \
  && for level in none sse2 ssse3 avx2 sse4.1 sse4.2; \
     do \
       ./configure --prefix="${GMAP_PREFIX}" --with-simd-level="${level}"; \
       make; \
       make check; \
       make install; \
     done \
  && add_runtime_dep libbz2-1.0 perl zlib1g


FROM "${IMAGE}"

ARG GMAP_VERSION
ARG GMAP_PREFIX_ARG="/opt/gmap/${GMAP_VERSION}"
ENV GMAP_PREFIX="${GMAP_PREFIX_ARG}"

ENV PATH "${GMAP_PREFIX}/bin:${PATH}"

COPY --from=builder "${GMAP_PREFIX}" "${GMAP_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/gmap.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
