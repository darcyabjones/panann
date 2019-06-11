ARG IMAGE

FROM ${IMAGE} as builder

ARG GENOMETOOLS_VERSION
ARG GENOMETOOLS_URL="http://genometools.org/pub/genometools-${GENOMETOOLS_VERSION}.tar.gz"
ARG GENOMETOOLS_PREFIX_ARG="/opt/genometools/${GENOMETOOLS_VERSION}"
ENV GENOMETOOLS_PREFIX="${GENOMETOOLS_PREFIX_ARG}"

WORKDIR /tmp
RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       libcairo2-dev \
       libpango1.0-dev \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget -O genometools.tar.gz "${GENOMETOOLS_URL}" \
  && tar zxf genometools.tar.gz \
  && rm genometools.tar.gz \
  && cd genometools*/ \
  && make \
  && make prefix="${GENOMETOOLS_PREFIX}" install

# Runtime requires
# libcairo2 \
# libpango-1.0-0 \
# libpangocairo-1.0-0 \


FROM ${IMAGE}

RUN  apt-get update \
  && apt-get install -y \
       libcairo2 \
       libpango-1.0-0 \
       libpangocairo-1.0-0 \
  && rm -rf /var/lib/apt/lists/*

ARG GENOMETOOLS_VERSION
ARG GENOMETOOLS_PREFIX_ARG="/opt/genometools/${GENOMETOOLS_VERSION}"
ENV GENOMETOOLS_PREFIX="${GENOMETOOLS_PREFIX_ARG}"

COPY --from=builder "${GENOMETOOLS_PREFIX}" "${GENOMETOOLS_PREFIX}"

ENV PATH="${GENOMETOOLS_PREFIX}/bin:${PATH}"
ENV INCLUDE="${GENOMETOOLS_PREFIX}/include:${INCLUDE}"
ENV CPATH="${GENOMETOOLS_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${GENOMETOOLS_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${GENOMETOOLS_PREFIX}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${GENOMETOOLS_PREFIX}/lib:${LD_RUN_PATH}"
