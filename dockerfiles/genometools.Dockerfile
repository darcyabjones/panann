ARG DEBIAN_VERSION="stretch-20190228-slim"
FROM debian:${DEBIAN_VERSION} as builder

ENV GENOMETOOLS_PREFIX="/opt/genometools"
ARG GENOMETOOLS_VERSION="1.5.10"
ENV GENOMETOOLS_URL="http://genometools.org/pub/genometools-${GENOMETOOLS_VERSION}.tar.gz"


WORKDIR /tmp
RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       libcairo2-dev \
       libpango1.0-dev \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget ${GENOMETOOLS_URL} \
  && tar zxf genometools*.tar.gz \
  && rm genometools*.tar.gz \
  && cd genometools*/ \
  && make \
  && make prefix="${GENOMETOOLS_PREFIX}" install


FROM debian:${DEBIAN_VERSION}

RUN  apt-get update \
  && apt-get install -y \
       libcairo2 \
       libpango-1.0-0 \
       libpangocairo-1.0-0 \
  && rm -rf /var/lib/apt/lists/*

ENV GENOMETOOLS_PREFIX="/opt/genometools"

COPY --from=builder ${GENOMETOOLS_PREFIX} ${GENOMETOOLS_PREFIX}


ENV PATH="${GENOMETOOLS_PREFIX}/bin:${PATH}"
ENV INCLUDE="${GENOMETOOLS_PREFIX}/include:${INCLUDE}"
ENV CPATH="${GENOMETOOLS_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${GENOMETOOLS_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${GENOMETOOLS_PREFIX}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${GENOMETOOLS_PREFIX}/lib:${LD_RUN_PATH}"
