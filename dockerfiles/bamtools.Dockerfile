ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder


RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       cmake \
       git \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

ENV BAMTOOLS_PREFIX="/opt/bamtools"
ARG BAMTOOLS_VERSION="v2.5.1"
ENV BAMTOOLS_REPO="https://github.com/pezmaster31/bamtools.git"

WORKDIR /tmp
RUN  git clone ${BAMTOOLS_REPO} . \
  && git fetch --tags \
  && git checkout tags/${BAMTOOLS_VERSION}

WORKDIR /tmp/build
RUN cmake -DCMAKE_INSTALL_PREFIX=${BAMTOOLS_PREFIX} .. \
  && make \
  && make install


FROM debian:${DEBIAN_VERSION}

LABEL maintainer="darcy.ab.jones@gmail.com"

RUN  apt-get update \
  && apt-get install -y \
       zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV BAMTOOLS_PREFIX="/opt/bamtools"

COPY --from=builder ${BAMTOOLS_PREFIX} ${BAMTOOLS_PREFIX}

ENV PATH="${BAMTOOLS_PREFIX}/bin:${PATH}"
ENV LIBRARY_PATH="${LD_LIBRARY_PATH}:${BAMTOOLS_PREFIX}/lib"
ENV CPATH="${CPATH}:${BAMTOOLS_PREFIX}/include"

ENTRYPOINT ["bamtools"]
CMD ["--help"]
