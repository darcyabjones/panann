ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

ENV STAR_PREFIX="/opt/star"
ARG STAR_VERSION="2.7.0e"
ENV STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz"

WORKDIR /tmp
RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       wget \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && wget ${STAR_URL} \
  && tar zxf ${STAR_VERSION}.tar.gz \
  && cd STAR-${STAR_VERSION}/source \
  && make STAR STARlong install \
  && mkdir ${STAR_PREFIX} \
  && mv ../bin/STAR ../bin/STARlong ${STAR_PREFIX}


FROM debian:${DEBIAN_VERSION}

ENV STAR_PREFIX="/opt/star"
COPY --from=builder ${STAR_PREFIX} ${STAR_PREFIX}

RUN  apt-get update \
  && apt-get install -y libgomp1 zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV PATH="${STAR_PREFIX}:${PATH}"
