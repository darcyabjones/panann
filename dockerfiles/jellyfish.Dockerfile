ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

ENV JELLYFISH_PREFIX="/opt/jellyfish"
ARG JELLYFISH_VERSION="2.2.10"
ENV JELLYFISH_URL="https://github.com/gmarcais/Jellyfish/releases/download/v${JELLYFISH_VERSION}/jellyfish-${JELLYFISH_VERSION}.tar.gz"

# TODO figure out if swig bindings necessary for any tools depending
# on jellyfish.

WORKDIR /tmp
RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget ${JELLYFISH_URL} \
  && tar zxf jellyfish*.tar.gz \
  && rm jellyfish*.tar.gz \
  && cd jellyfish*/ \
  && ./configure --prefix=${JELLYFISH_PREFIX} \
  && make \
  && make install


FROM debian:${DEBIAN_VERSION}

ENV JELLYFISH_PREFIX="/opt/jellyfish"

COPY --from=builder ${JELLYFISH_PREFIX} ${JELLYFISH_PREFIX}

ENV PATH="${JELLYFISH_PREFIX}:${PATH}"
ENV INCLUDE="${JELLYFISH_DIR}/include:${INCLUDE}"
ENV CPATH="${JELLYFISH_DIR}/include:${CPATH}"
ENV LIBRARY_PATH="${JELLYFISH_DIR}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${JELLYFISH_DIR}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${JELLYFISH_DIR}/lib:${LD_RUN_PATH}"
