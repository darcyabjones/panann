ARG IMAGE

FROM "${IMAGE}" as builder

ARG JELLYFISH_VERSION
ARG JELLYFISH_URL="https://github.com/gmarcais/Jellyfish/releases/download/v${JELLYFISH_VERSION}/jellyfish-${JELLYFISH_VERSION}.tar.gz"
ARG JELLYFISH_PREFIX_ARG="/opt/jellyfish/${JELLYFISH_VERSION}"
ENV JELLYFISH_PREFIX="${JELLYFISH_PREFIX_ARG}"

# TODO figure out if swig bindings necessary for any tools depending
# on jellyfish.

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       build-essential \
       ca-certificates \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O jellyfish.tar.gz "${JELLYFISH_URL}" \
  && tar zxf jellyfish.tar.gz \
  && rm jellyfish.tar.gz \
  && cd jellyfish*/ \
  && ./configure --prefix="${JELLYFISH_PREFIX}" \
  && make \
  && make install

# No runtime dependencies!!

FROM ${IMAGE}

ARG JELLYFISH_VERSION
ARG JELLYFISH_PREFIX_ARG="/opt/jellyfish/${JELLYFISH_VERSION}"
ENV JELLYFISH_PREFIX="${JELLYFISH_PREFIX_ARG}"
LABEL jellyfish.version="${JELLYFISH_VERSION}"

ENV PATH="${JELLYFISH_PREFIX}/bin:${PATH}"
ENV INCLUDE="${JELLYFISH_DIR}/include:${INCLUDE}"
ENV CPATH="${JELLYFISH_DIR}/include:${CPATH}"
ENV LIBRARY_PATH="${JELLYFISH_DIR}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${JELLYFISH_DIR}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${JELLYFISH_DIR}/lib:${LD_RUN_PATH}"

COPY --from=builder "${JELLYFISH_PREFIX}" "${JELLYFISH_PREFIX}"
