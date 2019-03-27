ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

ENV JELLYFISH_PREFIX="/opt/jellyfish"
ARG JELLYFISH_VERSION="2.2.10"
ENV JELLYFISH_URL="https://github.com/gmarcais/Jellyfish/releases/download/v${JELLYFISH_VERSION}/jellyfish-${JELLYFISH_VERSION}.tar.gz"


WORKDIR /tmp
RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget ${JELLYFISH_URL} \
  && tar zxf jellyfish*.tar.gz \
  && cd jellyfish*/
