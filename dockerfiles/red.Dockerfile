ARG DEBIAN_VERSION="buster-20190228-slim"

# Need buster because of GCC8 dependency.

FROM debian:${DEBIAN_VERSION} as builder

ENV RED_PREFIX="/opt/red"
ARG RED_VERSION="v2.0"
ENV RED_REPO="https://github.com/TulsaBioinformaticsToolsmith/Red.git"

WORKDIR /tmp
RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       git \
  && rm -rf /var/lib/apt/lists/* \
  && git clone ${RED_REPO} . \
  && git fetch --tags \
  && git checkout tags/${RED_VERSION} \
  && cd src_2.0 \
  && sed -i 's/CXX = g++-8/CXX = g++/g' Makefile \
  && make bin \
  && make \
  && mkdir ${RED_PREFIX} \
  && mv ../bin/Red ${RED_PREFIX}


FROM debian:${DEBIAN_VERSION}
ENV RED_PREFIX="/opt/red"

COPY --from=builder ${RED_PREFIX} ${RED_PREFIX}

RUN  apt-get update \
  && apt-get install -y libgomp1 \
  && rm -rf /var/lib/apt/lists/*

ENV PATH="${RED_PREFIX}:${PATH}"
