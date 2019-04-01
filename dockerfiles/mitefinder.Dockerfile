ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

ENV MITEFINDER_PREFIX="/opt/mitefinder"
ARG MITEFINDER_VERSION="833754b0ff1899e8cb0f260d6d5011d3583b3012"
ENV MITEFINDER_REPO="https://github.com/screamer/miteFinder.git"

WORKDIR /tmp
RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       git \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp
RUN  git clone ${MITEFINDER_REPO} . \
  && git checkout ${MITEFINDER_VERSION} \
  && make \
  && mkdir ${MITEFINDER_PREFIX} \
  && mkdir ${MITEFINDER_PREFIX}/bin \
  && mv miteFinder ${MITEFINDER_PREFIX}/bin \
  && mv profile ${MITEFINDER_PREFIX}


FROM debian:${DEBIAN_VERSION}

ENV MITEFINDER_PREFIX="/opt/mitefinder"
COPY --from=builder ${MITEFINDER_PREFIX} ${MITEFINDER_PREFIX}

ENV PATH="${MITEFINDER_PREFIX}/bin"
ENV MITEFINDER_PROFILE="${MITEFINDER_PREFIX}/profile/pattern_scoring.txt"
