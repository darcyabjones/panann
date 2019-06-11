ARG IMAGE="darcyabjones/base"

FROM ${IMAGE} as builder

RUN  apt-get update \
  && apt-get install -y \
      build-essential \
      libtbb-dev \
      git \
      zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

# Runtime
#  libtbb2
#  zlib1g
#  perl - needs locales set
#  python2

ARG BOWTIE2_TAG
ARG BOWTIE2_PREFIX_ARG="/opt/bowtie2/${BOWTIE2_TAG}"
ARG BOWTIE2_REPO="https://github.com/BenLangmead/bowtie2.git"

ENV BOWTIE2_PREFIX="${BOWTIE2_PREFIX_ARG}"

WORKDIR /tmp/bowtie2
RUN  git clone ${BOWTIE2_REPO} . \
  && git fetch --tags \
  && git checkout "tags/${BOWTIE2_TAG}" \
  && make \
  && mkdir -p "${BOWTIE2_PREFIX}" \
  && mkdir -p "${BOWTIE2_PREFIX}/bin" \
  && find . -maxdepth 1 -type f -executable -exec mv {} "${BOWTIE2_PREFIX}/bin" \;


FROM ${IMAGE}

ARG BOWTIE2_TAG
ARG BOWTIE2_PREFIX_ARG="/opt/bowtie2/${BOWTIE2_TAG}"

ENV BOWTIE2_PREFIX="${BOWTIE2_PREFIX_ARG}"


RUN  apt-get update \
  && apt-get install -y --no-install-recommends \
       perl \
       libtbb2 \
       python \
       zlib1g \
  && rm -rf /var/lib/apt/lists/*

COPY --from=builder "${BOWTIE2_PREFIX}" "${BOWTIE2_PREFIX}"
ENV PATH="${PATH}:${BOWTIE2_PREFIX}/bin"
