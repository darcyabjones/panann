ARG IMAGE="darcyabjones/base"

FROM "${IMAGE}" as builder

ARG BOWTIE2_TAG
ARG BOWTIE2_PREFIX_ARG="/opt/bowtie2/${BOWTIE2_TAG}"
ARG BOWTIE2_REPO="https://github.com/BenLangmead/bowtie2.git"

ENV BOWTIE2_PREFIX="${BOWTIE2_PREFIX_ARG}"

WORKDIR /tmp/bowtie2
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       build-essential \
       ca-certificates \
       libtbb-dev \
       git \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${BOWTIE2_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${BOWTIE2_TAG}" \
  && make \
  && mkdir -p "${BOWTIE2_PREFIX}" \
  && mkdir -p "${BOWTIE2_PREFIX}/bin" \
  && find . -maxdepth 1 -type f -executable -exec mv {} "${BOWTIE2_PREFIX}/bin" \; \
  && add_runtime_dep libtbb2 perl python zlib1g

# CA cert stuff sometimes required for git clone https


FROM "${IMAGE}"

ARG BOWTIE2_TAG
ARG BOWTIE2_PREFIX_ARG="/opt/bowtie2/${BOWTIE2_TAG}"
ENV BOWTIE2_PREFIX="${BOWTIE2_PREFIX_ARG}"
ENV PATH="${PATH}:${BOWTIE2_PREFIX}/bin"

COPY --from=builder "${BOWTIE2_PREFIX}" "${BOWTIE2_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/bowtie2.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
