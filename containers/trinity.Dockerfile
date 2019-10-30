ARG IMAGE
ARG JELLYFISH_IMAGE
ARG BOWTIE2_IMAGE
ARG SALMON_IMAGE
ARG HTSLIB_IMAGE

FROM "${JELLYFISH_IMAGE}" as jellyfish_builder
FROM "${BOWTIE2_IMAGE}" as bowtie2_builder
FROM "${SALMON_IMAGE}" as salmon_builder
FROM "${HTSLIB_IMAGE}" as htslib_builder

FROM "${IMAGE}" as trinity_builder


ARG TRINITY_TAG
ARG TRINITY_REPO="https://github.com/trinityrnaseq/trinityrnaseq.git"
ARG TRINITY_PREFIX_ARG="/opt/trinity/${TRINITY_VERSION}"
ENV TRINITY_PREFIX="${TRINITY_PREFIX_ARG}"
ENV TRINITY_HOME="${TRINITY_PREFIX}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       cmake \
       git \
       python \
       rsync \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${TRINITY_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${TRINITY_TAG}" \
  && make \
  && make plugins \
  && mkdir -p "${TRINITY_PREFIX}" \
  && rsync -av --exclude='.*' /tmp/* "${TRINITY_PREFIX}" \
  && add_runtime_dep default-jre-headless python python-numpy perl zlib1g

# CA cert stuff sometime required for git clone https
# make install just calls a python script that does the Rsync bit.

# Runtime requires
    #  Possibly some perl modules

FROM "${IMAGE}"

ARG TRINITY_TAG
ARG TRINITY_PREFIX_ARG="/opt/trinity/${TRINITY_TAG}"
ENV TRINITY_PREFIX="${TRINITY_PREFIX_ARG}"
ENV TRINITY_HOME="${TRINITY_PREFIX}"
LABEL trinity.version="${TRINITY_TAG}"

ENV PATH "${TRINITY_PREFIX}:${PATH}"

COPY --from=trinity_builder "${TRINITY_PREFIX}" "${TRINITY_PREFIX}"
COPY --from=trinity_builder "${APT_REQUIREMENTS_FILE}" /build/apt/trinity.txt


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

# Jellyfish has no runtime deps
COPY --from=jellyfish_builder "${JELLYFISH_PREFIX}" "${JELLYFISH_PREFIX}"


ARG BOWTIE2_TAG
ARG BOWTIE2_PREFIX_ARG="/opt/bowtie2/${BOWTIE2_TAG}"
ENV BOWTIE2_PREFIX="${BOWTIE2_PREFIX_ARG}"
LABEL bowtie2.version="${BOWTIE2_TAG}"

ENV PATH="${PATH}:${BOWTIE2_PREFIX}/bin"

COPY --from=bowtie2_builder "${BOWTIE2_PREFIX}" "${BOWTIE2_PREFIX}"
COPY --from=bowtie2_builder "${APT_REQUIREMENTS_FILE}" /build/apt/bowtie2.txt

ARG SALMON_TAG
ARG SALMON_PREFIX_ARG="/opt/salmon/${SALMON_TAG}"
ENV SALMON_PREFIX="${SALMON_PREFIX_ARG}"
LABEL salmon.version="${SALMON_TAG}"

ENV PATH "${SALMON_PREFIX}/bin:${PATH}"
ENV LD_LIBRARY_PATH "${SALMON_PREFIX}/lib:${LD_LIBRARY_PATH}"

COPY --from=salmon_builder "${SALMON_PREFIX}" "${SALMON_PREFIX}"
COPY --from=salmon_builder "/build/apt-requirements.txt" /build/apt/salmon.txt


ARG HTSLIB_TAG
ARG SAMTOOLS_TAG
ARG HTSLIB_PREFIX_ARG="/opt/htslib/${HTSLIB_TAG}"
ARG SAMTOOLS_PREFIX_ARG="/opt/samtools/${SAMTOOLS_TAG}"
ENV HTSLIB_PREFIX="${HTSLIB_PREFIX_ARG}"
ENV SAMTOOLS_PREFIX="${SAMTOOLS_PREFIX_ARG}"
LABEL htslib.version="${HTSLIB_TAG}"
LABEL samtools.version="${SAMTOOLS_TAG}"

ENV PATH="${SAMTOOLS_PREFIX}/bin:${BCFTOOLS_PREFIX}/bin:${HTSLIB_PREFIX}/bin:${PATH}"
ENV CPATH="${HTSLIB_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LD_LIBRARY_PATH}"

COPY --from=htslib_builder "${HTSLIB_PREFIX}" "${HTSLIB_PREFIX}"
COPY --from=htslib_builder "${SAMTOOLS_PREFIX}" "${SAMTOOLS_PREFIX}"
COPY --from=htslib_builder "${APT_REQUIREMENTS_FILE}" /build/apt/htslib.txt


RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"

WORKDIR /
