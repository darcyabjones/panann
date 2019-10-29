ARG IMAGE
ARG AUGUSTUS_IMAGE
ARG HTSLIB_IMAGE
ARG BAMTOOLS_IMAGE
ARG SPALN_IMAGE

FROM "${AUGUSTUS_IMAGE}" as augustus_builder
FROM "${HTSLIB_IMAGE}" as htslib_builder
FROM "${BAMTOOLS_IMAGE}" as bamtools_builder
FROM "${SPALN_IMAGE}" as spaln_builder


FROM "${IMAGE}" as braker_builder

ARG BRAKER_COMMIT="e117150b8ad66ecf7cd5828c7f7fe476a4a8c191"
ARG BRAKER_REPO="https://github.com/Gaius-Augustus/BRAKER.git"
ARG BRAKER_PREFIX_ARG="/opt/braker/${BRAKER_COMMIT}"
ENV BRAKER_PREFIX="${BRAKER_PREFIX_ARG}"

# CPAN package is in debian unstable still
# liblogger-simple-perl

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       ca-certificates \
       git \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${BRAKER_REPO}" "${BRAKER_PREFIX}" \
  && cd "${BRAKER_PREFIX}" \
  && git checkout "${BRAKER_COMMIT}" \
  && rm -rf -- docs example \
  && add_runtime_dep \
       ncbi-blast+ \
       libhash-merge-perl \
       libmodule-load-conditional-perl \
       libparallel-forkmanager-perl \
       libscalar-util-numeric-perl \
       libyaml-perl \
       perl \
       python3 \
       python3-biopython


FROM "${IMAGE}"

ARG BRAKER_COMMIT="e117150b8ad66ecf7cd5828c7f7fe476a4a8c191"
ARG BRAKER_PREFIX_ARG="/opt/braker/${BRAKER_COMMIT}"
ENV BRAKER_PREFIX="${BRAKER_PREFIX_ARG}"
LABEL braker.version="${BRAKER_COMMIT}"

ENV PATH="${BRAKER_PREFIX}/scripts:${PATH}"

COPY --from=braker_builder "${BRAKER_PREFIX}" "${BRAKER_PREFIX}"
COPY --from=braker_builder "${APT_REQUIREMENTS_FILE}" /build/apt/braker.txt

ARG AUGUSTUS_COMMIT
ARG AUGUSTUS_PREFIX_ARG="/opt/augustus/${AUGUSTUS_COMMIT}"
ENV AUGUSTUS_PREFIX="${AUGUSTUS_PREFIX_ARG}"
ENV AUGUSTUS_CONFIG_PATH="${AUGUSTUS_PREFIX}/config"
LABEL augustus.version="${AUGUSTUS_COMMIT}"

ENV PATH="${PATH}:${AUGUSTUS_PREFIX}/bin:${AUGUSTUS_PREFIX}/scripts"

COPY --from=augustus_builder "${AUGUSTUS_PREFIX}" "${AUGUSTUS_PREFIX}"
COPY --from=augustus_builder "${APT_REQUIREMENTS_FILE}" /build/apt/augustus.txt


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


ARG BAMTOOLS_TAG
ARG BAMTOOLS_PREFIX_ARG="/opt/bamtools/${BAMTOOLS_TAG}"
ENV BAMTOOLS_PREFIX="${BAMTOOLS_PREFIX_ARG}"
LABEL bamtools.version="${BAMTOOLS_TAG}"

ENV PATH="${BAMTOOLS_PREFIX}/bin:${PATH}"
ENV LIBRARY_PATH="${LD_LIBRARY_PATH}:${BAMTOOLS_PREFIX}/lib"
ENV CPATH="${CPATH}:${BAMTOOLS_PREFIX}/include"

COPY --from=bamtools_builder "${BAMTOOLS_PREFIX}" "${BAMTOOLS_PREFIX}"
COPY --from=bamtools_builder "${APT_REQUIREMENTS_FILE}" /build/apt/bamtools.txt


ARG SPALN_TAG
ARG SPALN_PREFIX_ARG="/opt/spaln/${SPALN_TAG}"
ENV SPALN_PREFIX="${SPALN_PREFIX_ARG}"
ENV ALN_TAB="${SPALN_PREFIX}/table"
ENV ALN_DBS="${SPALN_PREFIX}/seqdb"
LABEL spaln.version="${SPALN_TAG}"

ENV PATH="${SPALN_PREFIX}/bin:${SPALN_PREFIX}/perl:${PATH}"

COPY --from=spaln_builder "${SPALN_PREFIX}" "${SPALN_PREFIX}"
COPY --from=spaln_builder "${APT_REQUIREMENTS_FILE}" /build/apt/spaln.txt

# CPAN package is in debian unstable still
# liblogger-simple-perl

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && apt-get install -y --no-install-recommends \
       cpanminus \
       make \
       tar \
  && rm -rf /var/lib/apt/lists/* \
  && cpanm Logger::Simple \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
