ARG IMAGE
ARG BRAKER_IMAGE
ARG AUGUSTUS_IMAGE
ARG HTSLIB_IMAGE
ARG BAMTOOLS_IMAGE
ARG SPALN_IMAGE
ARG GENEMARKES_IMAGE

FROM "${BRAKER_IMAGE}" as braker_builder
FROM "${AUGUSTUS_IMAGE}" as augustus_builder
FROM "${HTSLIB_IMAGE}" as htslib_builder
FROM "${BAMTOOLS_IMAGE}" as bamtools_builder
FROM "${SPALN_IMAGE}" as spaln_builder
FROM "${GENEMARKES_IMAGE}" as genemarkes_builder


FROM "${IMAGE}"

ARG BRAKER_COMMIT
ARG BRAKER_PREFIX_ARG="/opt/braker/${BRAKER_COMMIT}"
ENV BRAKER_PREFIX="${BRAKER_PREFIX_ARG}"

ENV PATH="${BRAKER_PREFIX}/scripts:${PATH}"

COPY --from=braker_builder "${BRAKER_PREFIX}" "${BRAKER_PREFIX}"
COPY --from=braker_builder "${APT_REQUIREMENTS_FILE}" /build/apt/braker.txt


ARG AUGUSTUS_COMMIT
ARG AUGUSTUS_PREFIX_ARG="/opt/augustus/${AUGUSTUS_COMMIT}"
ENV AUGUSTUS_PREFIX="${AUGUSTUS_PREFIX_ARG}"
ENV AUGUSTUS_CONFIG_PATH="${AUGUSTUS_PREFIX}/config"

ENV PATH="${PATH}:${AUGUSTUS_PREFIX}/bin:${AUGUSTUS_PREFIX}/scripts"

COPY --from=augustus_builder "${AUGUSTUS_PREFIX}" "${AUGUSTUS_PREFIX}"
COPY --from=augustus_builder "${APT_REQUIREMENTS_FILE}" /build/apt/augustus.txt


ARG HTSLIB_TAG
ARG SAMTOOLS_TAG
ARG HTSLIB_PREFIX_ARG="/opt/htslib/${HTSLIB_TAG}"
ARG SAMTOOLS_PREFIX_ARG="/opt/samtools/${SAMTOOLS_TAG}"
ENV HTSLIB_PREFIX="${HTSLIB_PREFIX_ARG}"
ENV SAMTOOLS_PREFIX="${SAMTOOLS_PREFIX_ARG}"

ENV PATH="${SAMTOOLS_PREFIX}/bin:${BCFTOOLS_PREFIX}/bin:${HTSLIB_PREFIX}/bin:${PATH}"
ENV CPATH="${HTSLIB_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LD_LIBRARY_PATH}"

COPY --from=htslib_builder "${HTSLIB_PREFIX}" "${HTSLIB_PREFIX}"
COPY --from=htslib_builder "${SAMTOOLS_PREFIX}" "${SAMTOOLS_PREFIX}"
COPY --from=htslib_builder "${APT_REQUIREMENTS_FILE}" /build/apt/htslib.txt


ARG BAMTOOLS_TAG="v2.5.1"
ARG BAMTOOLS_PREFIX_ARG="/opt/bamtools/${BAMTOOLS_TAG}"
ENV BAMTOOLS_PREFIX="${BAMTOOLS_PREFIX_ARG}"

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
ENV PATH="${SPALN_PREFIX}/bin:${SPALN_PREFIX}/perl:${PATH}"

COPY --from=spaln_builder "${SPALN_PREFIX}" "${SPALN_PREFIX}"
COPY --from=spaln_builder "${APT_REQUIREMENTS_FILE}" /build/apt/spaln.txt


ARG GENEMARKES_VERSION
ARG GENEMARKES_PREFIX_ARG="/opt/genemarkes/${GENEMARKES_VERSION}"
ENV GENEMARKES_PREFIX="${GENEMARKES_PREFIX_ARG}"

COPY --from=genemarkes_builder "${GENEMARKES_PREFIX}" "${GENEMARKES_PREFIX}"
COPY --from=genemarkes_builder "${GENEMARKES_PREFIX}/gm_key" /root/.gm_key
COPY --from=genemarkes_builder "${APT_REQUIREMENTS_FILE}" /build/apt/genemarkes.txt

ENV PATH="${GENEMARKES_PREFIX}:${PATH}"

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
