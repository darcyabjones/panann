ARG IMAGE
ARG GMAP_IMAGE
ARG FASTA_IMAGE
ARG HTSLIB_IMAGE

FROM "${GMAP_IMAGE}" as gmap_builder
FROM "${FASTA_IMAGE}" as fasta_builder
FROM "${HTSLIB_IMAGE}" as htslib_builder

FROM "${IMAGE}" as builder

ARG PASA_TAG
ARG PASA_REPO="https://github.com/PASApipeline/PASApipeline.git"
ARG PASA_PREFIX_ARG="/opt/pasa/${PASA_TAG}"
ENV PASA_PREFIX="${PASA_PREFIX_ARG}"
ENV PASA_HOME="${PASA_PREFIX}"

WORKDIR /tmp/pasa
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       git \
       wget \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${PASA_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${PASA_TAG}" \
  && git submodule init \
  && git submodule update \
  && make \
  && mkdir -p "${PASA_PREFIX}" \
  && cp -r bin "${PASA_PREFIX}/bin" \
  && cp -r scripts "${PASA_PREFIX}/scripts" \
  && cp -r schema "${PASA_PREFIX}/schema" \
  && cp -r pasa_conf "${PASA_PREFIX}/pasa_conf" \
  && cp -r misc_utilities "${PASA_PREFIX}/misc_utilities" \
  && rm -rf -- "${PASA_PREFIX}/misc_utilities/deprecated" \
  && cp -r PerlLib "${PASA_PREFIX}/PerlLib" \
  && cp -r PyLib "${PASA_PREFIX}/PyLib" \
  && mkdir -p "${PASA_PREFIX}/pasa-plugins" \
  && cp -r pasa-plugins/transdecoder "${PASA_PREFIX}/pasa-plugins/transdecoder" \
  && rm -r "${PASA_PREFIX}/pasa-plugins/transdecoder/sample_data" \
  && cp -r SAMPLE_HOOKS "${PASA_PREFIX}/SAMPLE_HOOKS" \
  && cp Launch_PASA_pipeline.pl "${PASA_PREFIX}/bin" \
  && cp Launch_PASA_pipeline.pl "${PASA_PREFIX}" \
  && wget -O "${PASA_PREFIX}/bin/blat" http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat \
  && chmod 755 "${PASA_PREFIX}/bin/blat" \
  && add_runtime_dep liburi-perl libdbd-sqlite3-perl perl python sqlite3 zlib1g

# Uncomment this to keep the tests around.
# && cp -r sample_data "${PASA_PREFIX}/sample_data" \

FROM "${IMAGE}"

ARG PASA_TAG
ARG PASA_PREFIX_ARG="/opt/pasa/${PASA_TAG}"
ENV PASA_PREFIX="${PASA_PREFIX_ARG}"
ENV PASA_HOME="${PASA_PREFIX}"

ENV PATH "${PASA_PREFIX}/bin:${PASA_PREFIX}/misc_utilities:${PASA_PREFIX}/scripts:${PATH}"

COPY --from=builder "${PASA_PREFIX}" "${PASA_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/pasa.txt

ARG GMAP_VERSION
ARG GMAP_PREFIX_ARG="/opt/gmap/${GMAP_VERSION}"
ENV GMAP_PREFIX="${GMAP_PREFIX_ARG}"

ENV PATH "${GMAP_PREFIX}/bin:${PATH}"

COPY --from=gmap_builder "${GMAP_PREFIX}" "${GMAP_PREFIX}"
COPY --from=gmap_builder "${APT_REQUIREMENTS_FILE}" /build/apt/gmap.txt

ARG FASTA_VERSION
ARG FASTA_PREFIX_ARG="/opt/fasta/${FASTA_VERSION}"
ENV FASTA_PREFIX="${FASTA_PREFIX_ARG}"

ENV PATH="${FASTA_PREFIX}/bin:${FASTA_PREFIX}/psisearch2:${PATH}"

COPY --from=fasta_builder "${FASTA_PREFIX}" "${FASTA_PREFIX}"
COPY --from=fasta_builder "${APT_REQUIREMENTS_FILE}" /build/apt/fasta.txt

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


RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
