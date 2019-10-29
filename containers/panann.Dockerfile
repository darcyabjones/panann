ARG IMAGE
ARG AEGEAN_IMAGE
ARG AUGUSTUS_IMAGE
ARG BAMTOOLS_IMAGE
ARG BEDTOOLS_IMAGE
ARG BOWTIE2_IMAGE
ARG BRAKER_IMAGE
ARG BUSCO_IMAGE
ARG CODINGQUARRY_IMAGE
ARG DEEPSIG_IMAGE
ARG EVM_IMAGE
ARG EXONERATE_IMAGE
ARG FASTA_IMAGE
ARG GEMOMA_IMAGE
ARG GENOMETOOLS_IMAGE
ARG GFFPAL_IMAGE
ARG GMAP_IMAGE
ARG HTSLIB_IMAGE
ARG JELLYFISH_IMAGE
ARG MMSEQS_IMAGE
ARG PASA_IMAGE
ARG SALMON_IMAGE
ARG SEQRENAMER_IMAGE
ARG SPALN_IMAGE
ARG STAR_IMAGE
ARG STRINGTIE_IMAGE
ARG TRINITY_IMAGE

FROM "${AEGEAN_IMAGE}" as aegean_builder
FROM "${AUGUSTUS_IMAGE}" as augustus_builder
FROM "${BAMTOOLS_IMAGE}" as bamtools_builder
FROM "${BEDTOOLS_IMAGE}" as bedtools_builder
FROM "${BOWTIE2_IMAGE}" as bowtie2_builder
FROM "${BRAKER_IMAGE}" as braker_builder
FROM "${BUSCO_IMAGE}" as busco_builder
FROM "${CODINGQUARRY_IMAGE}" as codingquarry_builder
FROM "${DEEPSIG_IMAGE}" as deepsig_builder
FROM "${EVM_IMAGE}" as evm_builder
FROM "${EXONERATE_IMAGE}" as exonerate_builder
FROM "${FASTA_IMAGE}" as fasta_builder
FROM "${GEMOMA_IMAGE}" as gemoma_builder
FROM "${GENOMETOOLS_IMAGE}" as genometools_builder
FROM "${GFFPAL_IMAGE}" as gffpal_builder
FROM "${GMAP_IMAGE}" as gmap_builder
FROM "${HTSLIB_IMAGE}" as htslib_builder
FROM "${JELLYFISH_IMAGE}" as jellyfish_builder
FROM "${MMSEQS_IMAGE}" as mmseqs_builder
FROM "${PASA_IMAGE}" as pasa_builder
FROM "${SALMON_IMAGE}" as salmon_builder
FROM "${SEQRENAMER_IMAGE}" as seqrenamer_builder
FROM "${SPALN_IMAGE}" as spaln_builder
FROM "${STAR_IMAGE}" as star_builder
FROM "${STRINGTIE_IMAGE}" as stringtie_builder
FROM "${TRINITY_IMAGE}" as trinity_builder

FROM "${IMAGE}"


ARG AEGEAN_VERSION
ARG AEGEAN_PREFIX_ARG="/opt/aegean/${AEGEAN_VERSION}"
ENV AEGEAN_PREFIX="${AEGEAN_PREFIX_ARG}"
LABEL aegean.version="${AEGEAN_VERSION}"

ENV PATH="${AEGEAN_PREFIX}/bin:${PATH}"
ENV INCLUDE="${AEGEAN_PREFIX}/include:${INCLUDE}"
ENV CPATH="${AEGEAN_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${AEGEAN_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${AEGEAN_PREFIX}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${AEGEAN_PREFIX}/lib:${LD_RUN_PATH}"

COPY --from=aegean_builder "${AEGEAN_PREFIX}" "${AEGEAN_PREFIX}"
COPY --from=aegean_builder "${APT_REQUIREMENTS_FILE}" /build/apt/aegean.txt


ARG AUGUSTUS_COMMIT
ARG AUGUSTUS_PREFIX_ARG="/opt/augustus/${AUGUSTUS_COMMIT}"
ENV AUGUSTUS_PREFIX="${AUGUSTUS_PREFIX_ARG}"

ENV PATH="${AUGUSTUS_PREFIX}/bin:${AUGUSTUS_PREFIX}/scripts:${PATH}"
ENV AUGUSTUS_CONFIG_PATH="${AUGUSTUS_PREFIX}/config"

COPY --from=augustus_builder "${AUGUSTUS_PREFIX}" "${AUGUSTUS_PREFIX}"
COPY --from=augustus_builder "${APT_REQUIREMENTS_FILE}" /build/apt/augustus.txt


ARG BAMTOOLS_TAG
ARG BAMTOOLS_PREFIX_ARG="/opt/bamtools/${BAMTOOLS_TAG}"
ENV BAMTOOLS_PREFIX="${BAMTOOLS_PREFIX_ARG}"
LABEL bamtools.version="${BAMTOOLS_TAG}"

ENV PATH="${BAMTOOLS_PREFIX}/bin:${PATH}"
ENV LIBRARY_PATH="${LD_LIBRARY_PATH}:${BAMTOOLS_PREFIX}/lib"
ENV CPATH="${CPATH}:${BAMTOOLS_PREFIX}/include"

COPY --from=bamtools_builder "${BAMTOOLS_PREFIX}" "${BAMTOOLS_PREFIX}"
COPY --from=bamtools_builder "${APT_REQUIREMENTS_FILE}" /build/apt/bamtools.txt


ARG BEDTOOLS_VERSION
ARG BEDTOOLS_PREFIX_ARG="/opt/bedtools/${BEDTOOLS_VERSION}"
ENV BEDTOOLS_PREFIX="${BEDTOOLS_PREFIX_ARG}"
LABEL bedtools.version="${BEDTOOLS_VERSION}"

ENV PATH "${BEDTOOLS_PREFIX}/bin:${PATH}"

COPY --from=bedtools_builder "${BEDTOOLS_PREFIX}" "${BEDTOOLS_PREFIX}"
COPY --from=bedtools_builder "${APT_REQUIREMENTS_FILE}" /build/apt/bedtools.txt


ARG BOWTIE2_TAG
ARG BOWTIE2_PREFIX_ARG="/opt/bowtie2/${BOWTIE2_TAG}"
ENV BOWTIE2_PREFIX="${BOWTIE2_PREFIX_ARG}"
ENV PATH="${PATH}:${BOWTIE2_PREFIX}/bin"

COPY --from=bowtie2_builder "${BOWTIE2_PREFIX}" "${BOWTIE2_PREFIX}"
COPY --from=bowtie2_builder "${APT_REQUIREMENTS_FILE}" /build/apt/bowtie2.txt


ARG BRAKER_COMMIT="e117150b8ad66ecf7cd5828c7f7fe476a4a8c191"
ARG BRAKER_PREFIX_ARG="/opt/braker/${BRAKER_COMMIT}"
ENV BRAKER_PREFIX="${BRAKER_PREFIX_ARG}"
LABEL braker.version="${BRAKER_COMMIT}"

ENV PATH="${BRAKER_PREFIX}/scripts:${PATH}"

COPY --from=braker_builder "${BRAKER_PREFIX}" "${BRAKER_PREFIX}"
COPY --from=braker_builder "${APT_REQUIREMENTS_FILE}" /build/apt/braker.txt


ARG BUSCO_COMMIT
ARG BUSCO_PREFIX_ARG="/opt/busco/${BUSCO_COMMIT}"
ENV BUSCO_PREFIX="${BUSCO_PREFIX_ARG}"

ENV PATH="${BUSCO_PREFIX}/scripts:${PATH}"
ENV PYTHONPATH="${PYTHONPATH}:${BUSCO_PREFIX}/lib/python3.5/site-packages"

COPY --from=busco_builder "${BUSCO_PREFIX}" "${BUSCO_PREFIX}"
COPY --from=busco_builder "${APT_REQUIREMENTS_FILE}" /build/apt/busco.txt


ARG CODINGQUARRY_VERSION
ARG CODINGQUARRY_URL="https://downloads.sourceforge.net/project/codingquarry/CodingQuarry_v2.0.tar.gz"
ARG CODINGQUARRY_PREFIX_ARG="/opt/codingquarry/${CODINGQUARRY_VERSION}"
ENV CODINGQUARRY_PREFIX="${CODINGQUARRY_PREFIX_ARG}"
ENV QUARRY_PATH="${CODINGQUARRY_PREFIX}/QuarryFiles"

ENV PATH "${CODINGQUARRY_PREFIX}/bin:${PATH}"

COPY --from=codingquarry_builder "${CODINGQUARRY_PREFIX}" "${CODINGQUARRY_PREFIX}"
COPY --from=codingquarry_builder "${APT_REQUIREMENTS_FILE}" /build/apt/codingquarry.txt


ARG DEEPSIG_COMMIT="69e01cb"
ARG DEEPSIG_PREFIX_ARG="/opt/deepsig/${DEEPSIG_COMMIT}"
ENV DEEPSIG_PREFIX="${DEEPSIG_PREFIX_ARG}"
ENV DEEPSIG_ROOT="${DEEPSIG_PREFIX}"
LABEL deepsig.version="${DEEPSIG_COMMIT}"

ARG TENSORFLOW_VERSION
ARG TENSORFLOW_PREFIX_ARG="/opt/tensorflow/${TENSORFLOW_VERSION}"
ENV TENSORFLOW_PREFIX="${TENSORFLOW_PREFIX_ARG}"
ENV TF_CPP_MIN_LOG_LEVEL=3
LABEL tensorflow.version="TENSORFLOW_VERSION"

ARG KERAS_VERSION
ARG KERAS_PREFIX_ARG="/opt/keras/${KERAS_VERSION}"
ENV KERAS_PREFIX="${KERAS_PREFIX_ARG}"
LABEL keras.version="KERAS_VERSION"

ENV PATH="${DEEPSIG_PREFIX}:${PATH}"
ENV PATH="${TENSORFLOW_PREFIX}:${PATH}"
ENV PYTHONPATH="${TENSORFLOW_PREFIX}/lib/python2.7/site-packages:${PYTHONPATH:-}"
ENV PYTHONPATH="${KERAS_PREFIX}/lib/python2.7/site-packages:${PYTHONPATH:-}"

COPY --from=deepsig_builder "${DEEPSIG_PREFIX}" "${DEEPSIG_PREFIX}"
COPY --from=deepsig_builder "${TENSORFLOW_PREFIX}" "${TENSORFLOW_PREFIX}"
COPY --from=deepsig_builder "${KERAS_PREFIX}" "${KERAS_PREFIX}"
COPY --from=deepsig_builder "${APT_REQUIREMENTS_FILE}" /build/apt/deepsig.txt


ARG EVM_COMMIT
ARG EVM_PREFIX_ARG
ENV EVM_PREFIX="${EVM_PREFIX_ARG}"

LABEL evm.version="${EVM_COMMIT}"

ENV PATH "${EVM_PREFIX}:${EVM_PREFIX}/EvmUtils:${EVM_PREFIX}/EvmUtils/misc:${PATH}"

COPY --from=evm_builder "${EVM_PREFIX}" "${EVM_PREFIX}"
COPY --from=evm_builder "${APT_REQUIREMENTS_FILE}" /build/apt/evm.txt


ARG EXONERATE_VERSION
ARG EXONERATE_PREFIX_ARG="/opt/exonerate/${EXONERATE_VERSION}"
ENV EXONERATE_PREFIX="${EXONERATE_PREFIX_ARG}"

LABEL exonerate.version="${EXONERATE_VERSION}"

ENV PATH "${EXONERATE_PREFIX}/bin:${PATH}"

COPY --from=exonerate_builder "${EXONERATE_PREFIX}" "${EXONERATE_PREFIX}"
COPY --from=exonerate_builder "${APT_REQUIREMENTS_FILE}" /build/apt/exonerate.txt


ARG FASTA_VERSION
ARG FASTA_PREFIX_ARG="/opt/fasta/${FASTA_VERSION}"
ENV FASTA_PREFIX="${FASTA_PREFIX_ARG}"
LABEL fasta.version="${FASTA_VERSION}"

ENV PATH="${FASTA_PREFIX}/bin:${FASTA_PREFIX}/psisearch2:${PATH}"

COPY --from=fasta_builder "${FASTA_PREFIX}" "${FASTA_PREFIX}"
COPY --from=fasta_builder "${APT_REQUIREMENTS_FILE}" /build/apt/fasta.txt


ARG GEMOMA_VERSION
ARG GEMOMA_PREFIX_ARG="/opt/gemoma/${GEMOMA_VERSION}"
ENV GEMOMA_PREFIX="${GEMOMA_PREFIX_ARG}"
ENV GEMOMA_JAR="${GEMOMA_PREFIX}/bin/GeMoMa-${GEMOMA_VERSION}.jar"
LABEL gemoma.version="${GEMOMA_VERSION}"

ENV PATH="${GEMOMA_PREFIX}/bin:${PATH}"

COPY --from=gemoma_builder "${GEMOMA_PREFIX}" "${GEMOMA_PREFIX}"
COPY --from=gemoma_builder "${APT_REQUIREMENTS_FILE}" /build/apt/gemoma.txt


ARG GENOMETOOLS_VERSION
ARG GENOMETOOLS_PREFIX_ARG="/opt/genometools/${GENOMETOOLS_VERSION}"
ENV GENOMETOOLS_PREFIX="${GENOMETOOLS_PREFIX_ARG}"

ENV PATH="${GENOMETOOLS_PREFIX}/bin:${PATH}"
ENV INCLUDE="${GENOMETOOLS_PREFIX}/include:${INCLUDE}"
ENV CPATH="${GENOMETOOLS_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${GENOMETOOLS_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${GENOMETOOLS_PREFIX}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${GENOMETOOLS_PREFIX}/lib:${LD_RUN_PATH}"

COPY --from=genometools_builder "${GENOMETOOLS_PREFIX}" "${GENOMETOOLS_PREFIX}"
COPY --from=genometools_builder "${APT_REQUIREMENTS_FILE}" /build/apt/genometools.txt


ARG GFFPAL_TAG
ARG GFFPAL_PREFIX_ARG="/opt/gffpal/${GFFPAL_TAG}"
ENV GFFPAL_PREFIX="${GFFPAL_PREFIX_ARG}"
LABEL gffpal.version="${GFFPAL_TAG}"

ENV PATH "${GFFPAL_PREFIX}/bin:${PATH}"
ENV PYTHONPATH "${GFFPAL_PREFIX}/lib/python3.7/site-packages:${PYTHONPATH}"

COPY --from=gffpal_builder "${GFFPAL_PREFIX}" "${GFFPAL_PREFIX}"
COPY --from=gffpal_builder "${APT_REQUIREMENTS_FILE}" /build/apt/gffpal.txt


ARG GMAP_VERSION
ARG GMAP_PREFIX_ARG="/opt/gmap/${GMAP_VERSION}"
ENV GMAP_PREFIX="${GMAP_PREFIX_ARG}"

ENV PATH "${GMAP_PREFIX}/bin:${PATH}"

COPY --from=gmap_builder "${GMAP_PREFIX}" "${GMAP_PREFIX}"
COPY --from=gmap_builder "${APT_REQUIREMENTS_FILE}" /build/apt/gmap.txt


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


ARG JELLYFISH_VERSION
ARG JELLYFISH_PREFIX_ARG="/opt/jellyfish/${JELLYFISH_VERSION}"
ENV JELLYFISH_PREFIX="${JELLYFISH_PREFIX_ARG}"

ENV PATH="${JELLYFISH_PREFIX}/bin:${PATH}"
ENV INCLUDE="${JELLYFISH_DIR}/include:${INCLUDE}"
ENV CPATH="${JELLYFISH_DIR}/include:${CPATH}"
ENV LIBRARY_PATH="${JELLYFISH_DIR}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${JELLYFISH_DIR}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${JELLYFISH_DIR}/lib:${LD_RUN_PATH}"

# Jellyfish has no runtime deps
COPY --from=jellyfish_builder "${JELLYFISH_PREFIX}" "${JELLYFISH_PREFIX}"


ARG MMSEQS_TAG
ARG MMSEQS_PREFIX_ARG
ENV MMSEQS_PREFIX="${MMSEQS_PREFIX_ARG}"
LABEL mmseqs.version="${MMSEQS_TAG}"

ENV PATH="${MMSEQS_PREFIX}/bin:${PATH}"

COPY --from=mmseqs_builder "${MMSEQS_PREFIX}" "${MMSEQS_PREFIX}"
COPY --from=mmseqs_builder "${APT_REQUIREMENTS_FILE}" /build/apt/mmseqs.txt


ARG PASA_TAG
ARG PASA_PREFIX_ARG
ENV PASA_PREFIX="${PASA_PREFIX_ARG}"
ENV PASA_HOME="${PASA_PREFIX}"
LABEL pasa.version="${PASA_TAG}"

ENV PATH "${PASA_PREFIX}:${PASA_PREFIX}/bin:${PASA_PREFIX}/misc_utilities:${PASA_PREFIX}/scripts:${PATH}"

COPY --from=pasa_builder "${PASA_PREFIX}" "${PASA_PREFIX}"
COPY --from=pasa_builder "${APT_REQUIREMENTS_FILE}" /build/apt/pasa.txt


ARG SALMON_TAG
ARG SALMON_PREFIX_ARG="/opt/salmon/${SALMON_TAG}"
ENV SALMON_PREFIX="${SALMON_PREFIX_ARG}"

ENV PATH "${SALMON_PREFIX}/bin:${PATH}"
ENV LD_LIBRARY_PATH "${SALMON_PREFIX}/lib:${LD_LIBRARY_PATH}"

COPY --from=salmon_builder "${SALMON_PREFIX}" "${SALMON_PREFIX}"
COPY --from=salmon_builder "/build/apt-requirements.txt" /build/apt/salmon.txt


ARG SEQRENAMER_TAG
ARG SEQRENAMER_PREFIX_ARG="/opt/seqrenamer/${SEQRENAMER_TAG}"
ENV SEQRENAMER_PREFIX="${SEQRENAMER_PREFIX_ARG}"
LABEL seqrenamer.version="${SEQRENAMER_TAG}"

ENV PATH "${SEQRENAMER_PREFIX}/bin:${PATH}"
ENV PYTHONPATH "${SEQRENAMER_PREFIX}/lib/python3.7/site-packages:${PYTHONPATH}"

COPY --from=seqrenamer_builder "${SEQRENAMER_PREFIX}" "${SEQRENAMER_PREFIX}"
COPY --from=seqrenamer_builder "${APT_REQUIREMENTS_FILE}" /build/apt/seqrenamer.txt


ARG SPALN_TAG
ARG SPALN_PREFIX_ARG="/opt/spaln/${SPALN_TAG}"
ENV SPALN_PREFIX="${SPALN_PREFIX_ARG}"
ENV ALN_TAB="${SPALN_PREFIX}/table"
ENV ALN_DBS="${SPALN_PREFIX}/seqdb"
ENV PATH="${SPALN_PREFIX}/bin:${SPALN_PREFIX}/perl:${PATH}"

COPY --from=spaln_builder "${SPALN_PREFIX}" "${SPALN_PREFIX}"
COPY --from=spaln_builder "${APT_REQUIREMENTS_FILE}" /build/apt/spaln.txt


ARG STAR_VERSION
ARG STAR_PREFIX_ARG="/opt/star/${STAR_VERSION}"
ENV STAR_PREFIX="${STAR_PREFIX_ARG}"
ENV PATH="${STAR_PREFIX}/bin:${PATH}"

COPY --from=star_builder "${STAR_PREFIX}" "${STAR_PREFIX}"
COPY --from=star_builder "${APT_REQUIREMENTS_FILE}" /build/apt/star.txt


ARG STRINGTIE_VERSION
ARG STRINGTIE_PREFIX_ARG="/opt/stringtie/${STRINGTIE_VERSION}"
ENV STRINGTIE_PREFIX="${STRINGTIE_PREFIX_ARG}"
ENV PATH "${STRINGTIE_PREFIX}/bin:${PATH}"

COPY --from=stringtie_builder "${STRINGTIE_PREFIX}" "${STRINGTIE_PREFIX}"
COPY --from=stringtie_builder "${APT_REQUIREMENTS_FILE}" /build/apt/stringtie.txt


ARG TRINITY_TAG
ARG TRINITY_PREFIX_ARG="/opt/trinity/${TRINITY_TAG}"
ENV TRINITY_PREFIX="${TRINITY_PREFIX_ARG}"
ENV TRINITY_HOME="${TRINITY_PREFIX}"

ENV PATH "${TRINITY_PREFIX}:${PATH}"

COPY --from=trinity_builder "${TRINITY_PREFIX}" "${TRINITY_PREFIX}"
COPY --from=trinity_builder "${APT_REQUIREMENTS_FILE}" /build/apt/trinity.txt

ENV PYTHONPATH=""

# liblogger-simple-perl is still in unstable.
# needed for braker.

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get upgrade \
  && dpkg --configure -a \
  && apt-get install --fix-broken -y --no-install-recommends python3-minimal \
  && apt_install_from_file /build/apt/*.txt \
  && apt-get install -y --no-install-recommends \
       cpanminus \
       make \
       tar \
  && rm -rf /var/lib/apt/lists/* \
  && apt-get autoremove \
  && cpanm Logger::Simple \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"

WORKDIR /
