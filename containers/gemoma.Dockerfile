ARG IMAGE
ARG HTSLIB_IMAGE
ARG MMSEQS_IMAGE

FROM "${HTSLIB_IMAGE}" as htslib_builder
FROM "${MMSEQS_IMAGE}" as mmseqs_builder

FROM "${IMAGE}" as gemoma_builder

ARG GEMOMA_VERSION
ARG GEMOMA_URL="http://www.jstacs.de/downloads/GeMoMa-1.6.1.zip"
ARG GEMOMA_PREFIX_ARG="/opt/gemoma/${GEMOMA_VERSION}"
ENV GEMOMA_PREFIX="${GEMOMA_PREFIX_ARG}"
ENV GEMOMA_JAR="${GEMOMA_PREFIX}/bin/GeMoMa-${GEMOMA_VERSION}.jar"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       unzip \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget "${GEMOMA_URL}" \
  && unzip "GeMoMa-${GEMOMA_VERSION}.zip" \
  && rm *.zip \
  && mkdir -p "${GEMOMA_PREFIX}/bin" \
  && cp "GeMoMa-${GEMOMA_VERSION}.jar" "${GEMOMA_PREFIX}/bin" \
  && ln -sf "${GEMOMA_PREFIX}/bin/GeMoMa-${GEMOMA_VERSION}.jar" "${GEMOMA_PREFIX}/bin/GeMoMa.jar" \
  && add_runtime_dep default-jre-headless


FROM "${IMAGE}"

ARG GEMOMA_VERSION
ARG GEMOMA_PREFIX_ARG="/opt/gemoma/${GEMOMA_VERSION}"
ENV GEMOMA_PREFIX="${GEMOMA_PREFIX_ARG}"
ENV GEMOMA_JAR="${GEMOMA_PREFIX}/bin/GeMoMa-${GEMOMA_VERSION}.jar"
LABEL gemoma.version="${GEMOMA_VERSION}"

ENV PATH="${GEMOMA_PREFIX}/bin:${PATH}"

COPY --from=gemoma_builder "${GEMOMA_PREFIX}" "${GEMOMA_PREFIX}"
COPY --from=gemoma_builder "${APT_REQUIREMENTS_FILE}" /build/apt/gemoma.txt


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

ARG MMSEQS_TAG
ARG MMSEQS_PREFIX_ARG
ENV MMSEQS_PREFIX="${MMSEQS_PREFIX_ARG}"
LABEL mmseqs.version="${MMSEQS_TAG}"

ENV PATH="${MMSEQS_PREFIX}/bin:${PATH}"

COPY --from=mmseqs_builder "${MMSEQS_PREFIX}" "${MMSEQS_PREFIX}"
COPY --from=mmseqs_builder "${APT_REQUIREMENTS_FILE}" /build/apt/mmseqs.txt

# Calling the jar file is necessary to create the config xml file.
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && java -jar "${GEMOMA_JAR}" CLI \
  && rm -rf -- GeMoMa_temp \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"


WORKDIR /
