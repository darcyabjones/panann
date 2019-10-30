ARG IMAGE
ARG HTSLIB_IMAGE

FROM "${HTSLIB_IMAGE}" as htslib_builder

FROM "${IMAGE}" as stringtie_builder

ARG STRINGTIE_VERSION
ARG STRINGTIE_URL="http://ccb.jhu.edu/software/stringtie/dl/stringtie-${STRINGTIE_VERSION}.tar.gz"
ARG STRINGTIE_PREFIX_ARG="/opt/stringtie/${STRINGTIE_VERSION}"
ENV STRINGTIE_PREFIX="${STRINGTIE_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       wget \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && wget -O stringtie.tar.gz "${STRINGTIE_URL}" \
  && tar -zxf stringtie.tar.gz \
  && cd stringtie* \
  && make release \
  && mkdir -p "${STRINGTIE_PREFIX}/bin" \
  && cp -r stringtie "${STRINGTIE_PREFIX}/bin" \
  && add_runtime_dep zlib1g


FROM "${IMAGE}"

ARG STRINGTIE_VERSION
ARG STRINGTIE_PREFIX_ARG="/opt/stringtie/${STRINGTIE_VERSION}"
ENV STRINGTIE_PREFIX="${STRINGTIE_PREFIX_ARG}"
LABEL stringtie.version="${STRINGTIE_VERSION}"

ENV PATH "${STRINGTIE_PREFIX}/bin:${PATH}"

COPY --from=stringtie_builder "${STRINGTIE_PREFIX}" "${STRINGTIE_PREFIX}"
COPY --from=stringtie_builder "${APT_REQUIREMENTS_FILE}" /build/apt/stringtie.txt

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
