ARG IMAGE
ARG HTSLIB_IMAGE

FROM "${HTSLIB_IMAGE}" as htslib_builder

FROM "${IMAGE}" as star_builder

ARG STAR_VERSION
ARG STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz"
ARG STAR_PREFIX_ARG="/opt/star/${STAR_VERSION}"
ENV STAR_PREFIX="${STAR_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       build-essential \
       ca-certificates \
       wget \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O star.tar.gz "${STAR_URL}" \
  && tar -zxf star.tar.gz \
  && cd STAR*/source \
  && make LDFLAGSextra=-flto CXXFLAGSextra="-flto -march=x86-64 -mtune=native" STAR STARlong \
  && mkdir -p "${STAR_PREFIX}/bin" \
  && mv ../bin/Linux_x86_64/STAR* "${STAR_PREFIX}/bin" \
  && add_runtime_dep libgomp1 zlib1g

# See if we need to provide build flags as a build time-option?
# STAR already seems to require x86-64 instructions, so i think current settings are safe.
# LDFLAGSextra=-flto CXXFLAGSextra="-flto -march=x86-64 -mtune=native"


FROM "${IMAGE}"

ARG STAR_VERSION
ARG STAR_PREFIX_ARG="/opt/star/${STAR_VERSION}"
ENV STAR_PREFIX="${STAR_PREFIX_ARG}"
ENV PATH="${STAR_PREFIX}/bin:${PATH}"
LABEL star.version="${STAR_VERSION}"

COPY --from=star_builder "${STAR_PREFIX}" "${STAR_PREFIX}"
COPY --from=star_builder "${APT_REQUIREMENTS_FILE}" /build/apt/star.txt

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
