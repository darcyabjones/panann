ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

RUN  apt-get update \
  && apt-get install -y \
       libbz2-1.0 \
       libcurl3 \
       libgsl2 \
       libncurses5 \
       libperl5.24 \
       libssl1.1 \
       lzma \
       zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV HTSLIB_PREFIX="/opt/htslib"
ENV BCFTOOLS_PREFIX="/opt/bcftools"
ENV SAMTOOLS_PREFIX="/opt/samtools"

COPY --from="darcyabjones/htslib" "${HTSLIB_PREFIX}" "${HTSLIB_PREFIX}"
COPY --from="darcyabjones/htslib" "${BCFTOOLS_PREFIX}" "${BCFTOOLS_PREFIX}"
COPY --from="darcyabjones/htslib" "${SAMTOOLS_PREFIX}" "${SAMTOOLS_PREFIX}"

ENV PATH="${PATH}:${SAMTOOLS_PREFIX}/bin:${BCFTOOLS_PREFIX}/bin:${HTSLIB_PREFIX}/bin"
ENV CPATH="${HTSLIB_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LD_LIBRARY_PATH}"


ENV VELVET_PREFIX="/opt/velvet"
COPY --from="darcyabjones/velvet" ${VELVET_PREFIX} ${VELVET_PREFIX}

RUN  apt-get update \
  && apt-get install -y libgomp1 zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV PATH="${VELVET_PREFIX}:${PATH}"


ENV JELLYFISH_PREFIX="/opt/jellyfish"

COPY --from="darcyabjones/jellyfish" ${JELLYFISH_PREFIX} ${JELLYFISH_PREFIX}

ENV PATH="${JELLYFISH_PREFIX}:${PATH}"
ENV INCLUDE="${JELLYFISH_DIR}/include:${INCLUDE}"
ENV CPATH="${JELLYFISH_DIR}/include:${CPATH}"
ENV LIBRARY_PATH="${JELLYFISH_DIR}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${JELLYFISH_DIR}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${JELLYFISH_DIR}/lib:${LD_RUN_PATH}"


ENV BWA_PREFIX="/opt/bwa"
COPY --from="darcyabjones/bwa" ${BWA_PREFIX} ${BWA_PREFIX}

RUN  apt-get update \
  && apt-get install -y perl zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV PATH="${BWA_PREFIX}:${PATH}"
