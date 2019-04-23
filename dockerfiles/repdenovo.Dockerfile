ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

ENV BAMTOOLS_PREFIX="/opt/bamtools"
COPY --from="darcyabjones/bamtools" ${BAMTOOLS_PREFIX} ${BAMTOOLS_PREFIX}

ENV PATH="${BAMTOOLS_PREFIX}/bin:${PATH}"
ENV LIBRARY_PATH="${LD_LIBRARY_PATH}:${BAMTOOLS_PREFIX}/lib"
ENV CPATH="${CPATH}:${BAMTOOLS_PREFIX}/include"

ENV REPDENOVO_PREFIX="/opt/repdenovo"
ARG REPDENOVO_VERSION="ffef0a1"
ENV REPDENOVO_REPO="https://github.com/Reedwarbler/REPdenovo.git"

RUN  apt-get update \
  && apt-get install -y \
    build-essential \
    git \
    zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*


# Removing static compile flag because gives weird warning.
# Move shebang to top of file because __author__ is at line 1.
# This prevents using without python executable first.
WORKDIR ${REPDENOVO_PREFIX}
RUN  git clone ${REPDENOVO_REPO} . \
  && git checkout ${REPDENOVO_VERSION} \
  && cd ${REPDENOVO_PREFIX}/TERefiner \
  && sed -i "s~../bamtools-master/include~${BAMTOOLS_PREFIX}/include/bamtools~" Makefile \
  && sed -i "s~../bamtools-master/lib~${BAMTOOLS_PREFIX}/lib~" Makefile \
  && sed -i 's/-static//' Makefile \
  && make \
  && cp TERefiner_1 ${REPDENOVO_PREFIX} \
  && cd ${REPDENOVO_PREFIX}/ContigsCompactor-v0.2.0/ContigsMerger \
  && make \
  && cp ContigsMerger ${REPDENOVO_PREFIX} \
  && chmod a+x ${REPDENOVO_PREFIX}/ContigsMerger ${REPDENOVO_PREFIX}/TERefiner_1 \
  && rm -rf -- ${REPDENOVO_PREFIX}/TERefiner ${REPDENOVO_PREFIX}/ContigsCompactor-v0.2.0 \
  && chmod a+x ${REPDENOVO_PREFIX}/main.py \
  && sed -i '/^#!/d' ${REPDENOVO_PREFIX}/main.py \
  && sed -i '1 i #!/usr/bin/env python' ${REPDENOVO_PREFIX}/main.py \
  && rm ${REPDENOVO_PREFIX}/*.txt ${REPDENOVO_PREFIX}/*.pdf


FROM debian:${DEBIAN_VERSION}

RUN  apt-get update \
  && apt-get install -y \
       libbz2-1.0 \
       libcurl3 \
       libgomp1 \
       libgsl2 \
       libncurses5 \
       libperl5.24 \
       libssl1.1 \
       lzma \
       perl \
       python \
       zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV HTSLIB_PREFIX="/opt/htslib"
ENV BCFTOOLS_PREFIX="/opt/bcftools"
ENV SAMTOOLS_PREFIX="/opt/samtools"
COPY --from="darcyabjones/htslib" "${HTSLIB_PREFIX}" "${HTSLIB_PREFIX}"
COPY --from="darcyabjones/htslib" "${BCFTOOLS_PREFIX}" "${BCFTOOLS_PREFIX}"
COPY --from="darcyabjones/htslib" "${SAMTOOLS_PREFIX}" "${SAMTOOLS_PREFIX}"

ENV VELVET_PREFIX="/opt/velvet"
COPY --from="darcyabjones/velvet" ${VELVET_PREFIX} ${VELVET_PREFIX}

ENV JELLYFISH_PREFIX="/opt/jellyfish"
COPY --from="darcyabjones/jellyfish" ${JELLYFISH_PREFIX} ${JELLYFISH_PREFIX}

ENV BWA_PREFIX="/opt/bwa"
COPY --from="darcyabjones/bwa" ${BWA_PREFIX} ${BWA_PREFIX}

ENV BAMTOOLS_PREFIX="/opt/bamtools"
COPY --from="darcyabjones/bamtools" ${BAMTOOLS_PREFIX} ${BAMTOOLS_PREFIX}

ENV REPDENOVO_PREFIX="/opt/repdenovo"
COPY --from=builder ${REPDENOVO_PREFIX} ${REPDENOVO_PREFIX}

ENV PATH="${PATH}:${SAMTOOLS_PREFIX}/bin:${BCFTOOLS_PREFIX}/bin:${HTSLIB_PREFIX}/bin:${JELLYFISH_PREFIX}:${VELVET_PREFIX}:${BWA_PREFIX}:${BAMTOOLS_PREFIX}/bin:${REPDENOVO_PREFIX}"

ENV LIBRARY_PATH="${JELLYFISH_DIR}/lib:${BAMTOOLS_PREFIX}/lib:${HTSLIB_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${LIBRARY_PATH}:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${LIBRARY_PATH}:${LD_RUN_PATH}"

ENV CPATH="${HTSLIB_PREFIX}/include:${BAMTOOLS_PREFIX}/include:${JELLYFISH_DIR}/include:${CPATH}"
ENV INCLUDE="${CPATH}:${INCLUDE}"


ENTRYPOINT ["main.py"]
CMD ["--help"]
