ARG IMAGE

FROM "${IMAGE}" as builder

ARG FASTA_VERSION
ARG FASTA_URL="http://faculty.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.8g.tar.gz"
ARG FASTA_PREFIX_ARG="/opt/fasta/${FASTA_VERSION}"
ENV FASTA_PREFIX="${FASTA_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget -O fasta.tar.gz "${FASTA_URL}" \
  && tar -zxf fasta.tar.gz \
  && cd fasta*/src \
  && make -f ../make/Makefile.linux64_sse2 all \
  && cd .. \
  && rm -f bin/README \
  && mkdir -p "${FASTA_PREFIX}" \
  && cp -r bin "${FASTA_PREFIX}/bin" \
  && cp -r conf "${FASTA_PREFIX}/conf" \
  && cp -r data "${FASTA_PREFIX}/data" \
  && cp -r psisearch2 "${FASTA_PREFIX}/psisearch2" \
  && ln -sf -- "${FASTA_PREFIX}/bin/fasta36" "${FASTA_PREFIX}/bin/fasta" \
  && add_runtime_dep perl python

# Symbolic link of fasta executable is necessary for pasa

from "${IMAGE}"

ARG FASTA_VERSION
ARG FASTA_PREFIX_ARG="/opt/fasta/${FASTA_VERSION}"
ENV FASTA_PREFIX="${FASTA_PREFIX_ARG}"
LABEL fasta.version="${FASTA_VERSION}"

ENV PATH="${FASTA_PREFIX}/bin:${FASTA_PREFIX}/psisearch2:${PATH}"

COPY --from=builder "${FASTA_PREFIX}" "${FASTA_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/fasta.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
