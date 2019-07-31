ARG IMAGE

FROM "${IMAGE}" as exonerate_builder

ARG EXONERATE_VERSION
ARG EXONERATE_URL="http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.4.0.tar.gz"
ARG EXONERATE_PREFIX_ARG="/opt/exonerate/${EXONERATE_VERSION}"
ENV EXONERATE_PREFIX="${EXONERATE_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       libglib2.0-dev \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O exonerate.tar.gz "${EXONERATE_URL}" \
  && tar -zxf exonerate.tar.gz \
  && cd exonerate*/ \
  && ./configure --prefix="${EXONERATE_PREFIX}" \
  && make \
  && make install \
  && add_runtime_dep libglib2.0


FROM "${IMAGE}"

ARG EXONERATE_VERSION
ARG EXONERATE_PREFIX_ARG="/opt/exonerate/${EXONERATE_VERSION}"
ENV EXONERATE_PREFIX="${EXONERATE_PREFIX_ARG}"

LABEL exonerate.version="${EXONERATE_VERSION}"

ENV PATH "${EXONERATE_PREFIX}/bin:${PATH}"

COPY --from=exonerate_builder "${EXONERATE_PREFIX}" "${EXONERATE_PREFIX}"
COPY --from=exonerate_builder "${APT_REQUIREMENTS_FILE}" /build/apt/exonerate.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
