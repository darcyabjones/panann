ARG IMAGE
ARG GENOMETOOLS_IMAGE

FROM "${GENOMETOOLS_IMAGE}" as genometools_builder

FROM "${IMAGE}" as builder

ARG AEGEAN_VERSION
ARG AEGEAN_URL="https://github.com/standage/AEGeAn/archive/v0.15.0.tar.gz"
ARG AEGEAN_PREFIX_ARG="/opt/aegean/${AEGEAN_VERSION}"
ENV AEGEAN_PREFIX="${AEGEAN_PREFIX_ARG}"

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

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       build-essential \
       ca-certificates \
       libcairo2-dev \
       libpango1.0-dev \
       python \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O aegean.tar.gz "${AEGEAN_URL}" \
  && tar zxf aegean.tar.gz \
  && rm aegean.tar.gz \
  && cd AEGeAn*/ \
  && sed -i "s~/usr/local/include/genometools~${GENOMETOOLS_PREFIX}/include/genometools~" Makefile \
  && make test \
  && make prefix="${AEGEAN_PREFIX}" install install-scripts \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"


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

COPY --from=builder "${AEGEAN_PREFIX}" "${AEGEAN_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/aegean.txt


ARG GENOMETOOLS_VERSION
ARG GENOMETOOLS_PREFIX_ARG="/opt/genometools/${GENOMETOOLS_VERSION}"
ENV GENOMETOOLS_PREFIX="${GENOMETOOLS_PREFIX_ARG}"
LABEL genometools.version="${GENOMETOOLS_VERSION}"

ENV PATH="${GENOMETOOLS_PREFIX}/bin:${PATH}"
ENV INCLUDE="${GENOMETOOLS_PREFIX}/include:${INCLUDE}"
ENV CPATH="${GENOMETOOLS_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${GENOMETOOLS_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${GENOMETOOLS_PREFIX}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${GENOMETOOLS_PREFIX}/lib:${LD_RUN_PATH}"

COPY --from=genometools_builder "${GENOMETOOLS_PREFIX}" "${GENOMETOOLS_PREFIX}"
COPY --from=genometools_builder "${APT_REQUIREMENTS_FILE}" /build/apt/genometools.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
