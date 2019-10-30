ARG IMAGE
ARG GENEMARKES_IMAGE
ARG SIGNALP5_IMAGE

FROM "${GENEMARKES_IMAGE}" as genemarkes_builder
FROM "${SIGNALP5_IMAGE}" as signalp5_builder


FROM "${IMAGE}"

ARG GENEMARKES_VERSION
ARG GENEMARKES_PREFIX_ARG="/opt/genemarkes/${GENEMARKES_VERSION}"
ENV GENEMARKES_PREFIX="${GENEMARKES_PREFIX_ARG}"
LABEL genemarkes.version="${GENEMARKES_VERSION}"

ENV PATH="${GENEMARKES_PREFIX}:${PATH}"

COPY --from=genemarkes_builder "${GENEMARKES_PREFIX}" "${GENEMARKES_PREFIX}"
COPY --from=genemarkes_builder "${APT_REQUIREMENTS_FILE}" /build/apt/genemarkes.txt


ARG SIGNALP5_VERSION
ARG SIGNALP5_PREFIX_ARG="/opt/signalp/${SIGNALP_VERSION}"
ARG SIGNALP5_TAR="sources/signalp-5.0.Linux.tar.gz"
ENV SIGNALP5_PREFIX="${SIGNALP5_PREFIX_ARG}"
ENV PATH="${SIGNALP5_PREFIX}/bin:${PATH}"
LABEL signalp5.version="${SIGNALP5_VERSION}"

COPY --from=signalp5_builder "${SIGNALP5_PREFIX}" "${SIGNALP5_PREFIX}"
COPY --from=signalp5_builder "${APT_REQUIREMENTS_FILE}" /build/apt/signalp5.txt


RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"

WORKDIR /
