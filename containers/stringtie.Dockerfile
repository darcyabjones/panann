ARG IMAGE

FROM "${IMAGE}" as builder

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
ENV  PATH "${STRINGTIE_PREFIX}/bin:${PATH}"

COPY --from=builder "${STRINGTIE_PREFIX}" "${STRINGTIE_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/stringtie.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
