ARG IMAGE

FROM "${IMAGE}" as builder

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
       wget \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
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

COPY --from=builder "${STAR_PREFIX}" "${STAR_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/star.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
