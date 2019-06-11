ARG IMAGE

FROM "${IMAGE}" as builder

ARG STAR_VERSION
ARG STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz"
ARG STAR_PREFIX_ARG="/opt/star/${STAR_VERSION}"
ENV STAR_PREFIX="${STAR_PREFIX_ARG}"

WORKDIR /tmp
RUN  apt-get update \
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
  && mv ../bin/Linux_x86_64/STAR* "${STAR_PREFIX}/bin"

# See if we can provide build flags as a build time-option?
# STAR already seems to require x86-64 instructions, so i think current settings are safe.
# LDFLAGSextra=-flto CXXFLAGSextra="-flto -march=x86-64 -mtune=native"

# Runtime requires
#  libgomp1
#  zlib1g-dev

FROM "${IMAGE}"

ARG STAR_VERSION
ARG STAR_PREFIX_ARG="/opt/star/${STAR_VERSION}"
ENV STAR_PREFIX="${STAR_PREFIX_ARG}"

COPY --from=builder "${STAR_PREFIX}" "${STAR_PREFIX}"

RUN  apt-get update \
  && apt-get install -y libgomp1 zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV PATH="${STAR_PREFIX}/bin:${PATH}"
