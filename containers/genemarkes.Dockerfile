ARG IMAGE

FROM "${IMAGE}"

ARG GENEMARKES_VERSION
ARG GENEMARKES_TAR
ARG GENEMARKES_KEY
ARG GENEMARKES_PREFIX_ARG="/opt/genemarkes/${GENEMARKES_VERSION}"
ENV GENEMARKES_PREFIX="${GENEMARKES_PREFIX_ARG}"
LABEL genemarkes.version="${GENEMARKES_VERSION}"

ENV PATH="${GENEMARKES_PREFIX}:${PATH}"

COPY "${GENEMARKES_TAR}" /tmp/genemark.tar.gz
COPY "${GENEMARKES_KEY}" /root/.gm_key


# CPAN package is in debian unstable still
# liblogger-simple-perl

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && add_runtime_dep \
       perl \
       libyaml-perl \
       libhash-merge-perl \
       libparallel-forkmanager-perl \
  && apt-get update \
  && apt_install_from_file "${APT_REQUIREMENTS_FILE}" \
  && apt-get install -y --no-install-recommends \
       cpanminus \
       make \
       tar \
  && rm -rf /var/lib/apt/lists/* \
  && cpanm Logger::Simple \
  && tar -zxf /tmp/genemark.tar.gz \
  && mkdir -p "${GENEMARKES_PREFIX}" \
  && cp /root/.gm_key "${GENEMARKES_PREFIX}/gm_key" \
  && cp -r gm_et_linux_64/gmes_petap/* "${GENEMARKES_PREFIX}" \
  && rm -rf -- gm_et_linux_64 genemark.tar.gz

WORKDIR /
