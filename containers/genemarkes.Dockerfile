ARG IMAGE

FROM "${IMAGE}"

ARG GENEMARKES_VERSION
ARG GENEMARKES_TAR
ARG GENEMARKES_KEY
ARG GENEMARKES_PREFIX_ARG="/opt/genemarkes/${GENEMARKES_VERSION}"
ENV GENEMARKES_PREFIX="${GENEMARKES_PREFIX_ARG}"

ENV PATH="${GENEMARKES_PREFIX}:${PATH}"

COPY "${GENEMARKES_TAR}" /tmp/genemark.tar.gz
COPY "${GENEMARKES_KEY}" /root/.gm_key

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && add_runtime_dep perl libyaml-perl libhash-merge-perl libparallel-forkmanager-perl cpanminus \
  && apt-get update \
  && apt_install_from_file "${APT_REQUIREMENTS_FILE}" \
  && apt-get install -y --no-install-recommends \
       make \
       tar \
  && rm -rf /var/lib/apt/lists/* \
  && cpanm Logger::Simple \
  && tar -zxf /tmp/genemark.tar.gz \
  && mkdir -p "${GENEMARKES_PREFIX}" \
  && cp -r gm_et_linux_64/gmes_petap/* "${GENEMARKES_PREFIX}" \
  && rm -rf -- gm_et_linux_64 genemark.tar.gz
