ARG IMAGE
ARG AUGUSTUS_IMAGE

FROM "${AUGUSTUS_IMAGE}" as augustus_builder

FROM "${IMAGE}" as builder

ARG BUSCO_COMMIT="1554283ab8ee7dd5b5290f4f748234f456c36e66"
ARG BUSCO_REPO="https://gitlab.com/ezlab/busco.git"
ARG BUSCO_PREFIX_ARG="/opt/busco/${BUSCO_COMMIT}"
ENV BUSCO_PREFIX="${BUSCO_PREFIX_ARG}"

ARG AUGUSTUS_COMMIT
ARG AUGUSTUS_PREFIX_ARG="/opt/augustus/${AUGUSTUS_COMMIT}"
ENV AUGUSTUS_PREFIX="${AUGUSTUS_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       ca-certificates \
       git \
       python3 \
       python3-distutils \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${BUSCO_REPO}" "${BUSCO_PREFIX}" \
  && cd "${BUSCO_PREFIX}" \
  && git checkout "${BUSCO_COMMIT}" \
  && sed -i 's~#!/usr/bin/env python~#!/usr/bin/env python3~g' scripts/run_BUSCO.py \
  && sed -i 's~#!/usr/bin/env python~#!/usr/bin/env python3~g' scripts/generate_plot.py \
  && cp config/config.ini.default config/config.ini \
  && sed -i "s~path = /home/osboxes/BUSCOVM/augustus/augustus-3.2.2/bin/~path = ${AUGUSTUS_PREFIX}/bin/~g" config/config.ini \
  && sed -i "s~path = /home/osboxes/BUSCOVM/augustus/augustus-3.2.2/scripts/~path = ${AUGUSTUS_PREFIX}/scripts/~g" config/config.ini \
  && sed -i 's~path = /home/osboxes/BUSCOVM/hmmer/hmmer-3.1b2-linux-intel-ia32/binaries/~path = /usr/bin/~g' config/config.ini \
  && python3 setup.py install --prefix=${PWD} \
  && rm -rf -- .git build src \
  && rm -rf -- *.pdf *.md setup.py CHANGELOG \
  && add_runtime_dep \
       hmmer \
       ncbi-blast+ \
       perl \
       python3 \
       r-base


FROM "${IMAGE}"

ARG BUSCO_COMMIT
ARG BUSCO_PREFIX_ARG="/opt/busco/${BUSCO_COMMIT}"
ENV BUSCO_PREFIX="${BUSCO_PREFIX_ARG}"
LABEL busco.version="${BUSCO_COMMIT}"

ENV PATH="${BUSCO_PREFIX}/scripts:${PATH}"
ENV PYTHONPATH="${PYTHONPATH}:${BUSCO_PREFIX}/lib/python3.7/site-packages"

COPY --from=builder "${BUSCO_PREFIX}" "${BUSCO_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/busco.txt

ARG AUGUSTUS_COMMIT
ARG AUGUSTUS_PREFIX_ARG="/opt/augustus/${AUGUSTUS_COMMIT}"
ENV AUGUSTUS_PREFIX="${AUGUSTUS_PREFIX_ARG}"
ENV AUGUSTUS_CONFIG_PATH="${AUGUSTUS_PREFIX}/config"
LABEL augustus.version="${AUGUSTUS_COMMIT}"

ENV PATH="${PATH}:${AUGUSTUS_PREFIX}/bin:${AUGUSTUS_PREFIX}/scripts"

COPY --from=augustus_builder "${AUGUSTUS_PREFIX}" "${AUGUSTUS_PREFIX}"
COPY --from=augustus_builder "${APT_REQUIREMENTS_FILE}" /build/apt/augustus.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
