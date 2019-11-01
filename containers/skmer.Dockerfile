ARG IMAGE
ARG JELLYFISH_IMAGE
ARG SEQTK_IMAGE

FROM "${JELLYFISH_IMAGE}" as jellyfish_builder
FROM "${SEQTK_IMAGE}" as seqtk_builder

FROM "${IMAGE}" as skmer_builder

ARG SKMER_COMMIT
ARG SKMER_REPO="https://github.com/shahab-sarmashghi/Skmer.git"
ARG SKMER_PREFIX_ARG="/opt/skmer/${SKMER_PREFIX}"
ENV SKMER_PREFIX="${SKMER_PREFIX_ARG}"

ARG MASH_VERSION
ARG MASH_URL="https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar"


WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       git \
       python3 \
       python3-pip \
       python3-setuptools \
       python3-numpy \
       python3-pandas \
       python3-scipy \
       python3-wheel \
       wget \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${SKMER_REPO}" . \
  && git fetch --tags \
  && git checkout "${SKMER_COMMIT}" \
  && pip3 install --prefix="${SKMER_PREFIX}" . \
  && add_python3_site "${SKMER_PREFIX}/lib/python3.7/site-packages" \
  && wget -c -O mash.tar "${MASH_URL}" \
  && tar -xf mash.tar \
  && mv mash-*/mash "${SKMER_PREFIX}/bin" \
  && add_runtime_dep \
       python3 \
       python3-biopython \
       python3-numpy \
       python3-pandas \
       python3-scipy


FROM "${IMAGE}"

ARG SKMER_COMMIT
ARG SKMER_PREFIX_ARG="/opt/skmer/${SKMER_PREFIX}"
ENV SKMER_PREFIX="${SKMER_PREFIX_ARG}"
LABEL skmer.version="${SKMER_COMMIT}"

ENV PATH "${SKMER_PREFIX}/bin:${PATH}"

COPY --from=skmer_builder "${SKMER_PREFIX}" "${SKMER_PREFIX}"
COPY --from=skmer_builder "${PYTHON3_SITE_PTH_FILE}" "${PYTHON3_SITE_DIR}/skmer.pth"
COPY --from=skmer_builder "${APT_REQUIREMENTS_FILE}" /build/apt/skmer.txt


ARG JELLYFISH_VERSION
ARG JELLYFISH_PREFIX_ARG="/opt/jellyfish/${JELLYFISH_VERSION}"
ENV JELLYFISH_PREFIX="${JELLYFISH_PREFIX_ARG}"
LABEL jellyfish.version="${JELLYFISH_VERSION}"

ENV PATH="${JELLYFISH_PREFIX}/bin:${PATH}"
ENV INCLUDE="${JELLYFISH_DIR}/include:${INCLUDE}"
ENV CPATH="${JELLYFISH_DIR}/include:${CPATH}"
ENV LIBRARY_PATH="${JELLYFISH_DIR}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${JELLYFISH_DIR}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${JELLYFISH_DIR}/lib:${LD_RUN_PATH}"

COPY --from=jellyfish_builder "${JELLYFISH_PREFIX}" "${JELLYFISH_PREFIX}"


ARG SEQTK_COMMIT
ARG SEQTK_PREFIX_ARG="/opt/seqtk/${SEQTK_COMMIT}"
ENV SEQTK_PREFIX="${SEQTK_PREFIX_ARG}"
LABEL seqtk.version="${SEQTK_COMMIT}"

ENV PATH "${SEQTK_PREFIX}/bin:${PATH}"

COPY --from=seqtk_builder "${SEQTK_PREFIX}" "${SEQTK_PREFIX}"
COPY --from=seqtk_builder "${APT_REQUIREMENTS_FILE}" /build/apt/seqtk.txt


RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat "${PYTHON3_SITE_DIR}/skmer.pth" >> "${PYTHON3_SITE_PTH_FILE}" \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"

WORKDIR /
