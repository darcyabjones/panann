ARG IMAGE="darcyabjones/base"

FROM "${IMAGE}" as builder

ARG DEEPSIG_COMMIT="69e01cb"
ARG DEEPSIG_PREFIX_ARG="/opt/deepsig/${DEEPSIG_COMMIT}"
ARG DEEPSIG_REPO="https://github.com/BolognaBiocomp/deepsig.git"
ENV DEEPSIG_PREFIX="${DEEPSIG_PREFIX_ARG}"
LABEL deepsig.version="${DEEPSIG_COMMIT}"

ARG TENSORFLOW_VERSION
ARG TENSORFLOW_PREFIX_ARG="/opt/tensorflow/${TENSORFLOW_VERSION}"
ENV TENSORFLOW_PREFIX="${TENSORFLOW_PREFIX_ARG}"
LABEL tensorflow.version="TENSORFLOW_VERSION"

ARG KERAS_VERSION
ARG KERAS_PREFIX_ARG="/opt/keras/${KERAS_VERSION}"
ENV KERAS_PREFIX="${KERAS_PREFIX_ARG}"
LABEL keras.version="KERAS_VERSION"

WORKDIR /opt
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       ca-certificates \
       git \
       python-pip \
       python-setuptools \
       python-wheel \
       python \
       python-biopython \
       python-numpy \
       python-scipy \
       python-six \
       python-yaml \
       python-h5py \
  && apt-get autoremove \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* \
  && mkdir -p "${DEEPSIG_PREFIX%/*}" \
  && git clone "${DEEPSIG_REPO}" "${DEEPSIG_PREFIX}" \
  && cd "${DEEPSIG_PREFIX}" \
  && git checkout "${DEEPSIG_COMMIT}" \
  && rm -rf -- .git \
  && python -m pip install --prefix="${TENSORFLOW_PREFIX}" tensorflow=="${TENSORFLOW_VERSION}" \
  && PYTHONPATH="${TENSORFLOW_PREFIX}/lib/python2.7/site-packages:${PYTHONPATH:-}" \
     python -m pip install --prefix="${KERAS_PREFIX}" keras=="${KERAS_VERSION}" \
  && add_runtime_dep \
       python \
       python-biopython \
       python-numpy \
       python-scipy \
       python-six \
       python-yaml \
       python-h5py \
       python-protobuf


FROM "${IMAGE}"

ARG DEEPSIG_COMMIT="69e01cb"
ARG DEEPSIG_PREFIX_ARG="/opt/deepsig/${DEEPSIG_COMMIT}"
ENV DEEPSIG_PREFIX="${DEEPSIG_PREFIX_ARG}"
ENV DEEPSIG_ROOT="${DEEPSIG_PREFIX}"
LABEL deepsig.version="${DEEPSIG_COMMIT}"

ARG TENSORFLOW_VERSION
ARG TENSORFLOW_PREFIX_ARG="/opt/tensorflow/${TENSORFLOW_VERSION}"
ENV TENSORFLOW_PREFIX="${TENSORFLOW_PREFIX_ARG}"
ENV TF_CPP_MIN_LOG_LEVEL=3
LABEL tensorflow.version="TENSORFLOW_VERSION"

ARG KERAS_VERSION
ARG KERAS_PREFIX_ARG="/opt/keras/${KERAS_VERSION}"
ENV KERAS_PREFIX="${KERAS_PREFIX_ARG}"
LABEL keras.version="KERAS_VERSION"

ENV PATH="${DEEPSIG_PREFIX}:${PATH}"
ENV PATH="${TENSORFLOW_PREFIX}:${PATH}"
ENV PYTHONPATH="${TENSORFLOW_PREFIX}/lib/python2.7/site-packages:${PYTHONPATH:-}"
ENV PYTHONPATH="${KERAS_PREFIX}/lib/python2.7/site-packages:${PYTHONPATH:-}"

COPY --from=builder "${DEEPSIG_PREFIX}" "${DEEPSIG_PREFIX}"
COPY --from=builder "${TENSORFLOW_PREFIX}" "${TENSORFLOW_PREFIX}"
COPY --from=builder "${KERAS_PREFIX}" "${KERAS_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/deepsig.txt


RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
