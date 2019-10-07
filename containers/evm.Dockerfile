ARG IMAGE


FROM "${IMAGE}" as evm_builder

ARG EVM_COMMIT
ARG EVM_REPO
ARG EVM_PREFIX_ARG
ENV EVM_PREFIX="${EVM_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       git \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${EVM_REPO}" . \
  && git checkout "${EVM_COMMIT}" \
  && mkdir -p "${EVM_PREFIX}" \
  && cp -r EvmUtils "${EVM_PREFIX}" \
  && cp -r PerlLib "${EVM_PREFIX}" \
  && cp evidence_modeler.pl "${EVM_PREFIX}" \
  && chmod a+x "${EVM_PREFIX}/evidence_modeler.pl" \
  && chmod a+x "${EVM_PREFIX}"/EvmUtils/*.pl \
  && chmod a+x "${EVM_PREFIX}"/EvmUtils/misc/*.pl \
  && rm -rf -- "${EVM_PREFIX}/EvmUtils/misc/example_data_files" \
  && add_runtime_dep libdbi-perl liburi-perl perl perl-tk

#libfindbin-libs-perl

FROM "${IMAGE}"

ARG EVM_COMMIT
ARG EVM_PREFIX_ARG
ENV EVM_PREFIX="${EVM_PREFIX_ARG}"

LABEL evm.version="${EVM_COMMIT}"

ENV PATH "${EVM_PREFIX}:${EVM_PREFIX}/EvmUtils:${EVM_PREFIX}/EvmUtils/misc:${PATH}"

COPY --from=evm_builder "${EVM_PREFIX}" "${EVM_PREFIX}"
COPY --from=evm_builder "${APT_REQUIREMENTS_FILE}" /build/apt/evm.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
