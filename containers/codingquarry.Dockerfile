ARG IMAGE

FROM "${IMAGE}" as builder

ARG CODINGQUARRY_VERSION
ARG CODINGQUARRY_URL="https://downloads.sourceforge.net/project/codingquarry/CodingQuarry_v2.0.tar.gz"
ARG CODINGQUARRY_PREFIX_ARG="/opt/codingquarry/${CODINGQUARRY_VERSION}"
ENV CODINGQUARRY_PREFIX="${CODINGQUARRY_PREFIX_ARG}"
ENV QUARRY_PATH="${CODINGQUARRY_PREFIX}/QuarryFiles"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
      build-essential \
      wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget -O codingquarry.tar.gz "${CODINGQUARRY_URL}" \
  && tar -zxf codingquarry.tar.gz \
  && cd CodingQuarry* \
  && make -f makefile \
  && mkdir -p "${CODINGQUARRY_PREFIX}/bin" \
  && cp -r CodingQuarry "${CODINGQUARRY_PREFIX}/bin" \
  && cp -r run_CQ-PM_stranded.sh "${CODINGQUARRY_PREFIX}/bin" \
  && cp -r run_CQ-PM_unstranded.sh "${CODINGQUARRY_PREFIX}/bin" \
  && cp -r run_CQ-PM_mine.sh "${CODINGQUARRY_PREFIX}/bin" \
  && cp -r CufflinksGTF_to_CodingQuarryGFF3.py "${CODINGQUARRY_PREFIX}/bin" \
  && cp -r QuarryFiles "${QUARRY_PATH}" \
  && add_runtime_dep libgomp1 python python-biopython gawk


FROM "${IMAGE}"
ARG CODINGQUARRY_VERSION
ARG CODINGQUARRY_URL="https://downloads.sourceforge.net/project/codingquarry/CodingQuarry_v2.0.tar.gz"
ARG CODINGQUARRY_PREFIX_ARG="/opt/codingquarry/${CODINGQUARRY_VERSION}"
ENV CODINGQUARRY_PREFIX="${CODINGQUARRY_PREFIX_ARG}"
ENV QUARRY_PATH="${CODINGQUARRY_PREFIX}/QuarryFiles"

ENV PATH "${CODINGQUARRY_PREFIX}/bin:${PATH}"

COPY --from=builder "${CODINGQUARRY_PREFIX}" "${CODINGQUARRY_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/codingquarry.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
