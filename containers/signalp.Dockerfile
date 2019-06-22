ARG IMAGE="darcyabjones/base"

FROM "${IMAGE}"

ARG SIGNALP_VERSION
ARG SIGNALP_PREFIX_ARG="/opt/signalp/${SIGNALP_VERSION}"
ARG SIGNALP_TAR="sources/signalp-5.0.Linux.tar.gz"
ENV SIGNALP_PREFIX="${SIGNALP_PREFIX_ARG}"


# Signalp needs to be called as the full path to the executable.
# So the shell script just wraps it so we can call it on the PATH
COPY "${SIGNALP_TAR}" /opt/signalp/signalp.tar.gz
WORKDIR /opt/signalp
RUN  set -eu \
  && . /build/base.sh \
  && tar -zxf signalp.tar.gz \
  && rm signalp.tar.gz \
  && mv signalp* "${SIGNALP_VERSION}" \
  && cd "${SIGNALP_PREFIX}" \
  && mv "${SIGNALP_PREFIX}/bin" "${SIGNALP_PREFIX}/exe" \
  && mkdir -p "${SIGNALP_PREFIX}/bin" \
  && echo "#!/usr/bin/env sh" > "${SIGNALP_PREFIX}/bin/signalp" \
  && echo '${SIGNALP_PREFIX}/exe/signalp $*' >> "${SIGNALP_PREFIX}/bin/signalp" \
  && chmod a+x "${SIGNALP_PREFIX}/bin/signalp" \
  && rm -rf -- "${SIGNALP_PREFIX}/signalp.1" "${SIGNALP_PREFIX}/test" \
  && apt-get update \
  && apt-get install -y gawk \

WORKDIR /
ENV PATH="${PATH}:${SIGNALP_PREFIX}/bin"
