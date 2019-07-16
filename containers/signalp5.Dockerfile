ARG IMAGE="darcyabjones/base"

FROM "${IMAGE}"

ARG SIGNALP5_VERSION
ARG SIGNALP5_PREFIX_ARG="/opt/signalp/${SIGNALP_VERSION}"
ARG SIGNALP5_TAR="sources/signalp-5.0.Linux.tar.gz"
ENV SIGNALP5_PREFIX="${SIGNALP5_PREFIX_ARG}"
ENV PATH="${SIGNALP5_PREFIX}/bin:${PATH}"
LABEL signalp5.version="${SIGNALP5_VERSION}"

# Signalp needs to be called as the full path to the executable.
# So the shell script just wraps it so we can call it on the PATH
COPY "${SIGNALP5_TAR}" /tmp/signalp.tar.gz
WORKDIR /tmp
RUN  set -eu \
  && . /build/base.sh \
  && tar -zxf signalp.tar.gz \
  && rm signalp.tar.gz \
  && mkdir "${SIGNALP5_PREFIX%/*}" \
  && mv signalp* "${SIGNALP5_PREFIX}" \
  && cd "${SIGNALP5_PREFIX}" \
  && mv "${SIGNALP5_PREFIX}/bin" "${SIGNALP5_PREFIX}/exe" \
  && mkdir -p "${SIGNALP5_PREFIX}/bin" \
  && echo "#!/usr/bin/env sh" > "${SIGNALP5_PREFIX}/bin/signalp" \
  && echo '${SIGNALP5_PREFIX}/exe/signalp $*' >> "${SIGNALP5_PREFIX}/bin/signalp" \
  && chmod a+x "${SIGNALP5_PREFIX}/bin/signalp" \
  && ln -sf "${SIGNALP5_PREFIX}/bin/signalp" "${SIGNALP5_PREFIX}/bin/signalp-${SIGNALP5_VERSION}" \
  && rm -rf -- "${SIGNALP5_PREFIX}/signalp.1" "${SIGNALP5_PREFIX}/test"

WORKDIR /
