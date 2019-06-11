ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION}

ENV DEBIAN_FRONTEND=noninteractive

# Set these with empty defaults to avoid using unset variables in path adds
ENV PATH "${PATH:-}"
ENV INCLUDE "${INCLUDE:-}"
ENV CPATH "${CPATH:-}"
ENV LIBRARY_PATH "${LIBRARY_PATH:-}"
ENV LD_LIBRARY_PATH "${LD_LIBRARY_PATH:-}"
ENV LD_RUN_PATH "${LD_RUN_PATH:-}"

ENV LANG "en_US.UTF-8"
ENV LANGUAGE "en_US.UTF-8"
ENV LC_ALL "en_US.UTF-8"

ENV BUILD_DIR "/build"
ENV ENV_FILE "${BUILD_DIR}/env.sh"
ENV ENV_DIR "${BUILD_DIR}/env"
ENV APT_REQUIREMENTS_FILE "${BUILD_DIR}/apt-requirements.txt"

COPY base.sh /build/base.sh

# This is the file that dash will source on startup
ENV ENV /build/base.sh

# This is the file that bash will source on startup for non-interactive.
ENV BASH_ENV "${ENV}"

RUN  apt-get update \
  && apt-get install -y --no-install-recommends \
       locales \
  && rm -rf /var/lib/apt/lists/* \
  && sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
  && dpkg-reconfigure --frontend=noninteractive locales \
  && update-locale LANG="${LANG}" \
  && touch "${ENV_FILE}" "${APT_REQUIREMENTS_FILE}" \
  && echo ". \${ENV}" >> /etc/bash.bashrc
