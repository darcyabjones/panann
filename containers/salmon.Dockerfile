ARG IMAGE

FROM "${IMAGE}" as builder

ARG SALMON_TAG="v0.13.1"
ARG SALMON_REPO="https://github.com/COMBINE-lab/salmon.git"
ARG SALMON_PREFIX_ARG="/opt/salmon/${SALMON_TAG}"
ENV SALMON_PREFIX="${SALMON_PREFIX_ARG}"

WORKDIR /tmp
RUN  echo "deb http://deb.debian.org/debian stretch-backports main" >> /etc/apt/sources.list \
  && apt-get update \
  && apt-get install -y \
       autoconf \
       build-essential \
       curl \
       git \
       libbz2-dev \
       libboost-all-dev \
       liblzma-dev \
       libtbb-dev \
       unzip \
       zlib1g-dev \
  && apt-get -t stretch-backports install -y cmake \
  && rm -rf /var/lib/apt/lists/* \
  && git clone "${SALMON_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${SALMON_TAG}" \
  && mkdir -p build \
  && cd build \
  && cmake -DCMAKE_INSTALL_PREFIX="${SALMON_PREFIX}" .. \
  && make \
  && make install \
  && make test

RUN  echo 'libtbb2 libbz2-1.0 zlib1g liblzma5' >> "${SALMON_PREFIX}/apt-requirements.txt" \
  && echo 'export PATH="${SALMON_PREFIX}/bin:${PATH}"' >> "${SALMON_PREFIX}/environment.sh" \
  && echo 'export LD_LIBRARY_PATH="${SALMON_PREFIX}/lib:${LD_LIBRARY_PATH}"' >> "${SALMON_PREFIX}/environment.sh"


FROM "${IMAGE}"

ARG SALMON_TAG
ARG SALMON_PREFIX_ARG="/opt/salmon/${SALMON_TAG}"
ENV SALMON_PREFIX="${SALMON_PREFIX_ARG}"

COPY --from=builder "${SALMON_PREFIX}" "${SALMON_PREFIX}"

RUN  apt-get update \
  && xargs -a "${SALMON_PREFIX}/apt-requirements.txt" -r -- apt-get install -y --no-install-recommends \
  && rm -rf /var/lib/apt/lists/*

ENV PATH="${SALMON_PREFIX}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${SALMON_PREFIX}/lib:${LD_LIBRARY_PATH}"
