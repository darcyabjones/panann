ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

ENV LTRDETECTOR_PREFIX="/opt/ltrdetector"
ARG LTRDETECTOR_VERSION="e50cfb7"
ENV LTRDETECTOR_REPO="https://github.com/TulsaBioinformaticsToolsmith/LtrDetector.git"

WORKDIR /tmp
RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       git \
  && rm -rf /var/lib/apt/lists/* \
  && git clone ${LTRDETECTOR_REPO} . \
  && git checkout ${LTRDETECTOR_VERSION} \
  && cd src \
  && make bin \
  && make tr -j \
  && mkdir ${LTRDETECTOR_PREFIX} \
  && mkdir ${LTRDETECTOR_PREFIX}/bin \
  && cp /tmp/bin/LtrDetector ${LTRDETECTOR_PREFIX}/bin \
  && cp /tmp/visualize.py ${LTRDETECTOR_PREFIX}/bin \
  && cp /tmp/requirements.txt ${LTRDETECTOR_PREFIX} \
  && chmod a+x ${LTRDETECTOR_PREFIX}/bin/visualize.py

FROM debian:${DEBIAN_VERSION}
ENV LTRDETECTOR_PREFIX="/opt/ltrdetector"

COPY --from=builder ${LTRDETECTOR_PREFIX} ${LTRDETECTOR_PREFIX}

RUN  apt-get update \
  && apt-get install -y libgomp1 \
  && rm -rf /var/lib/apt/lists/*

# Not installing python and dependencies because the visualise script
# doesn't do much.

ENV PATH="${LTRDETECTOR_PREFIX}/bin:${PATH}"
