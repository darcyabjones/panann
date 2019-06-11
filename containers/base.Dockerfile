ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION}

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG en_US.UTF-8

RUN  apt-get update \
  && apt-get install -y --no-install-recommends \
       locales \
  && rm -rf /var/lib/apt/lists/* \
  && sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
  && dpkg-reconfigure --frontend=noninteractive locales \
  && update-locale LANG=en_US.UTF-8
