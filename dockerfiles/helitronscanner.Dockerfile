FROM openjdk:8u212-jre-slim-stretch

LABEL maintainer="darcy.ab.jones@gmail.com"

ENV HELITRONSCANNER_PREFIX="/opt/helitronscanner"
ENV HELITRONSCANNER_URL="http://bo.csam.montclair.edu/du/assets/filesdb/HelitronScanner.zip"
# Alternative is at sourceforge.

# The man thing is required to get java installed on the minimal distro
WORKDIR ${HELITRONSCANNER_PREFIX}
RUN  apt-get update \
  && apt-get install -y \
       unzip \
       wget \
  && rm -rf -- /var/lib/apt/lists/* \
  && wget ${HELITRONSCANNER_URL} \
  && unzip *.zip \
  && mv HelitronScanner bin \
  && mkdir data \
  && mv TrainingSet/*.lcvs data \
  && rm *.zip \
  && rm -rf -- __MACOSX TrainingSet \
  && rm *.pdf

ENV HELITRONSCANNER="${HELITRONSCANNER_PREFIX}/bin/HelitronScanner.jar"
ENV HELITRONSCANNER_DATA="${HELITRONSCANNER_PREFIX}/data"
ENV PATH="${HELITRONSCANNER_PREFIX}/bin:${PATH}"

WORKDIR /
ENTRYPOINT ["java", "-jar", "/opt/helitronscanner/bin/HelitronScanner.jar"]
CMD ["-help"]
