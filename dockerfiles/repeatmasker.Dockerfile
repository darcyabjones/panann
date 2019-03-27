ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as nsegbuilder

ENV NSEG_PREFIX="/opt/nseg"
ENV NSEG_URL="ftp://ftp.ncbi.nih.gov/pub/seg/nseg/"

WORKDIR ${NSEG_PREFIX}
RUN  apt-get update \
  && apt-get install -y build-essential wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget -r -nd "${NSEG_URL}" \
  && make \
  && chmod a+x runnseg \
  && rm *.c *.h *.o makefile

ENV PATH="${PATH}:${NSEG_PREFIX}"


FROM debian:${DEBIAN_VERSION} as repscoutbuilder

ENV REPEATSCOUT_PREFIX="/opt/repeatscout"
ENV REPEATSCOUT_URL="http://www.repeatmasker.org/RepeatScout-1.0.5.tar.gz"
WORKDIR ${REPEATSCOUT_PREFIX}
RUN  apt-get update \
  && apt-get install -y build-essential wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget ${REPEATSCOUT_URL} \
  && tar zxf RepeatScout* \
  && rm RepeatScout*.tar.gz \
  && mv RepeatScout*/* ${REPEATSCOUT_PREFIX} \
  && rmdir RepeatScout*/ \
  && make \
  && rm *.c *.h *.o *.fa *.freq Makefile notebook README

ENV PATH="${PATH}:${REPEATSCOUT_PREFIX}"


FROM debian:${DEBIAN_VERSION} as reconbuilder

ENV RECON_PREFIX="/opt/recon"
ENV RECON_URL="http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz"
WORKDIR ${RECON_PREFIX}
RUN  apt-get update \
  && apt-get install -y build-essential wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget ${RECON_URL} \
  && tar zxf RECON* \
  && rm RECON*.tar.gz \
  && mv RECON*/* ./ \
  && rmdir RECON*/ \
  && cd src \
  && make install \
  && cd ${RECON_PREFIX} \
  && rm -rf -- src Demos \
  && sed -i "s~path = \"\"~path = \"${RECON_PREFIX}/bin\"~" scripts/recon.pl

ENV PATH="${PATH}:${RECON_PREFIX}/bin:${RECON_PREFIX}/scripts"


FROM debian:${DEBIAN_VERSION}

# Install blast and repeat modeller deps from other images.
ENV RMBLAST_PREFIX="/opt/rmblast"
COPY --from="darcyabjones/rmblast" "${RMBLAST_PREFIX}" "${RMBLAST_PREFIX}"

ENV NSEG_PREFIX="/opt/nseg"
COPY --from=nsegbuilder "${NSEG_PREFIX}" "${NSEG_PREFIX}"

ENV REPEATSCOUT_PREFIX="/opt/repeatscout"
COPY --from=repscoutbuilder "${REPEATSCOUT_PREFIX}" "${REPEATSCOUT_PREFIX}"

ENV RECON_PREFIX="/opt/recon"
COPY --from=reconbuilder "${RECON_PREFIX}" "${RECON_PREFIX}"


# Set some working environment variables.
ENV RM_LIB="/opt/rmlib"

ENV RMASK_PREFIX="/opt/repeatmasker"
ARG RMASK_VERSION="4-0-8"
ARG RMASK_URL="http://www.repeatmasker.org/RepeatMasker-open-${RMASK_VERSION}.tar.gz"

ENV RMOD_PREFIX="/opt/repeatmodeller"
ARG RMOD_VERSION="1.0.11"
ARG RMOD_URL="http://www.repeatmasker.org/RepeatModeler/RepeatModeler-open-${RMOD_VERSION}.tar.gz"

ENV TRF_PREFIX="/usr/local/bin"

# Install blast/rm dependencies, download trf, download repmodeler and repmask
WORKDIR /tmp
RUN  apt-get update \
  && apt-get install -y \
       hmmer \
       libblas3 \
       libbz2-1.0 \
       libcurl3-gnutls \
       libgmp10 \
       libgnutlsxx28 \
       libgomp1 \
       libjson-perl \
       liblapack3 \
       liblwp-useragent-determined-perl \
       libstdc++6 \
       libtext-soundex-perl \
       liburi-perl \
       libxml2 \
       libxslt1.1 \
       lzma \
       perl \
       wget \
       zlib1g \
  && rm -rf /var/lib/apt/lists/* \
  && wget http://tandem.bu.edu/trf/downloads/trf409.linux64 \
  && mv trf409.linux64 ${TRF_PREFIX} \
  && chmod a+x ${TRF_PREFIX}/trf409.linux64 \
  && ln -sf ${TRF_PREFIX}/trf409.linux64 ${TRF_PREFIX}/trf \
  && wget ${RMASK_URL} \
  && tar zxf RepeatMasker*.tar.gz \
  && rm RepeatMasker*.tar.gz \
  && mv RepeatMasker*/ ${RMASK_PREFIX} \
  && wget ${RMOD_URL} \
  && tar zxf RepeatModeler-*.tar.gz \
  && rm RepeatModeler*.tar.gz \
  && mv RepeatModeler*/ ${RMOD_PREFIX}

# Configure repeatmasker scripts
WORKDIR ${RMASK_PREFIX}
RUN  cp RepeatMaskerConfig.tmpl RepeatMaskerConfig.pm \
  && sed -i "s~TRF_PRGM\s*=\s*\"\"~TRF_PRGM = \"${TRF_PREFIX}/trf\"~" RepeatMaskerConfig.pm \
  && sed -i 's/DEFAULT_SEARCH_ENGINE\s*=\s*"crossmatch"/DEFAULT_SEARCH_ENGINE = "ncbi"/' RepeatMaskerConfig.pm \
  && sed -i 's~HMMER_DIR\s*=\s*"/usr/local/hmmer"~HMMER_DIR = "/usr/bin"~' RepeatMaskerConfig.pm \
  && sed -i "s~RMBLAST_DIR\s*=\s*\"/usr/local/rmblast\"~RMBLAST_DIR = \"${RMBLAST_PREFIX}/bin\"~" RepeatMaskerConfig.pm \
  && sed -i 's~"$REPEATMASKER_DIR/Libraries"~$ENV{'RM_LIB'}~' RepeatMaskerConfig.pm \
  && mv "${RMASK_PREFIX}/Libraries/RepeatPeps.lib" "${RMASK_PREFIX}" \
  && rm -rf -- "${RMASK_PREFIX}/Libraries" \
  && perl -i -0pe 's/^#\!.*perl.*/#\!\/usr\/bin\/env perl/g' \
       RepeatMasker \
       DateRepeats \
       ProcessRepeats \
       RepeatProteinMask \
       DupMasker \
       util/queryRepeatDatabase.pl \
       util/queryTaxonomyDatabase.pl \
       util/rmOutToGFF3.pl \
       util/rmToUCSCTables.pl


# NOTWORKING$REPEATMASKER_DIR/Libraries
# Configure repeatmodeler scripts
# Unsure if RECON_DIR needs to point to base dir or bin
WORKDIR ${RMOD_PREFIX}
RUN  cp RepModelConfig.pm.tmpl RepModelConfig.pm \
  && sed -i "s~REPEATMASKER_DIR\s*=\s*\"/usr/local/RepeatMasker\"~REPEATMASKER_DIR = \"${RMASK_PREFIX}\"~" RepModelConfig.pm \
  && sed -i 's~"$REPEATMASKER_DIR/Libraries/RepeatMasker.lib"~$ENV{"RM_LIB"} . "/RepeatMasker.lib"~' RepModelConfig.pm \
  && sed -i "s~RMBLAST_DIR\s*=\s*\"/usr/local/rmblast\"~RMBLAST_DIR = \"${RMBLAST_PREFIX}/bin\"~" RepModelConfig.pm \
  && sed -i "s~RECON_DIR\s*=\s*\"/usr/local/bin\"~RECON_DIR = \"${RECON_PREFIX}\"~" RepModelConfig.pm \
  && sed -i "s~NSEG_PRGM\s*=\s*\"/usr/local/bin/nseg\"~NSEG_PRGM = \"${NSEG_PREFIX}/nseg\"~" RepModelConfig.pm \
  && sed -i "s~RSCOUT_DIR\s*=\s*\"/usr/local/bin/\"~RSCOUT_DIR = \"${REPEATSCOUT_PREFIX}\"~" RepModelConfig.pm \
  && perl -i -0pe 's/^#\!.*/#\!\/usr\/bin\/env perl/g' \
       configure \
       BuildDatabase \
       Refiner \
       RepeatClassifier \
       RepeatModeler \
       TRFMask \
       util/dfamConsensusTool.pl \
       util/renameIds.pl \
       util/viewMSA.pl

COPY prep_repeatmasker_lib.sh /usr/local/bin
RUN chmod a+x /usr/local/bin/prep_repeatmasker_lib.sh

ENV PATH="${RMASK_PREFIX}:${RMASK_PREFIX}/util:${RMOD_PREFIX}:${RMOD_PREFIX}/util:${RMBLAST_PREFIX}/bin:${NSEG_PREFIX}:${REPEATSCOUT_PREFIX}:${RECON_PREFIX}/bin:${RECON_PREFIX}/scripts:${PATH}"
ENV CPATH="${RMBLAST_PREFIX}/include"
ENV LIBRARY_PATH="${RMBLAST_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${RMBLAST_PREFIX}/lib:${LD_LIBRARY_PATH}"
