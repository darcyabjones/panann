PWD=$(shell pwd)


docker: .repeatmasker .minimap .red .star .spaln .velvet .mitefinder

.rmblast: rmblast.Dockerfile
	docker build --file $< -t darcyabjones/rmblast .
	touch $@

.repeatmasker: repeatmasker.Dockerfile .rmblast
	docker build --file $< -t darcyabjones/repeatmasker .
	touch $@

.red: red.Dockerfile
	docker build --file $< -t darcyabjones/red .
	touch $@

.minimap: minimap2.Dockerfile
	docker build --file $< -t darcyabjones/minimap2 .
	touch $@

.spaln: spaln.Dockerfile
	docker build --file $< -t darcyabjones/spaln .
	touch $@

.star: star.Dockerfile
	docker build --file $< -t darcyabjones/star .
	touch $@

.velvet: velvet.Dockerfile
	docker build --file $< -t darcyabjones/velvet .
	touch $@

.jellyfish: jellyfish.Dockerfile
	docker build --file $< -t darcyabjones/jellyfish .
	touch $@

.mitefinder: mitefinder.Dockerfile
	docker build --file $< -t darcyabjones/mitefinder .
	touch $@

.bwa: bwa.Dockerfile
	docker build --file $< -t darcyabjones/bwa .
	touch $@

.htslib: htslib.Dockerfile
	docker build --file $< -t darcyabjones/htslib .
	touch $@

downloads/RepeatMaskerMetaData-20181026.tar.gz:
	mkdir -p downloads
	cd downloads \
	&& wget RepeatMaskerMetaData-20181026.tar.gz
	
downloads/Dfam.hmm.gz:
	mkdir -p downloads
	cd downloads \
	&& wget http://dfam.org/releases/Dfam_3.0/families/Dfam.hmm.gz

downloads/Dfam.embl.gz:
	mkdir -p downloads
	cd downloads \
	&& wget http://dfam.org/releases/Dfam_3.0/families/Dfam.embl.gz

downloads/RepBaseRepeatMaskerEdition-20181026.tar.gz:
	mkdir -p downloads
	cd downloads \
	&& wget \
	     --user $RB_USERNAME \
	     --password $RB_PASSWORD \
	     -O repeatmaskerlibraries.tar.gz \
         http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20181026.tar.gz


rmlib: downloads/RepeatMaskerMetaData-20181026.tar.gz downloads/Dfam.hmm.gz downloads/Dfam.embl.gz downloads/RepBaseRepeatMaskerEdition-20181026.tar.gz
	docker run \
	  --rm \
	  -v "$(PWD)/downloads:/downloads:rw" \
	  -v "$(PWD)/rmlib:/opt/rmlib:rw" \
	  -w /downloads \
	  darcyabjones/repeatmasker \
	  prep_repeatmasker_lib.sh \
	    RepBaseRepeatMaskerEdition-20181026.tar.gz \
	    RepeatMaskerMetaData-20181026.tar.gz \
	    Dfam.hmm.gz \
	    Dfam.embl.gz \
	    /opt/rmlib
