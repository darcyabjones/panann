PWD=$(shell pwd)
DEBIAN_VERSION=buster-20191014-slim
TAG=panann-v0.0.1
MAINTAINER=darcy.ab.jones@gmail.com
DOCKER_USER=darcyabjones
BASE_IMAGE=$(DOCKER_USER)/base:$(TAG)

INSTALL_BASE=/opt

define make_image_name
$(DOCKER_USER)/$(notdir $1):$(TAG)
endef

# SOFTWARE_VERSIONs

AEGEAN_VERSION=v0.15.0
AEGEAN_URL=https://github.com/standage/AEGeAn/archive/$(AEGEAN_VERSION).tar.gz
AEGEAN_PREFIX_ARG=$(INSTALL_BASE)/aegean/$(AEGEAN_VERSION)

AUGUSTUS_COMMIT=8b1b14a7489e4545e89c8725dc33268f6c2a9117
AUGUSTUS_REPO=https://github.com/Gaius-Augustus/Augustus.git
AUGUSTUS_PREFIX_ARG=$(INSTALL_BASE)/augustus/$(AUGUSTUS_COMMIT)

BAMTOOLS_TAG=v2.5.1
BAMTOOLS_REPO=https://github.com/pezmaster31/bamtools.git
BAMTOOLS_PREFIX_ARG=$(INSTALL_BASE)/bamtools/$(BAMTOOLS_TAG)

BEDTOOLS_VERSION=2.28.0
BEDTOOLS_URL=https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz
BEDTOOLS_PREFIX_ARG=$(INSTALL_BASE)/bedtools/$(BEDTOOLS_VERSION)

BOWTIE2_TAG=v2.3.5.1
BOWTIE2_REPO=https://github.com/BenLangmead/bowtie2.git
BOWTIE2_PREFIX_ARG=$(INSTALL_BASE)/bowtie2/$(BOWTIE2_TAG)

BRAKER_COMMIT=e117150b8ad66ecf7cd5828c7f7fe476a4a8c191
BRAKER_REPO=https://github.com/Gaius-Augustus/BRAKER.git
BRAKER_PREFIX_ARG=$(INSTALL_BASE)/braker/$(BRAKER_COMMIT)

BUSCO_COMMIT=1554283ab8ee7dd5b5290f4f748234f456c36e66
BUSCO_REPO=https://gitlab.com/ezlab/busco.git
BUSCO_PREFIX_ARG=$(INSTALL_BASE)/busco/$(BUSCO_COMMIT)

CODINGQUARRY_VERSION=v2.0
CODINGQUARRY_URL=https://downloads.sourceforge.net/project/codingquarry/CodingQuarry_v2.0.tar.gz
CODINGQUARRY_PREFIX_ARG=$(INSTALL_BASE)/codingquarry/$(CODINGQUARRY_VERSION)

DEEPSIG_COMMIT=69e01cb
DEEPSIG_PREFIX_ARG=$(INSTALL_BASE)/deepsig/$(DEEPSIG_COMMIT)
DEEPSIG_REPO=https://github.com/BolognaBiocomp/deepsig.git

EVM_COMMIT=73350ce
EVM_REPO=https://github.com/EVidenceModeler/EVidenceModeler.git
EVM_PREFIX_ARG=$(INSTALL_BASE)/evm/$(EVM_COMMIT)

EXONERATE_VERSION=2.4.0
EXONERATE_URL=http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.4.0.tar.gz
EXONERATE_PREFIX_ARG=$(INSTALL_BASE)/exonerate/$(EXONERATE_VERSION)

FASTA_VERSION=36.3.8g
FASTA_URL=http://faculty.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.8g.tar.gz
FASTA_PREFIX_ARG=$(INSTALL_BASE)/fasta/$(FASTA_VERSION)

GEMOMA_VERSION=1.6.1
GEMOMA_URL=http://www.jstacs.de/downloads/GeMoMa-1.6.1.zip
GEMOMA_PREFIX_ARG=$(INSTALL_BASE)/gemoma/$(GEMOMA_VERSION)

GENEMARKES_VERSION=v.4.38
GENEMARKES_TAR=sources/gm_et_linux_64.tar.gz
GENEMARKES_KEY=sources/gm_key_64
GENEMARKES_PREFIX_ARG=$(INSTALL_BASE)/genemarkes/$(GENEMARKES_VERSION)

GENOMETOOLS_VERSION=1.5.10
GENOMETOOLS_URL=http://genometools.org/pub/genometools-$(GENOMETOOLS_VERSION).tar.gz
GENOMETOOLS_PREFIX_ARG=$(INSTALL_BASE)/genometools/$(GENOMETOOLS_VERSION)

GFFPAL_TAG=35753d8
GFFPAL_REPO=https://github.com/darcyabjones/gffpal.git
GFFPAL_PREFIX_ARG=$(INSTALL_BASE)/gffpal/$(GFFPAL_TAG)

GMAP_VERSION=2019-05-12
GMAP_URL=http://research-pub.gene.com/gmap/src/gmap-gsnap-$(GMAP_VERSION).tar.gz
GMAP_PREFIX_ARG=$(INSTALL_BASE)/gmap/$(GMAP_VERSION)

HTSLIB_TAG=1.9
BCFTOOLS_TAG=1.9
SAMTOOLS_TAG=1.9
HTSLIB_REPO=https://github.com/samtools/htslib.git
BCFTOOLS_REPO=https://github.com/samtools/bcftools.git
SAMTOOLS_REPO=https://github.com/samtools/samtools.git
HTSLIB_PREFIX_ARG=$(INSTALL_BASE)/htslib/$(HTSLIB_TAG)
BCFTOOLS_PREFIX_ARG=$(INSTALL_BASE)/bcftools/$(BCFTOOLS_TAG)
SAMTOOLS_PREFIX_ARG=$(INSTALL_BASE)/samtools/$(SAMTOOLS_TAG)

JELLYFISH_VERSION=2.2.10
JELLYFISH_URL=https://github.com/gmarcais/Jellyfish/releases/download/v$(JELLYFISH_VERSION)/jellyfish-$(JELLYFISH_VERSION).tar.gz
JELLYFISH_PREFIX_ARG=$(INSTALL_BASE)/jellyfish/$(JELLYFISH_VERSION)

KERAS_VERSION=2.2.4
KERAS_PREFIX_ARG=$(INSTALL_BASE)/keras/$(KERAS_VERSION)

MMSEQS_TAG=9-d36de
MMSEQS_REPO=https://github.com/soedinglab/MMseqs2.git
MMSEQS_PREFIX_ARG=$(INSTALL_BASE)/mmseqs/$(MMSEQS_TAG)

PASA_TAG=pasa-v2.3.3
PASA_REPO=https://github.com/PASApipeline/PASApipeline.git
PASA_PREFIX_ARG=$(INSTALL_BASE)/pasa/$(PASA_TAG)

SALMON_TAG=v0.13.1
SALMON_REPO=https://github.com/COMBINE-lab/salmon.git
SALMON_PREFIX_ARG=$(INSTALL_BASE)/salmon/$(SALMON_TAG)

SEQRENAMER_TAG=v0.0.1
SEQRENAMER_REPO=https://github.com/darcyabjones/seqrenamer.git
SEQRENAMER_PREFIX_ARG=$(INSTALL_BASE)/seqrenamer/$(SEQRENAMER_TAG)

SIGNALP5_VERSION=5.0
SIGNALP5_PREFIX_ARG=$(INSTALL_BASE)/signalp/$(SIGNALP5_VERSION)
SIGNALP5_TAR=sources/signalp-5.0.Linux.tar.gz

SPALN_TAG=Ver.2.3.3
SPALN_REPO=https://github.com/ogotoh/spaln.git
SPALN_PREFIX_ARG=$(INSTALL_BASE)/spaln/$(SPALN_TAG)

STAR_VERSION=2.7.0e
STAR_URL=https://github.com/alexdobin/STAR/archive/$(STAR_VERSION).tar.gz
STAR_PREFIX_ARG=$(INSTALL_BASE)/star/$(STAR_VERSION)

STRINGTIE_VERSION=1.3.6
STRINGTIE_URL=http://ccb.jhu.edu/software/stringtie/dl/stringtie-$(STRINGTIE_VERSION).tar.gz
STRINGTIE_PREFIX_ARG=$(INSTALL_BASE)/stringtie/$(STRINGTIE_VERSION)

TENSORFLOW_VERSION=1.5.0
TENSORFLOW_PREFIX_ARG=$(INSTALL_BASE)/tensorflow/$(TENSORFLOW_VERSION)

TRINITY_TAG=Trinity-v2.8.4
TRINITY_REPO=https://github.com/trinityrnaseq/trinityrnaseq.git
TRINITY_PREFIX_ARG=$(INSTALL_BASE)/trinity/$(TRINITY_TAG)



GENERAL=htslib python3 genometools aegean seqrenamer gffpal evm
ALIGNMENT=spaln gmap mmseqs
RNASEQ=trinity star stringtie
ANNOTATION=augustus codingquarry deepsig busco pasa gemoma braker
OPEN=$(GENERAL) $(ALIGNMENT) $(RNASEQ) $(ANNOTATION)

PROPRIETARY=genemarkes signalp5


all: docker/all singularity/all

# Docker targets
docker: docker/panann
docker/all: docker/general docker/alignment docker/rnaseq docker/annotation
	@echo "WARNING: proprietary software will not be built by default."
	@echo "WARNING: call make docker/proprietary explicitly to create them."

docker/general: $(addprefix docker/,$(GENERAL))
docker/alignment: $(addprefix docker/,$(ALIGNMENT))
docker/rnaseq: $(addprefix docker/,$(RNASEQ))
docker/annotation: $(addprefix docker/,$(ANNOTATION))

# Proprietary software won't be built by default.
docker/proprietary: $(addprefix docker/,$(PROPRIETARY))

# Singularity targets
singularity: singularity/panann.sif
singularity/all: singularity/general singularity/alignment singularity/rnaseq singularity/annotation
	@echo "WARNING: proprietary software will not be built by default."
	@echo "WARNING: call make singularity/proprietary explicitly to create them."
singularity/general: $(addprefix singularity/,$(addsuffix .sif,$(GENERAL)))
singularity/alignment: $(addprefix singularity/,$(addsuffix .sif,$(ALIGNMENT)))
singularity/rnaseq: $(addprefix singularity/,$(addsuffix .sif,$(RNASEQ)))
singularity/annotation: $(addprefix singularity/,$(addsuffix .sif,$(ANNOTATION)))

# Proprietary software won't be built by default.
singularity/proprietary: $(addprefix singularity/,$(addsuffix .sif,$(PROPRIETARY)))

# Pull containers for dockerhub
dockerpull: dockerpull/panann
dockerpull/all: dockerpull/general dockerpull/alignment dockerpull/rnaseq dockerpull/annotation
dockerpull/general: $(addprefix dockerpull/,$(GENERAL))
dockerpull/alignment: $(addprefix dockerpull/,$(ALIGNMENT))
dockerpull/rnaseq: $(addprefix dockerpull/,$(RNASEQ))
dockerpull/annotation: $(addprefix dockerpull/,$(ANNOTATION))


# Utils
clean:
	docker stop $(shell docker ps -aq) || true
	docker rm $(shell docker ps -a -q) || true
	docker rmi $(shell docker images -q) || true
	singularity cache clean --force
	rm -rf -- docker

tidy:
	# Just removes untagged images to free up space
	singularity cache clean --force
	docker image prune -f

# Generic rules
singularity/%.sif: docker/%
	@mkdir -p "$(dir $@)"
	@mkdir -p tmp tmp/cache
	SINGULARITY_TMPDIR="${PWD}/tmp" SINGULARITY_CACHEDIR="${PWD}/tmp/cache" singularity build \
	  --force \
	  --tmpdir "${PWD}/tmp" \
	  "$@" \
	  "docker-daemon://$(call make_image_name,$(notdir $(basename $@)))"

docker/%.tar.gz: docker/%
	@mkdir -p "$(dir $@)"
	docker save "$(call make_image_name,$(basename $(basename $@)))" | gzip > $@

dockerpull/%:
	docker pull "$(call make_image_name,$(notdir $@))"
	@mkdir -p "$(dir $@)"
	@touch docker/$(notdir $@)


docker/base: base.Dockerfile
	docker build \
	  --build-arg DEBIAN_VERSION=$(DEBIAN_VERSION) \
	  --build-arg MAINTAINER="$(MAINTAINER)" \
	  --build-arg PIPELINE_VERSION="$(TAG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


# Basics


docker/python3: python3.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


#

docker/aegean: aegean.Dockerfile docker/base docker/genometools
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg GENOMETOOLS_IMAGE="$(call make_image_name,genometools)" \
	  --build-arg GENOMETOOLS_VERSION="$(GENOMETOOLS_VERSION)" \
	  --build-arg GENOMETOOLS_PREFIX_ARG="$(GENOMETOOLS_PREFIX_ARG)" \
	  --build-arg AEGEAN_VERSION="$(AEGEAN_VERSION)" \
	  --build-arg AEGEAN_URL="$(AEGEAN_URL)" \
	  --build-arg AEGEAN_PREFIX_ARG="$(AEGEAN_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/augustus: augustus.Dockerfile docker/base docker/htslib
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg AUGUSTUS_COMMIT="$(AUGUSTUS_COMMIT)" \
	  --build-arg AUGUSTUS_REPO="$(AUGUSTUS_REPO)" \
	  --build-arg AUGUSTUS_PREFIX_ARG="$(AUGUSTUS_PREFIX_ARG)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_REPO="$(HTSLIB_REPO)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_REPO="$(SAMTOOLS_REPO)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/bamtools: bamtools.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg BAMTOOLS_TAG="$(BAMTOOLS_TAG)" \
	  --build-arg BAMTOOLS_REPO="$(BAMTOOLS_REPO)" \
	  --build-arg BAMTOOLS_PREFIX_ARG="$(BAMTOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/bedtools: bedtools.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg BEDTOOLS_VERSION="$(BEDTOOLS_VERSION)" \
	  --build-arg BEDTOOLS_URL="$(BEDTOOLS_URL)" \
	  --build-arg BEDTOOLS_PREFIX_ARG="$(BEDTOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/bowtie2: bowtie2.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg BOWTIE2_TAG="$(BOWTIE2_TAG)" \
	  --build-arg BOWTIE2_REPO="$(BOWTIE2_REPO)" \
	  --build-arg BOWTIE2_PREFIX_ARG="$(BOWTIE2_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/braker: braker.Dockerfile docker/base docker/augustus docker/htslib docker/bamtools docker/spaln
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg BRAKER_COMMIT="$(BRAKER_COMMIT)" \
	  --build-arg BRAKER_REPO="$(BRAKER_REPO)" \
	  --build-arg BRAKER_PREFIX_ARG="$(BRAKER_PREFIX_ARG)" \
	  --build-arg AUGUSTUS_IMAGE="$(call make_image_name,augustus)" \
	  --build-arg AUGUSTUS_COMMIT="$(AUGUSTUS_COMMIT)" \
	  --build-arg AUGUSTUS_PREFIX_ARG="$(AUGUSTUS_PREFIX_ARG)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --build-arg BAMTOOLS_IMAGE="$(call make_image_name,bamtools)" \
	  --build-arg BAMTOOLS_TAG="$(BAMTOOLS_TAG)" \
	  --build-arg BAMTOOLS_PREFIX_ARG="$(BAMTOOLS_PREFIX_ARG)" \
	  --build-arg SPALN_IMAGE="$(call make_image_name,spaln)" \
	  --build-arg SPALN_TAG="$(SPALN_TAG)" \
	  --build-arg SPALN_PREFIX_ARG="$(SPALN_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/busco: busco.Dockerfile docker/base docker/augustus
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg BUSCO_COMMIT="$(BUSCO_COMMIT)" \
	  --build-arg BUSCO_REPO="$(BUSCO_REPO)" \
	  --build-arg BUSCO_PREFIX_ARG="$(BUSCO_PREFIX_ARG)" \
	  --build-arg AUGUSTUS_IMAGE="$(call make_image_name,augustus)" \
	  --build-arg AUGUSTUS_COMMIT="$(AUGUSTUS_COMMIT)" \
	  --build-arg AUGUSTUS_PREFIX_ARG="$(AUGUSTUS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/codingquarry: codingquarry.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg CODINGQUARRY_VERSION="$(CODINGQUARRY_VERSION)" \
	  --build-arg CODINGQUARRY_URL="$(CODINGQUARRY_URL)" \
	  --build-arg CODINGQUARRY_PREFIX_ARG="$(CODINGQUARRY_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/deepsig: deepsig.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg DEEPSIG_COMMIT="$(DEEPSIG_COMMIT)" \
	  --build-arg DEEPSIG_REPO="$(DEEPSIG_REPO)" \
	  --build-arg DEEPSIG_PREFIX_ARG="$(DEEPSIG_PREFIX_ARG)" \
	  --build-arg TENSORFLOW_VERSION="$(TENSORFLOW_VERSION)" \
	  --build-arg TENSORFLOW_PREFIX_ARG="$(TENSORFLOW_PREFIX_ARG)" \
	  --build-arg KERAS_VERSION="$(KERAS_VERSION)" \
	  --build-arg KERAS_PREFIX_ARG="$(KERAS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/evm: evm.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg EVM_COMMIT="$(EVM_COMMIT)" \
	  --build-arg EVM_REPO="$(EVM_REPO)" \
	  --build-arg EVM_PREFIX_ARG="$(EVM_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/exonerate: exonerate.Dockerfile docker/base docker/htslib
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg EXONERATE_VERSION="$(EXONERATE_VERSION)" \
	  --build-arg EXONERATE_URL="$(EXONERATE_URL)" \
	  --build-arg EXONERATE_PREFIX_ARG="$(EXONERATE_PREFIX_ARG)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/fasta: fasta.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg FASTA_VERSION="$(FASTA_VERSION)" \
	  --build-arg FASTA_URL="$(FASTA_URL)" \
	  --build-arg FASTA_PREFIX_ARG="$(FASTA_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/gemoma: gemoma.Dockerfile docker/base docker/mmseqs docker/htslib
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg GEMOMA_VERSION="$(GEMOMA_VERSION)" \
	  --build-arg GEMOMA_URL="$(GEMOMA_URL)" \
	  --build-arg GEMOMA_PREFIX_ARG="$(GEMOMA_PREFIX_ARG)" \
	  --build-arg MMSEQS_IMAGE="$(call make_image_name,mmseqs)" \
	  --build-arg MMSEQS_TAG="$(MMSEQS_TAG)" \
	  --build-arg MMSEQS_PREFIX_ARG="$(MMSEQS_PREFIX_ARG)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/genemarkes: genemarkes.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg GENEMARKES_VERSION="$(GENEMARKES_VERSION)" \
	  --build-arg GENEMARKES_TAR="$(GENEMARKES_TAR)" \
	  --build-arg GENEMARKES_KEY="$(GENEMARKES_KEY)" \
	  --build-arg GENEMARKES_PREFIX_ARG="$(GENEMARKES_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/genometools: genometools.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg GENOMETOOLS_VERSION="$(GENOMETOOLS_VERSION)" \
	  --build-arg GENOMETOOLS_URL="$(GENOMETOOLS_URL)" \
	  --build-arg GENOMETOOLS_PREFIX_ARG="$(GENOMETOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/gffpal: gffpal.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg GFFPAL_TAG="$(GFFPAL_TAG)" \
	  --build-arg GFFPAL_REPO="$(GFFPAL_REPO)" \
	  --build-arg GFFPAL_PREFIX_ARG="$(GFFPAL_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/gmap: gmap.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg GMAP_VERSION="$(GMAP_VERSION)" \
	  --build-arg GMAP_URL="$(GMAP_URL)" \
	  --build-arg GMAP_PREFIX_ARG="$(GMAP_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/htslib: htslib.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_REPO="$(HTSLIB_REPO)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_REPO="$(SAMTOOLS_REPO)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --build-arg BCFTOOLS_TAG="$(BCFTOOLS_TAG)" \
	  --build-arg BCFTOOLS_REPO="$(BCFTOOLS_REPO)" \
	  --build-arg BCFTOOLS_PREFIX_ARG="$(BCFTOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/jellyfish: jellyfish.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg JELLYFISH_VERSION="$(JELLYFISH_VERSION)" \
	  --build-arg JELLYFISH_URL="$(JELLYFISH_URL)" \
	  --build-arg JELLYFISH_PREFIX_ARG="$(JELLYFISH_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/mmseqs: mmseqs.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg MMSEQS_TAG="$(MMSEQS_TAG)" \
	  --build-arg MMSEQS_REPO="$(MMSEQS_REPO)" \
	  --build-arg MMSEQS_PREFIX_ARG="$(MMSEQS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/pasa: pasa.Dockerfile docker/base docker/gmap docker/fasta docker/htslib
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg PASA_TAG="$(PASA_TAG)" \
	  --build-arg PASA_REPO="$(PASA_REPO)" \
	  --build-arg PASA_PREFIX_ARG="$(PASA_PREFIX_ARG)" \
	  --build-arg GMAP_IMAGE="$(call make_image_name,gmap)" \
	  --build-arg GMAP_VERSION="$(GMAP_VERSION)" \
	  --build-arg GMAP_PREFIX_ARG="$(GMAP_PREFIX_ARG)" \
	  --build-arg FASTA_IMAGE="$(call make_image_name,fasta)" \
	  --build-arg FASTA_VERSION="$(FASTA_VERSION)" \
	  --build-arg FASTA_PREFIX_ARG="$(FASTA_PREFIX_ARG)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/salmon: salmon.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg SALMON_TAG="$(SALMON_TAG)" \
	  --build-arg SALMON_REPO="$(SALMON_REPO)" \
	  --build-arg SALMON_PREFIX_ARG="$(SALMON_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/seqrenamer: seqrenamer.Dockerfile docker/base docker/gffpal
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg SEQRENAMER_TAG="$(SEQRENAMER_TAG)" \
	  --build-arg SEQRENAMER_REPO="$(SEQRENAMER_REPO)" \
	  --build-arg SEQRENAMER_PREFIX_ARG="$(SEQRENAMER_PREFIX_ARG)" \
	  --build-arg GFFPAL_IMAGE="$(call make_image_name,gffpal)" \
	  --build-arg GFFPAL_TAG="$(GFFPAL_TAG)" \
	  --build-arg GFFPAL_PREFIX_ARG="$(GFFPAL_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/signalp5: signalp5.Dockerfile docker/base $(SIGNALP5_TAR)
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg SIGNALP5_VERSION="$(SIGNALP5_VERSION)" \
	  --build-arg SIGNALP5_TAR="$(SIGNALP5_TAR)" \
	  --build-arg SIGNALP5_PREFIX_ARG="$(SIGNALP5_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/spaln: spaln.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg SPALN_TAG="$(SPALN_TAG)" \
	  --build-arg SPALN_REPO="$(SPALN_REPO)" \
	  --build-arg SPALN_PREFIX_ARG="$(SPALN_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/star: star.Dockerfile docker/base docker/htslib
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg STAR_VERSION="$(STAR_VERSION)" \
	  --build-arg STAR_URL="$(STAR_URL)" \
	  --build-arg STAR_PREFIX_ARG="$(STAR_PREFIX_ARG)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/stringtie: stringtie.Dockerfile docker/base docker/htslib
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg STRINGTIE_VERSION="$(STRINGTIE_VERSION)" \
	  --build-arg STRINGTIE_URL="$(STRINGTIE_URL)" \
	  --build-arg STRINGTIE_PREFIX_ARG="$(STRINGTIE_PREFIX_ARG)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@

docker/trinity: docker/jellyfish docker/bowtie2 docker/htslib docker/salmon
docker/trinity: trinity.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg TRINITY_TAG="$(TRINITY_TAG)" \
	  --build-arg TRINITY_REPO="$(TRINITY_REPO)" \
	  --build-arg TRINITY_PREFIX_ARG="$(TRINITY_PREFIX_ARG)" \
	  --build-arg JELLYFISH_IMAGE="$(call make_image_name,jellyfish)" \
	  --build-arg JELLYFISH_VERSION="$(JELLYFISH_VERSION)" \
	  --build-arg JELLYFISH_PREFIX_ARG="$(JELLYFISH_PREFIX_ARG)" \
	  --build-arg BOWTIE2_IMAGE="$(call make_image_name,bowtie2)" \
	  --build-arg BOWTIE2_TAG="$(BOWTIE2_TAG)" \
	  --build-arg BOWTIE2_PREFIX_ARG="$(BOWTIE2_PREFIX_ARG)" \
	  --build-arg SALMON_IMAGE="$(call make_image_name,salmon)" \
	  --build-arg SALMON_TAG="$(SALMON_TAG)" \
	  --build-arg SALMON_PREFIX_ARG="$(SALMON_PREFIX_ARG)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


## The monolith!!!!
docker/panann: docker/aegean docker/augustus docker/bamtools
docker/panann: docker/bedtools docker/bowtie2 docker/braker
docker/panann: docker/busco docker/codingquarry docker/deepsig
docker/panann: docker/evm docker/exonerate docker/fasta
docker/panann: docker/gemoma docker/genometools docker/gffpal
docker/panann: docker/gmap docker/htslib docker/jellyfish
docker/panann: docker/mmseqs docker/pasa docker/python3 docker/salmon
docker/panann: docker/spaln docker/seqrenamer docker/star
docker/panann: docker/stringtie docker/trinity
docker/panann: panann.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg AEGEAN_IMAGE="$(call make_image_name,aegean)" \
	  --build-arg AEGEAN_VERSION="$(AEGEAN_VERSION)" \
	  --build-arg AEGEAN_PREFIX_ARG="$(AEGEAN_PREFIX_ARG)" \
	  --build-arg AUGUSTUS_IMAGE="$(call make_image_name,augustus)" \
	  --build-arg AUGUSTUS_COMMIT="$(AUGUSTUS_COMMIT)" \
	  --build-arg AUGUSTUS_PREFIX_ARG="$(AUGUSTUS_PREFIX_ARG)" \
	  --build-arg BAMTOOLS_IMAGE="$(call make_image_name,bamtools)" \
	  --build-arg BAMTOOLS_TAG="$(BAMTOOLS_TAG)" \
	  --build-arg BAMTOOLS_PREFIX_ARG="$(BAMTOOLS_PREFIX_ARG)" \
	  --build-arg BEDTOOLS_IMAGE="$(call make_image_name,bedtools)" \
	  --build-arg BEDTOOLS_VERSION="$(BEDTOOLS_VERSION)" \
	  --build-arg BEDTOOLS_PREFIX_ARG="$(BEDTOOLS_PREFIX_ARG)" \
	  --build-arg BOWTIE2_IMAGE="$(call make_image_name,bowtie2)" \
	  --build-arg BOWTIE2_TAG="$(BOWTIE2_TAG)" \
	  --build-arg BOWTIE2_PREFIX_ARG="$(BOWTIE2_PREFIX_ARG)" \
	  --build-arg BRAKER_IMAGE="$(call make_image_name,braker)" \
	  --build-arg BRAKER_COMMIT="$(BRAKER_COMMIT)" \
	  --build-arg BRAKER_PREFIX_ARG="$(BRAKER_PREFIX_ARG)" \
	  --build-arg BUSCO_IMAGE="$(call make_image_name,busco)" \
	  --build-arg BUSCO_COMMIT="$(BUSCO_COMMIT)" \
	  --build-arg BUSCO_PREFIX_ARG="$(BUSCO_PREFIX_ARG)" \
	  --build-arg CODINGQUARRY_IMAGE="$(call make_image_name,codingquarry)" \
	  --build-arg CODINGQUARRY_VERSION="$(CODINGQUARRY_VERSION)" \
	  --build-arg CODINGQUARRY_PREFIX_ARG="$(CODINGQUARRY_PREFIX_ARG)" \
	  --build-arg DEEPSIG_IMAGE="$(call make_image_name,deepsig)" \
	  --build-arg DEEPSIG_COMMIT="$(DEEPSIG_COMMIT)" \
	  --build-arg DEEPSIG_PREFIX_ARG="$(DEEPSIG_PREFIX_ARG)" \
	  --build-arg TENSORFLOW_VERSION="$(TENSORFLOW_VERSION)" \
	  --build-arg TENSORFLOW_PREFIX_ARG="$(TENSORFLOW_PREFIX_ARG)" \
	  --build-arg KERAS_VERSION="$(KERAS_VERSION)" \
	  --build-arg KERAS_PREFIX_ARG="$(KERAS_PREFIX_ARG)" \
	  --build-arg EVM_IMAGE="$(call make_image_name,evm)" \
	  --build-arg EVM_COMMIT="$(EVM_COMMIT)" \
	  --build-arg EVM_PREFIX_ARG="$(EVM_PREFIX_ARG)" \
	  --build-arg EXONERATE_IMAGE="$(call make_image_name,exonerate)" \
	  --build-arg EXONERATE_VERSION="$(EXONERATE_VERSION)" \
	  --build-arg EXONERATE_PREFIX_ARG="$(EXONERATE_PREFIX_ARG)" \
	  --build-arg FASTA_IMAGE="$(call make_image_name,fasta)" \
	  --build-arg FASTA_VERSION="$(FASTA_VERSION)" \
	  --build-arg FASTA_PREFIX_ARG="$(FASTA_PREFIX_ARG)" \
	  --build-arg GEMOMA_IMAGE="$(call make_image_name,gemoma)" \
	  --build-arg GEMOMA_VERSION="$(GEMOMA_VERSION)" \
	  --build-arg GEMOMA_PREFIX_ARG="$(GEMOMA_PREFIX_ARG)" \
	  --build-arg GENOMETOOLS_IMAGE="$(call make_image_name,genometools)" \
	  --build-arg GENOMETOOLS_VERSION="$(GENOMETOOLS_VERSION)" \
	  --build-arg GENOMETOOLS_PREFIX_ARG="$(GENOMETOOLS_PREFIX_ARG)" \
	  --build-arg GFFPAL_IMAGE="$(call make_image_name,gffpal)" \
	  --build-arg GFFPAL_TAG="$(GFFPAL_TAG)" \
	  --build-arg GFFPAL_PREFIX_ARG="$(GFFPAL_PREFIX_ARG)" \
	  --build-arg GMAP_IMAGE="$(call make_image_name,gmap)" \
	  --build-arg GMAP_VERSION="$(GMAP_VERSION)" \
	  --build-arg GMAP_PREFIX_ARG="$(GMAP_PREFIX_ARG)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --build-arg BCFTOOLS_TAG="$(BCFTOOLS_TAG)" \
	  --build-arg BCFTOOLS_PREFIX_ARG="$(BCFTOOLS_PREFIX_ARG)" \
	  --build-arg JELLYFISH_IMAGE="$(call make_image_name,jellyfish)" \
	  --build-arg JELLYFISH_VERSION="$(JELLYFISH_VERSION)" \
	  --build-arg JELLYFISH_PREFIX_ARG="$(JELLYFISH_PREFIX_ARG)" \
	  --build-arg MMSEQS_IMAGE="$(call make_image_name,mmseqs)" \
	  --build-arg MMSEQS_TAG="$(MMSEQS_TAG)" \
	  --build-arg MMSEQS_PREFIX_ARG="$(MMSEQS_PREFIX_ARG)" \
	  --build-arg PASA_IMAGE="$(call make_image_name,pasa)" \
	  --build-arg PASA_TAG="$(PASA_TAG)" \
	  --build-arg PASA_PREFIX_ARG="$(PASA_PREFIX_ARG)" \
	  --build-arg PYTHON3_IMAGE="$(call make_image_name,python3)" \
	  --build-arg SALMON_IMAGE="$(call make_image_name,salmon)" \
	  --build-arg SALMON_TAG="$(SALMON_TAG)" \
	  --build-arg SALMON_PREFIX_ARG="$(SALMON_PREFIX_ARG)" \
	  --build-arg SEQRENAMER_IMAGE="$(call make_image_name,seqrenamer)" \
	  --build-arg SEQRENAMER_TAG="$(SEQRENAMER_TAG)" \
	  --build-arg SEQRENAMER_PREFIX_ARG="$(SEQRENAMER_PREFIX_ARG)" \
	  --build-arg SPALN_IMAGE="$(call make_image_name,spaln)" \
	  --build-arg SPALN_TAG="$(SPALN_TAG)" \
	  --build-arg SPALN_PREFIX_ARG="$(SPALN_PREFIX_ARG)" \
	  --build-arg STAR_IMAGE="$(call make_image_name,star)" \
	  --build-arg STAR_VERSION="$(STAR_VERSION)" \
	  --build-arg STAR_PREFIX_ARG="$(STAR_PREFIX_ARG)" \
	  --build-arg STRINGTIE_IMAGE="$(call make_image_name,stringtie)" \
	  --build-arg STRINGTIE_VERSION="$(STRINGTIE_VERSION)" \
	  --build-arg STRINGTIE_PREFIX_ARG="$(STRINGTIE_PREFIX_ARG)" \
	  --build-arg TRINITY_IMAGE="$(call make_image_name,trinity)" \
	  --build-arg TRINITY_TAG="$(TRINITY_TAG)" \
	  --build-arg TRINITY_PREFIX_ARG="$(TRINITY_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/panann-plus: panann-plus.Dockerfile docker/panann docker/genemarkes docker/signalp5
	docker build \
	  --build-arg IMAGE="$(call make_image_name,panann)" \
	  --build-arg SIGNALP5_IMAGE="$(call make_image_name,signalp5)" \
	  --build-arg SIGNALP5_VERSION="$(SIGNALP5_VERSION)" \
	  --build-arg SIGNALP5_PREFIX_ARG="$(SIGNALP5_PREFIX_ARG)" \
	  --build-arg GENEMARKES_IMAGE="$(call make_image_name,genemarkes)" \
	  --build-arg GENEMARKES_VERSION="$(GENEMARKES_VERSION)" \
	  --build-arg GENEMARKES_PREFIX_ARG="$(GENEMARKES_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@





# SIBELIAZ_COMMIT=c9e2257
# SIBELIAZ_REPO=https://github.com/medvedevgroup/SibeliaZ.git
# SIBELIAZ_PREFIX_ARG=$(INSTALL_BASE)/sibeliaz/$(SIBELIAZ_COMMIT)
#
# SKMER_COMMIT=d6f7c23
# SKMER_REPO=https://github.com/shahab-sarmashghi/Skmer.git
# SKMER_PREFIX_ARG=$(INSTALL_BASE)/skmer/$(SKMER_COMMIT)
# SKMER_PREFIX=$(SKMER_PREFIX_ARG)
#
# MASH_VERSION=v2.2
# MASH_URL=https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar

#MINIMAP2_TAG=v2.16
#K8_VERSION=0.2.4
#MINIMAP2_REPO=https://github.com/lh3/minimap2.git
#K8_URL=https://github.com/attractivechaos/k8/releases/download/v$(K8_VERSION)/k8-$(K8_VERSION).tar.bz2
#MINIMAP2_PREFIX_ARG=$(INSTALL_BASE)/minimap2/$(MINIMAP_TAG)
#K8_PREFIX_ARG=$(INSTALL_BASE)/k8/$(K8_VERSION)
#
#
#SEQTK_COMMIT=ca4785c
#SEQTK_REPO=https://github.com/lh3/seqtk.git
#SEQTK_PREFIX_ARG=$(INSTALL_BASE)/seqtk/$(SEQTK_COMMIT)
#
#docker/r: r.Dockerfile docker/base
#	docker build \
#	  --build-arg IMAGE="$(BASE_IMAGE)" \
#	  --file $< \
#	  -t "$(call make_image_name,$@)" \
#	  .
#	@mkdir -p "$(dir $@)"
#	@touch $@
#
#
#docker/sibeliaz: sibeliaz.Dockerfile docker/base
#	docker build \
#	  --build-arg IMAGE="$(BASE_IMAGE)" \
#	  --build-arg SIBELIAZ_COMMIT="$(SIBELIAZ_COMMIT)" \
#	  --build-arg SIBELIAZ_REPO="$(SIBELIAZ_REPO)" \
#	  --build-arg SIBELIAZ_PREFIX_ARG="$(SIBELIAZ_PREFIX_ARG)" \
#	  --file $< \
#	  -t "$(call make_image_name,$@)" \
#	  .
#	@mkdir -p "$(dir $@)"
#	@touch $@
#
#
#
#
#docker/skmer: skmer.Dockerfile docker/base docker/jellyfish docker/seqtk
#	docker build \
#	  --build-arg IMAGE="$(BASE_IMAGE)" \
#	  --build-arg SKMER_COMMIT="$(SKMER_COMMIT)" \
#	  --build-arg SKMER_REPO="$(SKMER_REPO)" \
#	  --build-arg SKMER_PREFIX_ARG="$(SKMER_PREFIX_ARG)" \
#	  --build-arg MASH_VERSION="$(MASH_VERSION)" \
#	  --build-arg MASH_URL="$(MASH_URL)" \
#	  --build-arg JELLYFISH_IMAGE="$(call make_image_name,jellyfish)" \
#	  --build-arg JELLYFISH_VERSION="$(JELLYFISH_VERSION)" \
#	  --build-arg JELLYFISH_PREFIX_ARG="$(JELLYFISH_PREFIX_ARG)" \
#	  --build-arg SEQTK_IMAGE="$(call make_image_name,seqtk)" \
#	  --build-arg SEQTK_COMMIT="$(SEQTK_COMMIT)" \
#	  --build-arg SEQTK_PREFIX_ARG="$(SEQTK_PREFIX_ARG)" \
#	  --file $< \
#	  -t "$(call make_image_name,$@)" \
#	  .
#	@mkdir -p "$(dir $@)"
#	@touch $@
#
#
#docker/minimap2: minimap2.Dockerfile docker/base
#	docker build \
#	  --build-arg IMAGE="$(BASE_IMAGE)" \
#	  --build-arg MINIMAP2_TAG="$(MINIMAP2_TAG)" \
#	  --build-arg MINIMAP2_REPO="$(MINIMAP2_REPO)" \
#	  --build-arg MINIMAP2_PREFIX_ARG="$(MINIMAP2_PREFIX_ARG)" \
#	  --build-arg K8_VERSION="$(K8_VERSION)" \
#	  --build-arg K8_URL="$(K8_URL)" \
#	  --build-arg K8_PREFIX_ARG="$(K8_PREFIX_ARG)" \
#	  --file $< \
#	  -t "$(call make_image_name,$@)" \
#	  .
#	@mkdir -p "$(dir $@)"
#	@touch $@
#
#
#docker/seqtk: seqtk.Dockerfile docker/base
#	docker build \
#	  --build-arg IMAGE="$(BASE_IMAGE)" \
#	  --build-arg SEQTK_COMMIT="$(SEQTK_COMMIT)" \
#	  --build-arg SEQTK_REPO="$(SEQTK_REPO)" \
#	  --build-arg SEQTK_PREFIX_ARG="$(SEQTK_PREFIX_ARG)" \
#	  --file $< \
#	  -t "$(call make_image_name,$@)" \
#	  .
#	@mkdir -p "$(dir $@)"
#	@touch $@
