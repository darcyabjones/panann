PWD=$(shell pwd)
DEBIAN_VERSION=stretch-20190228-slim
TAG=panann-v0.0.1
DOCKER_USER=darcyabjones
BASE_IMAGE=$(DOCKER_USER)/base:$(TAG)

INSTALL_BASE=/opt

define make_image_name
$(DOCKER_USER)/$(notdir $1):$(TAG)
endef

# SOFTWARE_VERSIONs
JELLYFISH_VERSION=2.2.10
JELLYFISH_URL=https://github.com/gmarcais/Jellyfish/releases/download/v$(JELLYFISH_VERSION)/jellyfish-$(JELLYFISH_VERSION).tar.gz
JELLYFISH_PREFIX_ARG=$(INSTALL_BASE)/jellyfish/$(JELLYFISH_VERSION)

BOWTIE2_TAG=v2.3.5.1
BOWTIE2_REPO=https://github.com/BenLangmead/bowtie2.git
BOWTIE2_PREFIX_ARG=$(INSTALL_BASE)/bowtie2/$(BOWTIE2_TAG)

SALMON_TAG=v0.13.1
SALMON_REPO=https://github.com/COMBINE-lab/salmon.git
SALMON_PREFIX_ARG=$(INSTALL_BASE)/salmon/$(SALMON_TAG)

TRINITY_TAG=Trinity-v2.8.4
TRINITY_REPO=https://github.com/trinityrnaseq/trinityrnaseq.git
TRINITY_PREFIX_ARG=$(INSTALL_BASE)/trinity/$(TRINITY_TAG)

STAR_VERSION=2.7.0e
STAR_URL=https://github.com/alexdobin/STAR/archive/$(STAR_VERSION).tar.gz
STAR_PREFIX_ARG=$(INSTALL_BASE)/star/$(STAR_VERSION)

STRINGTIE_VERSION=1.3.6
STRINGTIE_URL=http://ccb.jhu.edu/software/stringtie/dl/stringtie-$(STRINGTIE_VERSION).tar.gz
STRINGTIE_PREFIX_ARG=$(INSTALL_BASE)/stringtie/$(STRINGTIE_VERSION)

HTSLIB_TAG=1.9
BCFTOOLS_TAG=1.9
SAMTOOLS_TAG=1.9
HTSLIB_REPO=https://github.com/samtools/htslib.git
BCFTOOLS_REPO=https://github.com/samtools/bcftools.git
SAMTOOLS_REPO=https://github.com/samtools/samtools.git
HTSLIB_PREFIX_ARG=$(INSTALL_BASE)/htslib/$(HTSLIB_TAG)
BCFTOOLS_PREFIX_ARG=$(INSTALL_BASE)/bcftools/$(BCFTOOLS_TAG)
SAMTOOLS_PREFIX_ARG=$(INSTALL_BASE)/samtools/$(SAMTOOLS_TAG)

BAMTOOLS_TAG=v2.5.1
BAMTOOLS_REPO=https://github.com/pezmaster31/bamtools.git
BAMTOOLS_PREFIX_ARG=$(INSTALL_BASE)/bamtools/$(BAMTOOLS_TAG)


MINIMAP2_TAG=v2.16
K8_VERSION=0.2.4
MINIMAP2_REPO=https://github.com/lh3/minimap2.git
K8_URL=https://github.com/attractivechaos/k8/releases/download/v$(K8_VERSION)/k8-$(K8_VERSION).tar.bz2
MINIMAP2_PREFIX_ARG=$(INSTALL_BASE)/minimap2/$(MINIMAP_TAG)
K8_PREFIX_ARG=$(INSTALL_BASE)/k8/$(K8_VERSION)

SPALN_TAG=Ver.2.3.3
SPALN_REPO=https://github.com/ogotoh/spaln.git
SPALN_PREFIX_ARG=$(INSTALL_BASE)/spaln/$(SPALN_TAG)

GMAP_VERSION=2019-05-12
GMAP_URL=http://research-pub.gene.com/gmap/src/gmap-gsnap-$(GMAP_VERSION).tar.gz
GMAP_PREFIX_ARG=$(INSTALL_BASE)/gmap/$(GMAP_VERSION)

FASTA_VERSION=36.3.8g
FASTA_URL=http://faculty.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.8g.tar.gz
FASTA_PREFIX_ARG=$(INSTALL_BASE)/fasta/$(FASTA_VERSION)

MMSEQS_TAG=7-4e23d
MMSEQS_REPO=https://github.com/soedinglab/MMseqs2.git
MMSEQS_PREFIX_ARG=$(INSTALL_BASE)/mmseqs/$(MMSEQS_TAG)


GENOMETOOLS_VERSION=1.5.10
GENOMETOOLS_URL=http://genometools.org/pub/genometools-$(GENOMETOOLS_VERSION).tar.gz
GENOMETOOLS_PREFIX_ARG=$(INSTALL_BASE)/genometools/$(GENOMETOOLS_VERSION)


AUGUSTUS_COMMIT=8b1b14a7489e4545e89c8725dc33268f6c2a9117
AUGUSTUS_REPO=https://github.com/Gaius-Augustus/Augustus.git
AUGUSTUS_PREFIX_ARG=$(INSTALL_BASE)/augustus/$(AUGUSTUS_COMMIT)

CODINGQUARRY_VERSION=v2.0
CODINGQUARRY_URL=https://downloads.sourceforge.net/project/codingquarry/CodingQuarry_v2.0.tar.gz
CODINGQUARRY_PREFIX_ARG=$(INSTALL_BASE)/codingquarry/$(CODINGQUARRY_VERSION)

BUSCO_COMMIT=1554283ab8ee7dd5b5290f4f748234f456c36e66
BUSCO_REPO=https://gitlab.com/ezlab/busco.git
BUSCO_PREFIX_ARG=$(INSTALL_BASE)/busco/$(BUSCO_COMMIT)

PASA_TAG=pasa-v2.3.3
PASA_REPO=https://github.com/PASApipeline/PASApipeline.git
PASA_PREFIX_ARG=$(INSTALL_BASE)/pasa/$(PASA_TAG)

BRAKER_COMMIT=e117150b8ad66ecf7cd5828c7f7fe476a4a8c191
BRAKER_REPO=https://github.com/Gaius-Augustus/BRAKER.git
BRAKER_PREFIX_ARG=$(INSTALL_BASE)/braker/$(BRAKER_COMMIT)

GEMOMA_VERSION=1.6.1
GEMOMA_URL=http://www.jstacs.de/downloads/GeMoMa-1.6.1.zip
GEMOMA_PREFIX_ARG=$(INSTALL_BASE)/gemoma/$(GEMOMA_VERSION)


GENEMARKES_VERSION=v.4.38
GENEMARKES_TAR=sources/gm_et_linux_64.tar.gz
GENEMARKES_KEY=sources/gm_key_64
GENEMARKES_PREFIX_ARG=$(INSTALL_BASE)/genemarkes/$(GENEMARKES_VERSION)


SIGNALP_VERSION=5.0
SIGNALP_PREFIX_ARG=$(INSTALL_BASE)/signalp/$(SIGNALP_VERSION)
SIGNALP_TAR=sources/signalp-5.0.Linux.tar.gz

# Repeat masking stuff
RED_TAG=v2.0
RED_REPO=https://github.com/TulsaBioinformaticsToolsmith/Red.git
RED_PREFIX_ARG=$(INSTALL_BASE)/red/$(RED_TAG)

RMBLAST_VERSION=2.9.0+
RMBLAST_URL=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-src.tar.gz
RMBLAST_PATCH_URL=http://www.repeatmasker.org/isb-2.9.0+-rmblast.patch.gz
RMBLAST_PREFIX_ARG=$(INSTALL_BASE)/rmblast/$(RMBLAST_VERSION)

TRF_VERSION=409
TRF_URL=http://tandem.bu.edu/trf/downloads/trf409.linux64
TRF_PREFIX_ARG=$(INSTALL_BASE)/trf/$(TRF_VERSION)

NSEG_URL=ftp://ftp.ncbi.nih.gov/pub/seg/nseg/
NSEG_PREFIX_ARG=$(INSTALL_BASE)/nseg

REPEATSCOUT_VERSION=1.0.5
REPEATSCOUT_URL=http://www.repeatmasker.org/RepeatScout-1.0.5.tar.gz
REPEATSCOUT_PREFIX_ARG=$(INSTALL_BASE)/repeatscout/$(REPEATSCOUT_VERSION)

RECON_VERSION=1.08
RECON_URL=http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz
RECON_PREFIX_ARG=$(INSTALL_BASE)/recon/$(RECON_VERSION)

RMASK_VERSION=4-0-9
RMASK_URL=http://www.repeatmasker.org/RepeatMasker-open-$(RMASK_VERSION).tar.gz
RMASK_PREFIX_ARG=$(INSTALL_BASE)/repeatmasker/$(RMASK_VERSION)

RMOD_VERSION=1.0.11
RMOD_URL=http://www.repeatmasker.org/RepeatModeler/RepeatModeler-open-$(RMOD_VERSION).tar.gz
RMOD_PREFIX_ARG=$(INSTALL_BASE)/repeatmodeller/$(RMOD_VERSION)

RM_LIB_ARG=/data/rmlib

RMASK_META_TAR=RepeatMaskerMetaData-20181026.tar.gz
RMASK_META_URL=http://www.repeatmasker.org/libraries/RepeatMaskerMetaData-20181026.tar.gz

REPBASE_TAR=RepBaseRepeatMaskerEdition-20181026.tar.gz
REPBASE_URL=http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20181026.tar.gz
REPBASE_USERNAME:=
REPBASE_PASSWORD:=


all: docker

docker: docker/general docker/alignment docker/rnaseq docker/annotation

# Singularity is distributed as a monolith
singularity: singularity/panann.sif

docker/general: docker/htslib docker/genometools docker/python3 docker/posix
docker/alignment: docker/minimap2 docker/spaln docker/gmap
docker/rnaseq: docker/trinity docker/star docker/stringtie
docker/annotation: docker/augustus docker/codingquarry docker/busco docker/pasa

# Proprietary software won't be built by default.
docker/proprietary: docker/genemarkes docker/signalp

clean:
	docker stop $(shell docker ps -aq) || true
	docker rm $(shell docker ps -a -q) || true
	docker rmi $(shell docker images -q) || true
	rm -rf -- docker

tidy:
	# Just removes untagged images to free up space
	#docker rmi -f $(shell docker images -f "dangling=true" -q) || true
	docker image prune -f

# Generic rules
singularity/%.sif: docker/%
	@mkdir -p "$(dir $@)"
	singularity build "$@" "docker-daemon://$(call make_image_name,$(notdir $(basename $@)))"

docker/%.tar.gz: docker/%
	@mkdir -p "$(dir $@)"
	docker save "$(call make_image_name,$(basename $(basename $@)))" | gzip > $@


docker/base: base.Dockerfile
	docker build \
	  --build-arg DEBIAN_VERSION=$(DEBIAN_VERSION) \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/python3: python3.Dockerfile
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/posix: posix.Dockerfile
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


# RNAseq stuff

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

docker/trinity: docker/jellyfish docker/bowtie2 docker/htslib docker/salmon
docker/trinity: trinity.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg JELLYFISH_IMAGE="$(call make_image_name,jellyfish)" \
	  --build-arg BOWTIE2_IMAGE="$(call make_image_name,bowtie2)" \
	  --build-arg SALMON_IMAGE="$(call make_image_name,salmon)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg TRINITY_TAG="$(TRINITY_TAG)" \
	  --build-arg TRINITY_REPO="$(TRINITY_REPO)" \
	  --build-arg TRINITY_PREFIX_ARG="$(TRINITY_PREFIX_ARG)" \
	  --build-arg JELLYFISH_VERSION="$(JELLYFISH_VERSION)" \
	  --build-arg JELLYFISH_PREFIX_ARG="$(JELLYFISH_PREFIX_ARG)" \
	  --build-arg BOWTIE2_TAG="$(BOWTIE2_TAG)" \
	  --build-arg BOWTIE2_PREFIX_ARG="$(BOWTIE2_PREFIX_ARG)" \
	  --build-arg SALMON_TAG="$(SALMON_TAG)" \
	  --build-arg SALMON_PREFIX_ARG="$(SALMON_PREFIX_ARG)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@

docker/star: star.Dockerfile docker/base docker/htslib
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg STAR_VERSION="$(STAR_VERSION)" \
	  --build-arg STAR_URL="$(STAR_URL)" \
	  --build-arg STAR_PREFIX_ARG="$(STAR_PREFIX_ARG)" \
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
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg STRINGTIE_VERSION="$(STRINGTIE_VERSION)" \
	  --build-arg STRINGTIE_URL="$(STRINGTIE_URL)" \
	  --build-arg STRINGTIE_PREFIX_ARG="$(STRINGTIE_PREFIX_ARG)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@

# General tools

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


# Aligners

docker/minimap2: minimap2.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg MINIMAP2_TAG="$(MINIMAP2_TAG)" \
	  --build-arg MINIMAP2_REPO="$(MINIMAP2_REPO)" \
	  --build-arg MINIMAP2_PREFIX_ARG="$(MINIMAP2_PREFIX_ARG)" \
	  --build-arg K8_VERSION="$(K8_VERSION)" \
	  --build-arg K8_URL="$(K8_URL)" \
	  --build-arg K8_PREFIX_ARG="$(K8_PREFIX_ARG)" \
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


# Annotation tools

docker/augustus: augustus.Dockerfile docker/base docker/htslib
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg AUGUSTUS_COMMIT="$(AUGUSTUS_COMMIT)" \
	  --build-arg AUGUSTUS_REPO="$(AUGUSTUS_REPO)" \
	  --build-arg AUGUSTUS_PREFIX_ARG="$(AUGUSTUS_PREFIX_ARG)" \
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

docker/pasa: pasa.Dockerfile docker/base docker/gmap docker/fasta docker/htslib
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg GMAP_IMAGE="$(call make_image_name,gmap)" \
	  --build-arg FASTA_IMAGE="$(call make_image_name,fasta)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg PASA_TAG="$(PASA_TAG)" \
	  --build-arg PASA_REPO="$(PASA_REPO)" \
	  --build-arg PASA_PREFIX_ARG="$(PASA_PREFIX_ARG)" \
	  --build-arg GMAP_VERSION="$(GMAP_VERSION)" \
	  --build-arg GMAP_PREFIX_ARG="$(GMAP_PREFIX_ARG)" \
	  --build-arg FASTA_VERSION="$(FASTA_VERSION)" \
	  --build-arg FASTA_PREFIX_ARG="$(FASTA_PREFIX_ARG)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
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

docker/gemoma: gemoma.Dockerfile docker/base docker/mmseqs docker/htslib
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg MMSEQS_IMAGE="$(call make_image_name,mmseqs)" \
	  --build-arg MMSEQS_TAG="$(MMSEQS_TAG)" \
	  --build-arg MMSEQS_PREFIX_ARG="$(MMSEQS_PREFIX_ARG)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --build-arg GEMOMA_VERSION="$(GEMOMA_VERSION)" \
	  --build-arg GEMOMA_URL="$(GEMOMA_URL)" \
	  --build-arg GEMOMA_PREFIX_ARG="$(GEMOMA_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/busco: busco.Dockerfile docker/base docker/augustus
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg AUGUSTUS_IMAGE="$(call make_image_name,augustus)" \
	  --build-arg BUSCO_COMMIT="$(BUSCO_COMMIT)" \
	  --build-arg BUSCO_REPO="$(BUSCO_REPO)" \
	  --build-arg BUSCO_PREFIX_ARG="$(BUSCO_PREFIX_ARG)" \
	  --build-arg AUGUSTUS_COMMIT="$(AUGUSTUS_COMMIT)" \
	  --build-arg AUGUSTUS_PREFIX_ARG="$(AUGUSTUS_PREFIX_ARG)" \
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

docker/braker_genemarkes: braker_genemarkes.Dockerfile docker/base docker/braker docker/augustus docker/htslib docker/bamtools docker/spaln docker/genemarkes
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg BRAKER_IMAGE="$(call make_image_name,braker)" \
	  --build-arg BRAKER_COMMIT="$(BRAKER_COMMIT)" \
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
	  --build-arg GENEMARKES_IMAGE="$(call make_image_name,genemarkes)" \
	  --build-arg GENEMARKES_VERSION="$(GENEMARKES_VERSION)" \
	  --build-arg GENEMARKES_PREFIX_ARG="$(GENEMARKES_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@

docker/signalp: signalp.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg SIGNALP_VERSION="$(SIGNALP_VERSION)" \
	  --build-arg SIGNALP_TAR="$(SIGNALP_TAR)" \
	  --build-arg SIGNALP_PREFIX_ARG="$(SIGNALP_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


## The monolith!!!!
docker/panann: panann.Dockerfile docker/base docker/jellyfish docker/bowtie2 docker/salmon docker/trinity docker/star docker/stringtie docker/htslib docker/genometools docker/minimap2 docker/spaln docker/gmap docker/augustus docker/codingquarry docker/busco
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg JELLYFISH_IMAGE="$(call make_image_name,jellyfish)" \
	  --build-arg JELLYFISH_VERSION="$(JELLYFISH_VERSION)" \
	  --build-arg JELLYFISH_PREFIX_ARG="$(JELLYFISH_PREFIX_ARG)" \
	  --build-arg BOWTIE2_IMAGE="$(call make_image_name,bowtie2)" \
	  --build-arg BOWTIE2_TAG="$(BOWTIE2_TAG)" \
	  --build-arg BOWTIE2_PREFIX_ARG="$(BOWTIE2_PREFIX_ARG)" \
	  --build-arg SALMON_IMAGE="$(call make_image_name,salmon)" \
	  --build-arg SALMON_TAG="$(SALMON_TAG)" \
	  --build-arg SALMON_PREFIX_ARG="$(SALMON_PREFIX_ARG)" \
	  --build-arg TRINITY_IMAGE="$(call make_image_name,trinity)" \
	  --build-arg TRINITY_TAG="$(TRINITY_TAG)" \
	  --build-arg TRINITY_PREFIX_ARG="$(TRINITY_PREFIX_ARG)" \
	  --build-arg STAR_IMAGE="$(call make_image_name,star)" \
	  --build-arg STAR_VERSION="$(STAR_VERSION)" \
	  --build-arg STAR_PREFIX_ARG="$(STAR_PREFIX_ARG)" \
	  --build-arg STRINGTIE_IMAGE="$(call make_image_name,stringtie)" \
	  --build-arg STRINGTIE_VERSION="$(STRINGTIE_VERSION)" \
	  --build-arg STRINGTIE_PREFIX_ARG="$(STRINGTIE_PREFIX_ARG)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --build-arg BCFTOOLS_TAG="$(BCFTOOLS_TAG)" \
	  --build-arg BCFTOOLS_PREFIX_ARG="$(BCFTOOLS_PREFIX_ARG)" \
	  --build-arg GENOMETOOLS_IMAGE="$(call make_image_name,genometools)" \
	  --build-arg GENOMETOOLS_VERSION="$(GENOMETOOLS_VERSION)" \
	  --build-arg GENOMETOOLS_PREFIX_ARG="$(GENOMETOOLS_PREFIX_ARG)" \
	  --build-arg MINIMAP2_IMAGE="$(call make_image_name,minimap2)" \
	  --build-arg MINIMAP2_TAG="$(MINIMAP2_TAG)" \
	  --build-arg MINIMAP2_PREFIX_ARG="$(MINIMAP2_PREFIX_ARG)" \
	  --build-arg K8_VERSION="$(K8_VERSION)" \
	  --build-arg K8_PREFIX_ARG="$(K8_PREFIX_ARG)" \
	  --build-arg SPALN_IMAGE="$(call make_image_name,spaln)" \
	  --build-arg SPALN_TAG="$(SPALN_TAG)" \
	  --build-arg SPALN_PREFIX_ARG="$(SPALN_PREFIX_ARG)" \
	  --build-arg GMAP_IMAGE="$(call make_image_name,gmap)" \
	  --build-arg GMAP_VERSION="$(GMAP_VERSION)" \
	  --build-arg GMAP_PREFIX_ARG="$(GMAP_PREFIX_ARG)" \
	  --build-arg AUGUSTUS_IMAGE="$(call make_image_name,augustus)" \
	  --build-arg AUGUSTUS_COMMIT="$(AUGUSTUS_COMMIT)" \
	  --build-arg AUGUSTUS_PREFIX_ARG="$(AUGUSTUS_PREFIX_ARG)" \
	  --build-arg CODINGQUARRY_IMAGE="$(call make_image_name,codingquarry)" \
	  --build-arg CODINGQUARRY_VERSION="$(CODINGQUARRY_VERSION)" \
	  --build-arg CODINGQUARRY_PREFIX_ARG="$(CODINGQUARRY_PREFIX_ARG)" \
	  --build-arg BUSCO_IMAGE="$(call make_image_name,busco)" \
	  --build-arg BUSCO_COMMIT="$(BUSCO_COMMIT)" \
	  --build-arg BUSCO_PREFIX_ARG="$(BUSCO_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


# Repeat masking stuff

docker/edge_base: base.Dockerfile
	docker build \
	  --build-arg DEBIAN_VERSION="buster-20190228-slim" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/red: red.Dockerfile docker/edge_base docker/base
	docker build \
	  --build-arg BUILD_IMAGE="$(call make_image_name,edge_base)" \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg RED_TAG="$(RED_TAG)" \
	  --build-arg RED_REPO="$(RED_REPO)" \
	  --build-arg RED_PREFIX_ARG="$(RED_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/rmblast: rmblast.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg RMBLAST_VERSION="$(RMBLAST_VERSION)" \
	  --build-arg RMBLAST_URL="$(RMBLAST_URL)" \
	  --build-arg RMBLAST_PATCH_URL="$(RMBLAST_PATCH_URL)" \
	  --build-arg RMBLAST_PREFIX_ARG="$(RMBLAST_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/repeatmasker: repeatmasker.Dockerfile docker/base docker/rmblast
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg RMBLAST_IMAGE="$(call make_image_name,rmblast)" \
	  --build-arg RMBLAST_VERSION="$(RMBLAST_VERSION)" \
	  --build-arg RMBLAST_PREFIX_ARG="$(RMBLAST_PREFIX_ARG)" \
	  --build-arg TRF_VERSION="$(TRF_VERSION)" \
	  --build-arg TRF_URL="$(TRF_URL)" \
	  --build-arg TRF_PREFIX_ARG="$(TRF_PREFIX_ARG)" \
	  --build-arg NSEG_URL="$(NSEG_URL)" \
	  --build-arg NSEG_PREFIX_ARG="$(NSEG_PREFIX_ARG)" \
	  --build-arg REPEATSCOUT_VERSION="$(REPEATSCOUT_VERSION)" \
	  --build-arg REPEATSCOUT_URL="$(REPEATSCOUT_URL)" \
	  --build-arg REPEATSCOUT_PREFIX_ARG="$(REPEATSCOUT_PREFIX_ARG)" \
	  --build-arg RECON_VERSION="$(RECON_VERSION)" \
	  --build-arg RECON_URL="$(RECON_URL)" \
	  --build-arg RECON_PREFIX_ARG="$(RECON_PREFIX_ARG)" \
	  --build-arg RMASK_VERSION="$(RMASK_VERSION)" \
	  --build-arg RMASK_URL="$(RMASK_URL)" \
	  --build-arg RMOD_VERSION="$(RMOD_VERSION)" \
	  --build-arg RMOD_URL="$(RMOD_URL)" \
	  --build-arg RMOD_PREFIX_ARG="$(RMOD_PREFIX_ARG)" \
	  --build-arg RM_LIB_ARG="$(RM_LIB_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


downloads/Dfam.hmm.gz:
	mkdir -p downloads
	cd downloads \
	&& wget http://dfam.org/releases/Dfam_3.0/families/Dfam.hmm.gz


downloads/Dfam.embl.gz:
	mkdir -p downloads
	cd downloads \
	&& wget http://dfam.org/releases/Dfam_3.0/families/Dfam.embl.gz


downloads/$(RMASK_META_TAR):
	mkdir -p downloads
	cd downloads \
	&& wget -O $@ $(RMASK_META_URL)


downloads/$(REPBASE_TAR):
	mkdir -p downloads
	cd downloads \
	&& wget \
	     --user $(REPBASE_USERNAME) \
	     --password $(REPBASE_PASSWORD) \
	     -O $@ \
	     $(REPBASE_URL)


rmlib: downloads/Dfam.hmm.gz downloads/Dfam.embl.gz
rmlib: downloads/$(RMASK_META_TAR) downloads/$(REPBASE_TAR)
rmlib: docker/repeatmasker
	docker run \
	  --rm \
	  -v "$(PWD)/downloads:/downloads:rw" \
	  -v "$(PWD)/rmlib:$(RM_LIB_ARG):rw" \
	  -w /downloads \
	  $(call make_image_name,repeatmasker) \
	  prep_repeatmasker_lib.sh \
	    $(REPBASE_TAR) \
	    $(RMASK_META_TAR) \
	    Dfam.hmm.gz \
	    Dfam.embl.gz \
	    $(RM_LIB_ARG)
