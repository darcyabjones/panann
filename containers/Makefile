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

MINIMAP_TAG=v2.16
K8_VERSION=0.2.4
MINIMAP_REPO=https://github.com/lh3/minimap2.git
K8_URL=https://github.com/attractivechaos/k8/releases/download/v$(K8_VERSION)/k8-$(K8_VERSION).tar.bz2
MINIMAP_PREFIX_ARG=$(INSTALL_BASE)/minimap2/$(MINIMAP_TAG)
K8_PREFIX_ARG=$(INSTALL_BASE)/k8/$(K8_VERSION)

SPALN_TAG=Ver.2.3.3
SPALN_REPO=https://github.com/ogotoh/spaln.git
SPALN_PREFIX_ARG=$(INSTALL_BASE)/spaln/$(SPALN_TAG)

GMAP_VERSION=2019-05-12
GMAP_URL=http://research-pub.gene.com/gmap/src/gmap-gsnap-$(GMAP_VERSION).tar.gz
GMAP_PREFIX_ARG=$(INSTALL_BASE)/gmap/$(GMAP_VERSION)

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


all: docker

docker: docker/general docker/alignment docker/rnaseq docker/annotation

# Singularity is distributed as a monolith
singularity: singularity/panann.sif

docker/general: docker/htslib docker/genometools
docker/alignment: docker/minimap2 docker/spaln docker/gmap
docker/rnaseq: docker/trinity docker/star docker/stringtie
docker/annotation: docker/augustus docker/codingquarry docker/busco

# Proprietary software won't be build by default.
docker/proprietary:

clean:
	docker stop $(shell docker ps -aq) || true
	docker rm $(shell docker ps -a -q) || true
	docker rmi $(shell docker images -q) || true
	rm -rf -- singularity
	rm -rf -- docker

tidy:
	# Just removes untagged images to free up space
	docker rmi $(shell docker images | grep "^<none>" | awk '{print $3}')

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

docker/bowtie2: bowtie2.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg BOWTIE2_TAG="$(BOWTIE2_TAG)" \
	  --build-arg BOWTIE2_REPO="$(BOWTIE2_REPO)" \
	  --build-arg BOWTIE2_PREFIX_ARG="$(BOWTIE2_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .

docker/salmon: salmon.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg SALMON_TAG="$(SALMON_TAG)" \
	  --build-arg SALMON_REPO="$(SALMON_REPO)" \
	  --build-arg SALMON_PREFIX_ARG="$(SALMON_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .

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

docker/star: star.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg STAR_VERSION="$(STAR_VERSION)" \
	  --build-arg STAR_URL="$(STAR_URL)" \
	  --build-arg STAR_PREFIX_ARG="$(STAR_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .

docker/stringtie: stringtie.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg STRINGTIE_VERSION="$(STRINGTIE_VERSION)" \
	  --build-arg STRINGTIE_URL="$(STRINGTIE_URL)" \
	  --build-arg STRINGTIE_PREFIX_ARG="$(STRINGTIE_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .

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

docker/genometools: genometools.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg GENOMETOOLS_VERSION="$(GENOMETOOLS_VERSION)" \
	  --build-arg GENOMETOOLS_URL="$(GENOMETOOLS_URL)" \
	  --build-arg GENOMETOOLS_PREFIX_ARG="$(GENOMETOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .

# Aligners

docker/minimap2: minimap2.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg MINIMAP_TAG="$(MINIMAP_TAG)" \
	  --build-arg MINIMAP_REPO="$(MINIMAP_REPO)" \
	  --build-arg MINIMAP_PREFIX_ARG="$(MINIMAP_PREFIX_ARG)" \
	  --build-arg K8_VERSION="$(K8_VERSION)" \
	  --build-arg K8_URL="$(K8_URL)" \
	  --build-arg K8_PREFIX_ARG="$(K8_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .

docker/spaln: spaln.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg SPALN_TAG="$(SPALN_TAG)" \
	  --build-arg SPALN_REPO="$(SPALN_REPO)" \
	  --build-arg SPALN_PREFIX_ARG="$(SPALN_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .

docker/gmap: gmap.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg GMAP_VERSION="$(GMAP_VERSION)" \
	  --build-arg GMAP_URL="$(GMAP_URL)" \
	  --build-arg GMAP_PREFIX_ARG="$(GMAP_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .

# Annotation tools

docker/augustus: augustus.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg AUGUSTUS_COMMIT="$(AUGUSTUS_COMMIT)" \
	  --build-arg AUGUSTUS_REPO="$(AUGUSTUS_REPO)" \
	  --build-arg AUGUSTUS_PREFIX_ARG="$(AUGUSTUS_PREFIX_ARG)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_REPO="$(HTSLIB_REPO)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_REPO="$(SAMTOOLS_REPO)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .

docker/codingquarry: codingquarry.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg CODINGQUARRY_VERSION="$(CODINGQUARRY_VERSION)" \
	  --build-arg CODINGQUARRY_URL="$(CODINGQUARRY_URL)" \
	  --build-arg CODINGQUARRY_PREFIX_ARG="$(CODINGQUARRY_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .

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
