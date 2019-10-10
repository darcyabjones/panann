# pclust

This folder contains the Dockerfiles to build containers to run the pipeline.
The software can be build individually or as a monolithic container, depending on your preference.

Prebuilt released/tagged versions are available at:

- https://hub.docker.com/r/darcyabjones/pann
- https://cloud.sylabs.io/library/_container/5d71f338cb093ee1142b8791


To pull those images:


```
singularity pull library://darcyabjones/default/panann:v0.0.1

# or

sudo docker pull darcyabjones/panann:v0.0.1
```

Replacing `v0.0.1` with whatever version you actually want.


The final image contains:

- [MMSeqs](https://github.com/soedinglab/MMseqs2)

Versions of individual software are available as environment variables and labels in the containers.


## Versioning

The container tags will correspond to a release of the panann pipeline.
They are intended to be paired so that the correct commands are used etc.
It is possible that two docker builds with different tags will be identical because the pipeline has changed but the software dependencies haven't.


## Optimisation

Because MMSeqs2 requires a processor capable of using at least SSE4 instructions, this pipeline also has that minimum requirement.
If you're unsure about your processors capability, on linux you can run the command `grep "sse\|avx" /proc/cpuinfo` to get a list of your CPUs capabilities.

We've done our best to dynamically handle different vectorisation instruction sets in the image.
There should be at least a version that works with SSE and one that works with AVX2, and the correct version should be selected depending on your own processor.

If you really want to squeeze every bit of performance out of whatever you have, you might need to compile the software yourself using the `-march=native` gcc compiler options (or equivalent for your compiler).


## Build system

Because there are a few different bits of software, the container is built locally and pushed rather than using auto-build tools.
Dockerfiles and the Makefile that coordinates building is available at <github.com/darcyabjones/panann/tree/master/containers>.

If you're uncomfortable about running unverified containers you might like to build it yourself.


To build the images you will need [Docker](https://docs.docker.com/install/), [Make](https://www.gnu.org/software/make/), and optionally [Singularity](https://sylabs.io/guides/latest/user-guide/) installed.
You will also need root permission for your computer.

The following commands should work for Linux, Mac, and Linux emulators/VMs (e.g. WSL or Cygwin).

```
# To build all docker images separately.
sudo make docker/all

# To build the monolithic docker image
sudo make docker/pannn

# To build a specific dockerfile
sudo make docker/augustus

# To build all singularity images.
sudo make singularity/all

# To build a specific image.
sudo make singularity/mmseqs.sif

# To build the monolithic image.
sudo make singularity/panann.sif

# To tidy the docker image registry and singularity cache.
# Note that this may also remove docker images created by other methods.
sudo make tidy

# To remove all docker images and singularity images from your computer.
# Will remove EVERYTHING, aka. the "help i've run out of root disk space" command.
sudo make clean
```

Singularity containers are placed in a `singularity` subdirectory with the extension `.sif`.


Note that singularity also leaves some potentially large tar files in tmp folders that wont be removed automatically.
On my computers (ubuntu/fedora) these are put in `/var/tmp/docker-*`.


## Proprietary software

Genemark-ES and SignalP both have restricted licences for non-academic users.
This means that they can't be distributed by Docker or singularity hub.

We use DeepSig as a default alternative to SignalP.
You can probably live without genemark annotations if you have a good Augustus model and some good RNAseq data.

To use GeneMark or SignalP you can create containers for them in this directory.

Download the source files for GeneMark-ES and SignalP 5, and the GeneMark key file into a `./sources` subfolder.

```bash
ls sources
gm_key_64 signalp-5.0.Linux.tar.gz gm_et_linux_64.tar.gz

sudo make docker/proprietary
sudo make singularity/proprietary
```

This will make the docker containers and the singularity containers, which you can now use in the pipeline.
