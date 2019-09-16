ARG IMAGE

FROM "${IMAGE}" as builder

ARG MASH_TAG
ARG MASH_REPO="https://github.com/lh3/minimap2.git"
ARG MASH_PREFIX_ARG="/opt/minimap2/${MINIMAP2_TAG}"
ENV MASH_PREFIX="${MINIMAP2_PREFIX_ARG}"
