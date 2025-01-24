FROM bioconductor/tidyverse:RELEASE_3_20

LABEL org.opencontainers.image.authors='Mykle L. Hoban <mykle.hoban@hawaii.edu>'

# install tidyverse
# RUN /rocker_scripts/install_tidyverse.sh

# install packages
RUN Rscript -e 'BiocManager::install("VariantAnnotation")'
RUN Rscript -e 'install.packages(c("janitor","optparse"))'
