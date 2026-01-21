# FROM ubuntu:24.04
ARG BASE_IMAGE=rocker/r-ver
ARG amd64_tag=4.4.2
FROM ${BASE_IMAGE}:${amd64_tag} AS base-amd64

WORKDIR /usr/local/src

RUN apt-get update && apt-get install -y git curl tar xz-utils libglpk40 libglpk-dev libopenblas-dev \
    liblapack-dev libcurl4-openssl-dev libxml2-dev r-base libssl-dev libatlas-base-dev gfortran libhdf5-dev wget xxd

RUN curl -o cellranger-9.0.0.tar.xz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.0.tar.gz?Expires=1737640025&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=aLgEYBc-Rf-KP535ipMPZ4~5Ejj9UAE4Eopdg0GHaXe6wqBVEixGQvksJnjca1fvTheWj8D0OyyLTZI9Mp2tn24V0rpKg3ZcsUUQhqKzwnvXtsJgoWser3S8m17pKR9-THf-V5qRs~G9naqTQ3XzcWbdkMDdncnDE9k6nxOgWV7mdDXhTnsSV3KciHvaQIdt3NR2cvEF5oTYOqTGs-I1STArR~SDyZ1nB6Wss7aXlPN5LJ72N8bHaQO4~hxh5Bf6sSvi-X0Z1laUz6A04tX0cmp34633iIUnWsp2VNNxX~cVwVPCyboQXSJmKffu0ORuRgbNNvYWJK0RBjUL2Xk~mg__" && tar -xvf cellranger-9.0.0.tar.xz && ln -s /usr/local/src/cellranger-9.0.0/bin/cellranger /usr/local/bin/cellranger

RUN mkdir reference && cd reference && curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz" && tar -xvf refdata-gex-GRCh38-2024-A.tar.gz

RUN R -e 'install.packages(c("Seurat", "hdf5r"))'
RUN R -e 'install.packages("remotes")' &&  R -e 'remotes::install_github("mojaveazure/seurat-disk")'

RUN cd /usr/local/src && wget https://github.com/alexdobin/STAR/archive/2.7.1a.tar.gz
RUN tar -xzf 2.7.1a.tar.gz && \
    cd STAR-2.7.1a/source && make STAR && make STAR CXXFLAGS_SIMD=sse


RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")' && R -e 'BiocManager::install(c("GSEABase", "GSVA", "singscore"), ask = FALSE)'


RUN R -e 'install.packages("devtools")' && \
    R -e 'devtools::install_version("crossmatch", version = "1.3.1", repos = "http://cran.us.r-project.org", ask = FALSE)' && \
    R -e 'devtools::install_version("multicross", version = "2.1.0", repos = "http://cran.us.r-project.org", ask = FALSE)' && \
    R -e 'devtools::install_github("jackbibby1/SCPA", ask = FALSE)'



RUN mkdir /usr/local/src/code
COPY ./00_run.sh ./01_cellranger.sh ./02_seurat.R /usr/local/src/code/
COPY ./reference/* /usr/local/src/reference/

WORKDIR /home/data

CMD [ "/bin/bash" ]

# ENTRYPOINT [ "bash", "/code_SM/run.sh" ]
