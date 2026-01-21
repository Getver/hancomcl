# FROM ubuntu:22.04
# docker run -it --rm sm:chip_0.2
FROM rocker/r-ver:4.4.0 AS base-amd64

USER root
WORKDIR /opt/

RUN apt-get -y -m update && apt-get install -y build-essential curl libidn12 libnet-perl perl-doc liblmdb-dev wget libsqlite3-dev cmake \
    libgomp1 libxml-simple-perl libjson-perl parallel vmtouch cpanminus python3-pip texlive-latex-base vim python2 git libbz2-dev \
    liblzma-dev libcurl4-openssl-dev zlib1g-dev libncurses5-dev autoconf automake && \
    rm -rf /var/lib/apt/lists/*  && \
    cpanm HTML::Entities

RUN mkdir -p /tmp
RUN cd /tmp && git clone https://github.com/samtools/htslib.git && cd htslib && git submodule update --init --recursive && make -j 4 && cd .. && cp ./htslib/tabix /usr/local/bin/ && cp ./htslib/bgzip /usr/local/bin/ && cp ./htslib/htsfile /usr/local/bin/ && strip /usr/local/bin/tabix; true && strip /usr/local/bin/bgzip; true && strip /usr/local/bin/htsfile; true && git clone https://github.com/samtools/samtools.git && cd samtools && make -j 4 && cp ./samtools /usr/local/bin/ && cd .. && strip /usr/local/bin/samtools; true 
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/tmp/htslib

RUN pip install biopython pandas

RUN mkdir -p /library_QC

# COPY library /library_QC/library
COPY tool /library_QC/tool
COPY advnorm-1.1 /library_QC/advnorm-1.1
WORKDIR /library_QC/tool/
# RUN wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_amd_avx2_20241111.zip && unzip plink2_linux_amd_avx2_20241111.zip

RUN R -e "install.packages('ggplot2')"
RUN R -e "install.packages('reshape2')"
RUN R -e "install.packages('R6')"
RUN R -e "install.packages('/library_QC/tool/SNPolisher_3.0.tar.gz',repos=NULL, type='source')"

RUN mkdir -p /code_SM
COPY run.sh KORV1.sh KORV2.sh PMDA.sh header.txt PangenomiX.sh HCLV.sh pharmacoFocus.sh make_table.py KBA_A.sh KBA_B.sh plot.R /code_SM/

WORKDIR /code_SM/


# sm:chip_0.1
# CMD [ "/bin/bash" ]

# sm:chip_0.2
ENTRYPOINT [ "bash", "/code_SM/run.sh" ]

