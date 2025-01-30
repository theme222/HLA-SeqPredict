FROM python:3.11
WORKDIR /app
COPY . .

# Download and install Miniconda
WORKDIR /app
# RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
#     bash miniconda.sh -b -p /opt/conda && \
#     rm miniconda.sh
# ENV PATH=/opt/conda/bin:$PATH
# RUN conda --version
# RUN conda config --add channels defaults &&\
#     conda config --add channels bioconda &&\
#     conda config --add channels conda-forge &&\
#     conda config --set channel_priority strict

# # Installing hla-la
# RUN conda install -y hla-la
# RUN mkdir /opt/conda/opt/hla-la/graphs &&\
#     cd /opt/conda/opt/hla-la/graphs &&\
#     wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz && \
#     tar -xvzf PRG_MHC_GRCh38_withIMGT.tar.gz && \
#     cd /opt/conda/opt/hla-la/src && \
#     wget https://www.dropbox.com/s/mnkig0fhaym43m0/reference_HLA_ASM.tar.gz && \
#     tar -xvzf reference_HLA_ASM.tar.gz && \
#     ../bin/HLA-LA --action prepareGraph --PRG_graph_dir ../graphs/PRG_MHC_GRCh38_withIMGT

RUN apt-get update && apt-get install -y ninja-build

# Installing python libraries
WORKDIR /app
RUN pip install --no-cache-dir -r requirements.txt

# # Unzipping reference
# WORKDIR /app/references
# RUN tar -xf chr6.tar.xz

# # Installing bwa
# WORKDIR /app
# RUN git clone https://github.com/lh3/bwa.git
# WORKDIR /app/bwa
# RUN make
# ENV PATH="/app/bwa:${PATH}"
# RUN chmod +x ./*

# # Installing samtools
# WORKDIR /app
# RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.20/samtools-1.20.tar.bz2/download
# RUN tar -xvf download
# WORKDIR /app/samtools-1.20
# RUN make
# ENV PATH="$PATH:/app/samtools-1.20"
# WORKDIR /app
# RUN rm -f download

# # Indexing reference (This takes like 3 minutes)
# WORKDIR /app/references
# RUN bwa index chr6.fa

RUN apt-get update && apt-get install -y \
    build-essential \
    uuid-dev \
    libgpgme-dev \
    squashfs-tools \
    libseccomp-dev \
    pkg-config \
    cryptsetup-bin

RUN wget https://github.com/sylabs/singularity/releases/download/v4.1.4/singularity-ce_4.1.4-jammy_amd64.deb

RUN apt-get install ./singularity-ce_4.1.4-jammy_amd64.deb -y

# Run
WORKDIR /app
CMD ["python", "./main.py"]