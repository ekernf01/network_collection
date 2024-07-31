FROM kenjikamimoto126/celloracle_ubuntu:0.18.0
# needed for apt-get to not hang
ENV DEBIAN_FRONTEND=noninteractive 
RUN apt-get update && apt-get install -y bowtie2  
RUN apt-get update && apt-get install -y samtools  
RUN apt-get update && apt-get install -y bedtools  
RUN apt-get update && apt-get install -y r-base  
RUN apt-get update && apt-get install -y bowtie2 
RUN apt-get update && apt-get install -y pigz 
RUN pip install MACS2
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz
RUN tar xzvf sratoolkit.2.11.3-ubuntu64.tar.gz
ENV PATH="${PATH}:/sratoolkit.2.11.3-ubuntu64/bin"
RUN pip install git+https://github.com/ekernf01/pereggrn_networks.git 
COPY reprocess_genga_atac.sh /
COPY generate_network.py /
COPY TruSeq3-PE-2.fa /
