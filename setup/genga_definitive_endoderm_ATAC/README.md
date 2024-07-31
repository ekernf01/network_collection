## Genga et al. endoderm differentiation ATAC analysis

This folder documents how we obtained a TF-target network via motif analysis of the Genga et al. endoderm directed differentiation ATAC-seq data. You can repeat our analysis by running a docker container using the following command. You may need to modify the output folder.

```sh
> mkdir from_to_docker
> docker run --rm -it --mount type=bind,source=$PWD/from_to_docker/,destination=/from_to_docker/ ekernf01/genga_endoderm_reanalysis_2024 
(in the container)> ./reprocess_genga_atac.sh
```

#### More detail 

This docker container is from an image built via the `Dockerfile` in this folder. The script `reprocess_genga_atac.sh` will download the fastq files, download hg38, align reads via bowtie2, call peaks via macs2, and run motif analysis via gimmemotifs via celloracle. This roughly follows the [original analysis](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3630246), but with hg38. The script will write all files (FQ, BAM, BED, and TSV network) to the output location you specified above. You can also view them or tinker with them inside the container. The container will be deleted upon exit, but your output files will persist.