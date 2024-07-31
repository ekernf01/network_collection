#!/bin/bash
docker build -t genga_endoderm_reanalysis_2024 .
docker tag genga_endoderm_reanalysis_2024 ekernf01/genga_endoderm_reanalysis_2024
docker login
docker push ekernf01/genga_endoderm_reanalysis_2024
output=/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/network_collection/setup/genga_definitive_endoderm_ATAC/from_to_docker/
mkdir $output
docker run --rm -it --mount type=bind,source=${output},destination=/from_to_docker ekernf01/genga_endoderm_reanalysis_2024 
