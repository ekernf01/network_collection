dir.create("not_ready/ANANSE/networks", recursive = T, showWarnings = F)
for( celltype in c("astrocyte",
                 "fibroblast",
                 "heart",
                 "hIPS",
                 "keratinocyte",
                 "liver",
                 "macrophage",
                 "osteoblast")){
  target = paste0("not_ready/ANANSE/networks/osteoblast", ".txt.gz")
  if(!file.exists(target)){
    curl::curl_download(paste0("https://zenodo.org/record/4809063/files/", celltype,  ".network.txt.gz?download=1"),
                        target)
  }
}

for( celltype in c(
  "adrenal_gland",
  "bone_marrow",
  "brain",
  "cervix",
  "colon",
  "esophagus",
  "heart",
  "intestine",
  "liver",
  "lung",
  "ovary",
  "pancreas",
  "prostate",
  "skeletal_muscle",
  "skin",
  "spleen",
  "stomach",
  "thymus")){
  target = paste0("not_ready/ANANSE/networks/osteoblast", ".txt.gz")
  if(!file.exists(target)){
    curl::curl_download(paste0("https://zenodo.org/record/4807802/files/", celltype,  ".network.txt.gz?download=1"), 
                        target)
  }
}

# TODO: filter by probability?? reformat as expected.


