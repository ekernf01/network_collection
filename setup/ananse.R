dir.create("not_ready/ANANSE/networks", recursive = T, showWarnings = F)
dir.create("not_ready/ANANSE_tissue/networks", recursive = T, showWarnings = F)
dir.create("not_ready/ANANSE_0.5/networks", recursive = T, showWarnings = F)
dir.create("not_ready/ANANSE_tissue_0.5/networks", recursive = T, showWarnings = F)

reformat_ananse = function(downloaded_file_name){
  tissue_name = basename(downloaded_file_name)
  network_name = basename(dirname(dirname(downloaded_file_name)))
  strippit = function(x) tools::file_path_sans_ext(tools::file_path_sans_ext(x))
  if(strippit(tissue_name) %in% strippit(list_subnetworks(grn_name = paste0(network_name, "_0.5")))){
    cat("Skipping ", tissue_name, "\n")
    return()
  } else {
    cat("Doing ", tissue_name, "\n")
  }
  X = read.table(downloaded_file_name, header = T)
  X %<>% tidyr::separate("tf_target", into = c("tf", "target"), sep = "_")
  colnames(X) = EXPECTED_GRN_COLNAMES
  X$target %<>% factor
  X$regulator %<>% factor
  write_grn_by_subnetwork(grn_name = network_name, subnetwork_name = tissue_name, grn_df = X)
  write_grn_by_subnetwork(grn_name = paste0(network_name, "_0.8"), subnetwork_name = tissue_name, grn_df = subset(X, weight>=0.8))
}
for( celltype in c("astrocyte",
                 "fibroblast",
                 "heart",
                 "hIPS",
                 "keratinocyte",
                 "liver",
                 "macrophage",
                 "osteoblast")){
  target = paste0("not_ready/ANANSE/networks/", celltype, ".txt.gz")
  if(!file.exists(target)){
    curl::curl_download(paste0("https://zenodo.org/record/4809063/files/", celltype,  ".network.txt.gz?download=1"),
                        target)
  }
  reformat_ananse(target)
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
  target = paste0("not_ready/ANANSE_tissue/networks/", celltype, ".txt.gz")
  if(!file.exists(target)){
    curl::curl_download(paste0("https://zenodo.org/record/4807802/files/", celltype,  ".network.txt.gz?download=1"), 
                        target, quiet = F)
  }
  reformat_ananse(target)
}


