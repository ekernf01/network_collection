dir.create("not_ready/STRING/networks", recursive = T, showWarnings = F)
downloaded_file_name = "not_ready/STRING/networks/v11.5"
if(!file.exists(downloaded_file_name)){
  curl::curl_download(paste0("https://stringdb-static.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz"),
                      downloaded_file_name)
}
tissue_name = basename(downloaded_file_name)
network_name = basename(dirname(dirname(downloaded_file_name)))
cat(tissue_name, "\n")
X = read.table(downloaded_file_name, header = T)
# Convert ensembl peptide ID's to gene symbols
X$protein1 %<>% gsub("^9606.", "", .)
X$protein2 %<>% gsub("^9606.", "", .)
library('biomaRt')
G_list <- getBM(filters= "ensembl_peptide_id", 
                attributes= c("ensembl_peptide_id","hgnc_symbol"),
                values=union(X$protein1, X$protein2),
                mart= useDataset("hsapiens_gene_ensembl", useMart("ensembl")))
X %<>% merge(G_list,by.x="protein1",by.y="ensembl_peptide_id")
X %<>% dplyr::rename(regulator = hgnc_symbol)
X %<>% merge(G_list,by.x="protein2",by.y="ensembl_peptide_id")
X %<>% dplyr::rename(target = hgnc_symbol)
X %<>% dplyr::rename(weight = combined_score)
X = X[EXPECTED_GRN_COLNAMES]
write_grn_by_subnetwork(network_name, tissue_name, grn_df = X)






