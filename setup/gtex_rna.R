#' Tidy up the formatting of the gtex coexpression networks. 
#' The files are originally named 1-35 and we want to label them by tissue.
#' Also, they are symmetric in nature and only 1/2 the edges are explicitly stored, so we symmetrize.
#' This should leave behind both .parquet and .txt files. You can delete all the txt's. 

fix_gtex_rna_filenames = function(){
  where = here::here("networks", "gtex_rna", "networks")
  names_current = list.files( where, full.names = T)
  numbers_current = names_current %>% basename %>% gsub("^.*_", "", .) %>% as.numeric
  tissues = read.table( sep = "\t", file = file.path(dirname(where), "gtex_rna_tissue_names.txt"))[[1]]
  numbers_tissues = tissues %>% gsub("^ ", "", .) %>% gsub(" .*$", "", .) %>% as.numeric
  names_desired = tissues %>% 
    gsub("^[0-9]* ", "", .) %>% gsub(" - | |\\(", "_", .) %>% gsub("\\)", "", .) %>%
    paste0(".txt") %>% file.path(where, .)
  
  names_current = names_current[order(numbers_current)]
  names_desired = names_desired[order(numbers_tissues)]
  
  file.rename(from = names_current, to = names_desired)
}

make_gtex_networks_symmetric = function(){
  for(tissue in list_subnetworks("gtex_rna")){
    X = load_grn_by_subnetwork( "gtex_rna", tissue)
    X = rbind(X, setNames(X[c(2,1,3)], colnames(X)))
    X = X[!duplicated(X),]
    write_grn_by_subnetwork( "gtex_rna", tissue, grn_df = X )
  }
}

fix_gtex_rna_filenames()
make_gtex_networks_symmetric()

