#' Tidy up the formatting of the gtex coexpression networks. 
#'
fix_gtex_rna = function(grn_name){
  where = here("networks", grn_name, "networks")
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
  lapply(names_desired, function(my_chubby_file) system(paste("gzip", my_chubby_file)))
}
fix_gtex_rna(grn_name = "gtex_rna")
