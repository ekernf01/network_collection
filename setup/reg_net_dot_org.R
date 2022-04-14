#' Tidy up the formatting of the ENCODE DHS networks. 
#'
fix_regulatorynetworks.org = function(grn_name){
  where = here("networks", grn_name, "networks")
  names_current = list.files( where, full.names = T) %>% file.path("genes.regulate.genes") 
  names_desired = names_current %>% dirname %>% paste0(".txt")
  file.rename(from = names_current, to = names_desired)
  lapply(names_desired, function(my_chubby_file) system(paste("gzip", my_chubby_file)))
}
fix_regulatorynetworks.org("regulatorynetworks.org_human")
fix_regulatorynetworks.org("regulatorynetworks.org_mouse")
fix_one = function(grn_name, subnetwork_name) {
  grn_df = load_grn_by_subnetwork(grn_name, subnetwork_name) %>% 
    extract(4:6) %>%
    # This changes the order from target, regulator to regulator, target.
    set_colnames(c("target", "regulator", "weight")) %>%
    extract(EXPECTED_GRN_COLNAMES)
  write_grn_by_subnetwork( grn_name, subnetwork_name %>% gsub(".txt.gz", ".csv", .), grn_df )
}
iterate_within_grn(grn_name="regulatorynetworks.org_mouse", fix_one)
