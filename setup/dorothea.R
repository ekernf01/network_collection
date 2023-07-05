BiocManager::install("dorothea", version = 3.16)
dir.create("not_ready/dorothea/networks", recursive = T, showWarnings = F)
X = dorothea::
X = X[EXPECTED_GRN_COLNAMES]
write_grn_by_subnetwork("dorothea", "all", grn_df = X)






