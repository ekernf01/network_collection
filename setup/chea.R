X = read.table(here("networks", "chea", "gene_attribute_matrix.txt"), sep = "\t", header = T, comment.char = "")
X = X[-(1:2), -(2:3)]
rownames(X) = X$X.
X = X[,-1 ]
X %<>% as.matrix
X %<>% Matrix(sparse = T)
corner(X, 5, 5)
dim(X) # 21220 genes, 199 TF's
MatrixToTripletNamed = function (Y){
  assertthat::are_equal(class(Y), "dgCMatrix")
  X = summary(Y)
  X$i = rownames(Y)[X$i]
  X$j = colnames(Y)[X$j]
  X
}
X %>% MatrixToTripletNamed %>% 
  set_colnames(c("target", "regulator", "weight")) %>%
  extract(EXPECTED_GRN_COLNAMES) %>%
  write_grn_by_subnetwork( "chea", "all", grn_df = X )  
