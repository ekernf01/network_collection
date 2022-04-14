library("cellnetr", lib.loc = "~/Desktop/programs/cellnet_packages/")

convert_cellnet_grn = function(name, file){
  cnProc <- cellnetr::utils_loadObject(file);
  dir.create.nice(here("data", name, "full"))
  dir.create.nice(here( "networks", name ) )
  for( tissue in names(cnProc$igGRNs)){
    Y = igraph::as_adj( igraph::upgrade_graph( cnProc$igGRNs[[tissue]] ) )
    colnames(X) = c("regulator", "target", "weight")
    write_grn_by_subnetwork( name, tissue, grn_df = X )
  }
}


convert_cellnet_grn(name = "cellnet_human_Hg1332", file = "~/Downloads/cnProc_Hg1332_062414.R")
convert_cellnet_grn(name = "cellnet_mouse_4302",   file = "~/Downloads/cnProc_mouse4302_062414.R")
convert_cellnet_grn(name = "cellnet_human_Hugene", file = "~/Downloads/cnProc_Hugene_062414.R")
convert_cellnet_grn(name = "cellnet_mouse_mogene", file = "~/Downloads/cnProc_mogene_062414.R")
