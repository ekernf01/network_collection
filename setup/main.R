### Daily setup
{
  library("Matrix")
  library("magrittr")
  library("parallel")
  PROJ_DIR = "~/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/networks"
  setwd(PROJ_DIR)
  source( "../load_networks/load_networks.R" )
  Sys.setenv(GRN_PATH = "networks")
  check_networks_location()
}

# Data setup is not fully automated, but these scripts do a lot of the work.
# We hope they will increase the transparency of our project despite it not
# being fully automated. 
source( "setup/cellnet.R")
source( "setup/chea.R")
source( "setup/fntm_humanbase.R")
source( "setup/gtex_rna.R")
source( "setup/reg_net_dot_org.R")
source( "setup/FANTOM4.R")
source( "setup/ananse.R")
