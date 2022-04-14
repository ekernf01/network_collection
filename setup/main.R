### Daily setup
{
  library("Matrix")
  library("magrittr")
  library("parallel")
  source( "../R/load_networks.R" )
  dr_here()
  PROJ_DIR = here()
  getwd()
  here()
}

# Data setup is not fully automated, but these scripts do a lot of the work.
# We hope they will increase the transparency of our project despite not
# being fully automated. 
source( "cellnet.R")
source( "chea.R")
source( "fntm_humanbase.R")
source( "gtex_rna.R")
source( "reg_net_dot_org.R")
