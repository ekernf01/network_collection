import pandas as pd 
import os 
EXPECTED_GRN_COLNAMES = ["regulator", "target", "weight"]

def load_grn_metadata( complete_only = True ):
  """Load metadata on GRN sources."""
  metadata_df = pd.read_csv(os.path.join(os.environ["GRN_PATH"], "published_networks.csv"))
  metadata_df.index = metadata_df["name"]
  if complete_only: 
    metadata_df = metadata_df.loc[metadata_df['is_ready'] == "yes",:]
  return metadata_df

def list_subnetworks(grn_name):
  """
  List all tissues available from a given source. 
  Parameters:
        - grn_name (string) source to list tissues from
  """
  return [f for f in os.listdir(os.path.join(os.environ["GRN_PATH"], grn_name, "networks")) if not f.startswith('.')]


def load_grn_by_subnetwork( grn_name, subnetwork_name, sort_by_weight = True ):
  """
  Load a given gene regulatory (sub-)network.
  Parameters:
        - grn_name (string): source to list tissues from
  """
  grn_location = os.path.join(os.environ["GRN_PATH"], grn_name, "networks", subnetwork_name)
  if not os.path.exists(grn_location):
    raise ValueError("Error locating grn! Please report this error or check os.environ['GRN_PATH'].\n")
  
  X = pd.read_parquet( grn_location ) 
  
  # add score of -1 if missing
  if(X.shape[1] == 2):
    X["weight"] = -1
  
  X.set_axis(EXPECTED_GRN_COLNAMES, axis = 1, inplace = True)
  # This saves mem for fat networks.
  X[EXPECTED_GRN_COLNAMES[0]].astype("category")
  X[EXPECTED_GRN_COLNAMES[1]].astype("category")
  return X 

def load_grn_all_subnetworks(grn_name):
  return pd.concat(
      [load_grn_by_subnetwork(grn_name, s) for s in list_subnetworks(grn_name)]
   )
    
def validate_grn(
  grn_name,
  subnetwork_name, 
  grn_df = None
):
  """Make sure grn conforms to the expected structure."""
  if grn_df is None:
      grn_df = load_grn_by_subnetwork( grn_name, subnetwork_name )
  if len(grn_df[:, "regulator"].unique()) > len(grn_df[:, "target"].unique()):
    warning(paste0(rn_name, " ", subnetwork_name, " has more regulators than targets!\n") )
  assert type(grn_df) == pd.DataFrame
  assert al( [grn_df.columns[i] == EXPECTED_GRN_COLNAMES[i] for i in range(3)])
  return True



