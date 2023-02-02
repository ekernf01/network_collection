import pandas as pd 
import duckdb
import os 
import gc
import numpy as np

EXPECTED_GRN_COLNAMES = ["regulator", "target", "weight"]


class LightNetwork:
  """
  The LightNetwork class is a lightweight, read-only interface to our database of weighted TF-target interactions.
  It supports just a few simple operations -- for now, just querying with a target and returning
  relevant TF's and weights. 
  Under the hood, it uses duckdb to query parquet files, so data are not loaded into RAM. 
  """
  def __init__(
    self,
    grn_name: str = None,
    subnetwork_names: list = [],
    files:list = [],
    df: pd.DataFrame = None,
  ):
    """Create a Network object. 

    Args:
        grn_name (str, optional): source of networks, e.g. "celloracle_human"
        subnetwork_name (list, optional): filename for individual network (usually networks are grouped by cell/tissue type).
        files (list, optional): List of absolute paths to parquet files. 
            These files will be queried in addition to the stuff already specified.
        df (pd.DataFrame, optional): This DF will be queried in addition to the stuff already specified.
    """
    if grn_name is None and len(subnetwork_names) > 0:
      raise ValueError("Cannot find subnetworks unless 'grn_name' is specified.")
    if grn_name is not None:
      if len(subnetwork_names) == 0:
        subnetwork_names = list_subnetworks(grn_name)
      files = files + [os.path.join(os.environ["GRN_PATH"], grn_name, "networks", s) for s in subnetwork_names]
    if len(files) == 0 and df is None:
      raise ValueError("Please provide at least one of 'grn_name' or 'files' of 'df'.")
    if df is not None:
      assert all(df.columns[0:3] == ['regulator', 'target', 'weight'] ), "Columns must have names 'regulator', 'target', 'weight'"
      assert len(df.columns)==3 or df.columns[3]=='cell_type', "If there are 4 columns, the 4th must be 'cell_type'."
    self.files = files
    self.df = df
    return

  def save(self, filename) -> None:
    if not filename.lower().endswith(".parquet"):
      raise ValueError("Filename must end with .parquet.")
    self.get_all().to_parquet(filename)
    return

  def get_all(self):
    results_from_parquet = pd.DataFrame()
    results_from_memory  = pd.DataFrame()
    if len(self.files) > 0:
      files_formatted = [f"'{file}'" for file in self.files]
      results_from_parquet = duckdb.query(
        " UNION ".join(
            [
              f"SELECT * FROM {file}" for file in files_formatted
            ]
          )
      ).df()
    if not self.df is None:
      con = duckdb.connect()
      df=self.df
      results_from_memory = con.execute(f"SELECT * FROM df").df()
    return pd.concat([results_from_parquet, results_from_memory])

  def get_regulators(self, target: str) -> pd.DataFrame: 
    """Return all records having a given target gene.

    Args:
        target (str): A target gene present in this network. 

    Returns:
        pd.DataFrame: All records with the given target gene
    """
    results_from_parquet = pd.DataFrame()
    results_from_memory  = pd.DataFrame()
    if len(self.files) > 0:
      files_formatted = [f"'{file}'" for file in self.files]
      results_from_parquet = duckdb.query(
        " UNION ".join(
            [
              f"SELECT * FROM {file} WHERE target = '{target}'" for file in files_formatted
            ]
          )
      ).df()
    if not self.df is None:
      con = duckdb.connect()
      df=self.df
      results_from_memory = con.execute(f"SELECT * FROM df WHERE target = '{target}'").df()
    return pd.concat([results_from_parquet, results_from_memory])
    
  def get_all_regulators(self) -> set: 
    """Return a set of all regulators

    Returns:
        set: All distinct regulators
    """
    return self.get_all_one_field("regulator")

  def get_all_one_field(self, field: str) -> set: 
    """Return a set of all regulators, or all targets, or all cell types.

    Returns:
        set: All distinct values listed in a given column from this network
    """
    assert field in {"regulator", "target", "cell_type"}, " Can only get unique vals for 'regulator', 'target', or 'cell_type' "
    results_from_parquet = pd.DataFrame(columns = [field])
    results_from_memory  = pd.DataFrame(columns = [field])
    if len(self.files) > 0:
      files_formatted = [f"'{file}'" for file in self.files]
      files_formatted = [
        file for file in files_formatted 
        if field in set(duckdb.query(f"SELECT * FROM {file} WHERE 1=0").df().columns)
      ]
      if len(files_formatted) > 0:
        results_from_parquet = duckdb.query(
          " UNION ".join(
              [
                f"SELECT DISTINCT {field} FROM {file}" for file in files_formatted
              ]
            )
        ).df()
    if not self.df is None and field in set(self.df.columns):
      con = duckdb.connect()
      df=self.df
      results_from_memory = con.execute(f"SELECT DISTINCT {field} FROM df").df()
    return set(pd.concat([results_from_parquet, results_from_memory])[field].unique())
    
  def get_num_edges(self) -> int: 
    """Return the number of available regulator-target connections.
    This assumes no overlap between edges recorded in 'files' and stuff provided in 'df'.

    Returns:
        int
    """
    results_from_parquet = 0
    results_from_memory  = 0
    if len(self.files) > 0:
      files_formatted = [f"'{file}'" for file in self.files]
      results_from_parquet = duckdb.query(
        "SELECT COUNT(*) FROM " + \
        "( " + \
          " UNION ".join(
              [
                f"SELECT * FROM {file}" for file in files_formatted
              ]
            ) + \
        " )"
      ).df().iloc[0,0]
    if not self.df is None:
      con = duckdb.connect()
      df=self.df
      results_from_memory = con.execute("SELECT COUNT(*) FROM df").df().iloc[0,0]
    return int(results_from_parquet + results_from_memory)
    
def load_grn_metadata( complete_only = True ):
  """Load metadata on GRN sources."""
  metadata_df = pd.read_csv(os.path.join(os.environ["GRN_PATH"], "published_networks.csv"))
  metadata_df.index = metadata_df["name"]
  if complete_only: 
    metadata_df = metadata_df.loc[metadata_df['is_ready'] == "yes",:]
  return metadata_df

def list_subnetworks(grn_name: str):
  """
  List all tissues available from a given source. 
  Parameters:
        - grn_name (string) source to list tissues from
  Return value: pd.DataFrame
  """
  return [f for f in os.listdir(os.path.join(os.environ["GRN_PATH"], grn_name, "networks")) if not f.startswith('.')]

def load_grn_by_subnetwork( grn_name: str, subnetwork_name: str ):
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
  
  X.set_axis(EXPECTED_GRN_COLNAMES, axis = 1, copy = False)
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
    raise ValueError("".join([grn_name, " ", subnetwork_name, " has more regulators than targets!\n"]) )
  assert type(grn_df) == pd.DataFrame
  assert al( [grn_df.columns[i] == EXPECTED_GRN_COLNAMES[i] for i in range(3)])
  return True




def networkEdgesToMatrix(networkEdges, regulatorColumn=0, targetColumn=1):
    """Reformat a network from a two-column dataframe to the way that celloracle needs its input."""
    X = pd.crosstab(networkEdges.iloc[:,targetColumn], networkEdges.iloc[:,regulatorColumn])
    del networkEdges
    gc.collect()
    X = 1.0*(X > 0)
    X = X.rename_axis('gene_short_name').reset_index()
    X = X.rename_axis('peak_id').reset_index()
    X = makeNetworkSparse(X, 0.0)
    gc.collect()
    return X

humanTFs = pd.read_csv("../accessory_data/humanTFs.csv")

def pivotNetworkWideToLong(network_wide: pd.DataFrame):
    """Convert from CellOracle's preferred format to a triplet format

    Args:
        network_wide (pd.DataFrame): GRN structure in CellOracle's usual format
    """
    network_long = pd.concat([
        pd.DataFrame({
            "regulator": tf,
            "target": network_wide.loc[network_wide[tf]==1, "gene_short_name"],
            "weight": 1,
        })
        for tf in network_wide.columns[2:]
    ])
    return network_long

def makeRandomNetwork(targetGenes, density = 0, seed = 0, TFs = humanTFs['HGNC symbol'] ):
    """Generate a random network formatted the way that celloracle needs its input."""
    np.random.seed(seed)
    X = pd.DataFrame(
            np.random.binomial(
                n = 1, 
                p=density,
                size=(
                    len(targetGenes), 
                    len(TFs)
                )
            ),
            columns = TFs, 
            index = targetGenes
        )
    X.rename_axis('gene_short_name', inplace=True)
    X.reset_index(inplace=True)
    X.rename_axis('peak_id', inplace=True)
    X.reset_index(inplace=True)
    # CellOracle's preferred format wastes gobs of memory unless you sparsify.
    X = makeNetworkSparse(X, round(density))
    gc.collect()
    return X


def makeNetworkSparse(X, defaultValue):
    """Save memory by making a sparse representation of a base network"""
    X.iloc[:,2:] = X.iloc[:,2:].astype(pd.SparseDtype("float", defaultValue))
    return X


def makeNetworkDense(X):
    """Undo makeNetworkSparse"""
    X.iloc[:, 2:] = np.array(X.iloc[:, 2:])   #undo sparse representation         
    return X

