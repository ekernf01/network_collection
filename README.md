### Using the networks

This is a collection of pre-built gene regulatory networks. We offer R and Python code to quickly read and write tissue subnetworks and metadata from this collection. The loader code has minimal third-party dependencies (just `magrittr` and `arrow` in R and just `pandas` in python) and is tested inside of the benchmarking conda environment. See the benchmarking repo for more detail on that environment. 

In R:

```
source("R/load_networks.R")
# Set this to point to the "networks" folder adjacent to this README. 
options(GRN_PATH = "networks")
# What networks are available?
View(load_grn_metadata())
# What tissues do they cover, or how many?
list_subnetworks("gtex_rna") %>% head
lapply(load_grn_metadata()[[1]], list_subnetworks) %>% sapply(length)
# Show me the edges for a tissue. 
load_grn_by_subnetwork("gtex_rna", "Adipose_Subcutaneous.txt.gz") %>% head
# Show me the edges for all tissues.
iterate_within_grn("gtex_rna", load_grn_by_subnetwork) %>% sapply(dim)
# Count the edges for all tissues for all networks. Takes a while to run.
# iterate_over_grns(load_grn_by_subnetwork) %>% sapply(nrow)
```

In Python:

```
# Set this to point to the "load_networks" folder inside the "networks" folder adjacent to this README. 
sys.path.append('path/to/load_networks/') 
import load_networks
# Set this to point to the "networks" folder adjacent to this README. 
os.environ["GRN_PATH"] = "networks"
# What networks are available?
load_networks.load_grn_metadata()
# What tissues do they cover, or how many?
load_networks.list_subnetworks("gtex_rna")
[ load_networks.list_subnetworks(n)[0] for n in load_networks.load_grn_metadata()['name'] ]
# Show me the edges for a tissue. 
load_networks.load_grn_by_subnetwork("gtex_rna", "Adipose_Subcutaneous.csv.gz").head()
# Show me the edges for all tissues in one network.
[load_networks.load_grn_by_subnetwork("gtex_rna", n).shape for n in load_networks.list_subnetworks('gtex_rna') ]
```

### Installation 

This collection is not yet set up for deployment to non-Eric users. Main obstacles:

- The R code is loose scripts, not packages. 
- The Python code is not pip-installable or conda-installable. But it's in this repo, and you can point sys.path.append to it.
- The networks themselves are too big to put on GitHub. But they are on Patrick's AWS at s3://cahanlab/eric.kernfeld/eric_laptop/research/projects/perturbation_prediction/cell_type_knowledge_transfer/networks/.

### Storage format

Metadata are stored in `networks/published_networks.csv`. `networks` also contains a set of published GRN's stored uniformly. Each network is stored as `<source_name>/networks/<subnetwork_name>.csv.gz`. (UPDATE: now `.parquet`.) The basic format looks like this. Edge weight of -1 means the network is unweighted.

    regulator,target,edge_weight
    Pou5f1,Sox2,1
    ...

The data used to have some differences: gzipped/not, comma/tab delimited, target-first versus regulator-first. These have all been worked out but you may see out of date descriptions in the original README files, which are still included with some networks.

### Source data 

Source URL's are given in the metadata file. The setup is not fully automated, but large parts of it are, and scripts are given in the `setup` folder (start with `main.R`). A few notes on specific networks:

- For CellOracle, I installed the package, which ships with a default base network. Then there's a setup script that does the formatting and gzipping. No point-and-click downloads needed.
- Networks from FNTM and Humanbase were subsetted to retain only edges with posterior probability over 50%. See setup script; no point-and-click downloads needed.
- The least straightforward setup was for the CellNet networks, but it's also mostly automated. Starting from the [r package page](http://pcahan1.github.io/cellnetr/), I copied these files:

      cnProc_Hg1332_062414.R
      cnProc_mogene_062414.R
      cnProc_Hugene_062414.R
      cnProc_mouse4302_062414.R
    
  I installed `cellnetr` using [this tip](https://groups.google.com/forum/#!topic/cellnet_r/pXHt2J6ZH6I) to solve the `affy` version problem. I then read the data and exported GRN's using `cellnet.R` in the setup folder.

### Layout

- `load_networks`: Code to access these networks
- `networks`: Networks stored in triplet format
- `not_ready`: Dataset that we may in the future process into triplet format
- `setup`: code we used to assemble and format this collection.