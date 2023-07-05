This is a collection of pre-built gene regulatory networks, accompanied by the code used to acquire and clean the data. This part of our [benchmarking project](https://github.com/ekernf01/perturbation_benchmarking).

### Using the networks

For the Python API, consult the [companion package](https://github.com/ekernf01/load_networks).

There is also an unfinished R API included as a script in this folder. Usage:

```r
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

### Installation 

The networks themselves are too big to put on GitHub. Download from Zenodo (DOI: 10.5281/zenodo.8071809).

### Storage format

Metadata are stored in `networks/published_networks.csv`. `networks` also contains a set of published GRN's stored uniformly. Each network is stored as `<source_name>/networks/<subnetwork_name>.csv.gz`. (UPDATE: now `.parquet`.) The basic format looks like this. Edge weight of -1 means the network is unweighted.

    regulator,target,weight
    Pou5f1,Sox2,1
    ...

The data used to have some differences: gzipped/not, comma/tab delimited, target-first versus regulator-first. These have all been worked out but you may see out of date descriptions in the original README files, which are still included with some networks.

### Source data 

Source URL's, citations, and descriptions are given in the [metadata file](https://github.com/ekernf01/network_collection/blob/main/networks/published_networks.csv). The setup is not fully automated, but large parts of it are, and scripts are given in the `setup` folder (start with `main.R`). A few notes on specific networks:

- For CellOracle, I installed the package, which ships with a default base network. Then there's a setup script that does the formatting and gzipping. No point-and-click downloads needed.
- Networks from FNTM and Humanbase were subsetted to retain only edges with posterior probability over 50%. See setup script; no point-and-click downloads needed.
- The least straightforward setup was for the CellNet networks, but it's also mostly automated. Starting from the [r package page](http://pcahan1.github.io/cellnetr/), I copied these files:

      cnProc_Hg1332_062414.R
      cnProc_mogene_062414.R
      cnProc_Hugene_062414.R
      cnProc_mouse4302_062414.R
    
  I installed `cellnetr` using [this tip](https://groups.google.com/forum/#!topic/cellnet_r/pXHt2J6ZH6I) to solve the `affy` version problem. I then read the data and exported GRN's using `cellnet.R` in the setup folder.

### Layout

- `networks`: Networks stored in triplet format
- `not_ready`: Dataset that we may in the future process into triplet format
- `setup`: code we used to assemble and format this collection.
