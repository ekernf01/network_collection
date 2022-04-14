library("dplyr")
report_size = function(x) { cat("Edges: ", nrow(x)) ; return(x) }

#' Filter on >=50% posterior probability. 
#'
filter_file = function( filename_unfiltered, filename_filtered, min_posterior_prob = 0.5 ){
  chunked::read_chunkwise(filename_unfiltered, chunk_size=1e4, format = "table", header = F ) %>%
    dplyr::filter(V3 >= min_posterior_prob) %>%
    chunked::write_chunkwise(filename_filtered)
}

#' Obtain and filter fntm and humanbase networks. 
#'
retrieve_networks = function(source_name, base_url, suffix = "_top.gz"){
  destination = here("networks", source_name)
  dir.create.nice(file.path(destination, "networks"))
  tissues = 
    read.table(file.path(destination, "tissues.tsv"), sep = "\t", header = T) %>% 
    subset(download == "yes", select = "tissue", drop = T)
  files = tissues %>% gsub("\t", "", .) %>% gsub(" |-", "_", .) %>% tolower %>% paste0(suffix)
  do_one = function( name ) {
    message("\n", name); 
    filename_unfiltered = file.path(destination, "networks", name)
    filename_filtered = filename_unfiltered %>% gsub("(\\.gz)*$", "_filtered.txt", .)
    filename_filtered_gzipped = paste0(filename_filtered, ".gz")
    if(file.exists(filename_filtered_gzipped)){
      cat("Filtered file found. Finishing.\n")
      return()
    }
    curl::curl_download(url = paste0(base_url, name ),
                        destfile = filename_unfiltered,
                        quiet = F)
    system( paste0("gunzip ", filename_unfiltered) )
    filename_unfiltered %<>% gsub(".gz", "", .)
    filter_file(filename_unfiltered, filename_filtered)
    file.remove(filename_unfiltered)
    system( paste0("gzip ", filename_filtered) )
  }
  lapply( files, function(name) try(do_one(name)) )
}

retrieve_networks(source_name = "humanbase", base_url = "https://s3-us-west-2.amazonaws.com/humanbase/networks/")
retrieve_networks(source_name = "fntm",      base_url = "http://fntm.princeton.edu/static//networks10/")

#' Convert gene IDs from Entrez to symbol
#' 
#' 
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(annotate)
convert_ids = function(name, species){
  if(species=="human"){
    db = "org.Hs.eg.db"
  } else if(species=="mouse"){
    db = "org.Mm.eg.db"
  } else {
    stop("Sorry, only have human or mouse.\n")
  }
  filename_unfiltered = file.path(destination, "networks", name)
  filename_filtered = filename_unfiltered %>% gsub("(\\.gz)*$", "_filtered.txt.gz", .)
  filename_converted =  filename_filtered %>% gsub("(\\.gz)*$", "_converted.csv", .)
  filename_filtered %>%
    gsub("(.gz)*$", ".gz", .) %>%
    read.csv(stringsAsFactors = F) %>%
    dplyr::mutate(V1 = getSYMBOL(as.character(V1), data=db),
                  V2 = getSYMBOL(as.character(V2), data=db)) %>% 
    write.table(filename_converted, sep = ",", col.names = F, row.names = F)
  system(paste("gzip ", filename_converted))
}

convert_ids_all = function(source_name, species, suffix = "_top.gz"){
  destination = here("networks", source_name)
  dir.create.nice(file.path(destination, "networks"))
  tissues = 
    read.table(file.path(destination, "tissues.tsv"), sep = "\t", header = T) %>% 
    subset(download == "yes", select = "tissue", drop = T)
  files = tissues %>% gsub("\t", "", .) %>% gsub(" |-", "_", .) %>% tolower %>% paste0(suffix)
  lapply(files, convert_ids, species = species)
}

convert_ids_all(source_name = "humanbase", species = "human")
convert_ids_all(source_name = "fntm",    , species = "mouse")

