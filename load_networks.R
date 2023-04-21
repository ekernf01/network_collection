library("magrittr") # for pipes

EXPECTED_GRN_COLNAMES = c("regulator", "target", "weight")
EXPECTED_OUTPUT_TYPES = c("tsv", "csv", "txt", "parquet")
options(GRN_PATH = 'networks')

#' Convert a function into a version of itself wrapped in a try block.
#'
wrapWithTry = function(INNERFUN){
  return(function(...){
    try(INNERFUN(...))
  })
}
  
#' Load source metadata
#' 
#' 
#' Returns a dataframe
load_grn_metadata = function( complete_only = T ){
  metadata_df = read.csv(file.path(getOption("GRN_PATH"), "published_networks.csv"), stringsAsFactors = F)
  rownames(metadata_df) = metadata_df[["name"]]
  if(complete_only){ metadata_df = subset(metadata_df, is_ready == "yes")}
  metadata_df
}

#' List all tissues available from a given source. 
#'
#' @param grn_name source to list tissues from
list_subnetworks = function(grn_name){
  list.files(file.path(getOption("GRN_PATH"), grn_name, "networks"))
}

#' Apply a function to all or many grns.
#' 
#' @param INNERFUN function accepting a single input corresponding to the name column in the source metadata.
#' 
iterate_over_grns = function( INNERFUN, omit = NULL, include = NULL, ...){
  X = load_grn_metadata(  complete_only = T )
  sources_use = X[["name"]] %>% setdiff(omit)
  if(!is.null(include)){
    sources_use %<>% intersect(include)
  }
  lapply(sources_use, iterate_within_grn, INNERFUN = wrapWithTry(INNERFUN), ...) %>% setNames(sources_use)
}

#' Apply a function to all subnetworks from a given source.
#' 
#' 
iterate_within_grn = function(grn_name, INNERFUN, test = F, omit = NULL, restrict_to_subnet = NULL, mc.cores = 14, ...){
  cat("\n\n", grn_name, "\n")
  available_subnetworks = list_subnetworks(grn_name)
  available_subnetworks %<>% setdiff(omit)
  if(!is.null(restrict_to_subnet)){
    available_subnetworks %<>% intersect(restrict_to_subnet)
  }
  if(length(available_subnetworks) == 0){
    warning(paste0("Skipping network ", grn_name, " due to zero overlap with allowed subnets."))
    return(NULL)
  }
  if(test){
    available_subnetworks = sample(available_subnetworks, 1)
  }
  parallel::mcmapply(wrapWithTry(INNERFUN),
                     grn_name, 
                     available_subnetworks, 
                     MoreArgs = list(...),
                     SIMPLIFY = F,
                     mc.cores = mc.cores) %>%
    setNames(available_subnetworks)
}

#' Load a given gene regulatory (sub-)network.
#' 
#' @param grn_name @param dge_name This function looks in file.path(getOption("GRN_PATH"), grn_name, "networks", subnetwork_name).
#' 
load_grn_by_subnetwork = function( grn_name, subnetwork_name, sort_by_weight = T, ... ){
  grn_location = file.path(getOption("GRN_PATH"), grn_name, "networks", subnetwork_name)
  if(is.null(grn_location)){
    stop("Error locating grn! Please report this error or check getOption('GRN_PATH').\n")
  }
  extension = tools::file_ext(gsub(".gz$", "", grn_location))
  if (extension %in% c("txt", "tsv")){
    X = read.table(grn_location, header = F, row.names = NULL, stringsAsFactors = F, ... )
  } else if (extension=="csv"){
    X = read.csv  (grn_location, header = F, row.names = NULL, stringsAsFactors = F, ... ) 
  } else if (extension=="parquet"){
    # TODO: handle .parquet.gz, which currently will stampede through uncaught.
    X = arrow::read_parquet(grn_location) 
  } else {
    stop(paste("Unknown format! Use one of: \n", paste(collapse = " ", EXPECTED_OUTPUT_TYPES), "\n (gzipped is ok for txt, csv, tsv)"))
  }
  
  # add score of -1 if missing
  if(ncol(X) == 2){
    X[[3]] = -1
  }
  
  # Fix meanings of columns
  first_column_meaning = load_grn_metadata()[grn_name, "first_column"]
  if(first_column_meaning %in% c("symmetric", EXPECTED_GRN_COLNAMES[[1]])){
    colnames(X) = EXPECTED_GRN_COLNAMES
  } else if(first_column_meaning == EXPECTED_GRN_COLNAMES[[2]]){
    colnames(X) = EXPECTED_GRN_COLNAMES[c(2, 1, 3)]
  } else {
    stop("Can't tell if target or regulator comes first! Fix the metadata.\n")
  }
  
  # Remove artifacts from setup
  if( all( unlist(X[1,]) == c("i", "j", "x") ) ){
    warning(c("In ", grn_name, " subnetwork ", subnetwork_name, ", removing a presumed header row i, j, x.\n"))
    X = X[-1, ]
  }
  
  if(sort_by_weight){
    X = X[order(X[[EXPECTED_GRN_COLNAMES[[3]]]], decreasing = T),]
  }
  
  return(X[EXPECTED_GRN_COLNAMES])
}

#' Make sure grn conforms to the expected structure.
#'
validate_grn = function( grn_name, subnetwork_name, grn_df = NULL, ... ){
  if(is.null(grn_df)){
    grn_df = load_grn_by_subnetwork( grn_name, subnetwork_name, ... )
  }
  get_num_unique = function(x) x %>% unique %>% length 
  num_unique = grn_df[1:2] %>% lapply(get_num_unique) 
  if(num_unique$regulator > num_unique$target){
    warning(paste0(grn_name, " ", subnetwork_name, " has more regulators than targets!\n") )
  }
  assertthat::assert_that(is.data.frame(grn_df))
  assertthat::are_equal(colnames(grn_df), EXPECTED_GRN_COLNAMES)
  return(T)
}

#' (Over)write a subnetwork file. Default behavior: write the same edges already present, but in a predictable format. 
#'
#' @param grn_df Dataframe with columns regulator, target, weight, in that order!
#'
write_grn_by_subnetwork = function( 
  grn_name,
  subnetwork_name, 
  grn_df = load_grn_by_subnetwork( grn_name, subnetwork_name ), 
  format_out = ".parquet"
){
  validate_grn(grn_name, subnetwork_name, grn_df)
  # The hardest thing here is if it's stored as target-first, but given to this function as 
  # regulator-first. Things will become out of sync. This issues a warning telling you 
  # to fix it by hand. 
  first_column_meaning = load_grn_metadata()[grn_name, "first_column"]
  if(first_column_meaning == EXPECTED_GRN_COLNAMES[[2]]){
    warning(paste0(
      "Network ",
      grn_name,
      " was previously stored as ",
      EXPECTED_GRN_COLNAMES[[2]], 
      "-first but is being written as ",
      EXPECTED_GRN_COLNAMES[[1]], 
      "-first. Please manually update metadata."
    ))
  }
  
  grn_location = file.path(getOption("GRN_PATH"), grn_name, "networks", subnetwork_name) 
  for(ot in EXPECTED_OUTPUT_TYPES){
    grn_location %<>% gsub("\\.gz", "", .)
    grn_location %<>% gsub(paste0("\\.", ot, "$"), format_out, .)
  }
  if(format_out==".csv.gz"){
    try({
      dir.create(dirname(grn_location), recursive = T, showWarnings = F)
      gz1 <- gzfile(grn_location, "w")
      write.table(grn_df, gz1, col.names = F, row.names = F, quote = F, sep=",")
      close(gz1)
    })
  } else if(format_out==".parquet"){
    try({
      dir.create(dirname(grn_location), recursive = T, showWarnings = F)
      arrow::write_parquet(grn_df, grn_location)
    })
  } else {
    stop("Unknown output format. Use '.csv.gz' or '.parquet.'")
  }
}

retrieve_targets_all = function(grn_name, subnetwork_name, query_regulators, target_only = T, ... ){
  cat(grn_name, " ", subnetwork_name, "\n")
  targets_predicted = load_grn_by_subnetwork(grn_name, subnetwork_name, ...) 
  if(target_only){
    targets_predicted %<>% subset(toupper(regulator) %in% toupper(query_regulators), select = "target", drop = T)
  } else {
    targets_predicted %<>% subset(toupper(regulator) %in% toupper(query_regulators))
  }
  return(targets_predicted)
}

retrieve_regulators_all = function(grn_name, subnetwork_name, query_targets, ... ){
  cat(grn_name, " ", subnetwork_name, "\n")
  regulators_predicted = load_grn_by_subnetwork(grn_name, subnetwork_name, ...) %>% 
    subset(toupper(target) == toupper(query_targets), select = "regulator", drop = T)
  return(regulators_predicted)
}

count_regulator_intersections = function(grn_name, subnetwork_name, query_targets, normalize = F, do_sort ){
  cat(grn_name, " ", subnetwork_name, "\n")
  # Accept either a character vector or a list of character vectors 
  if(is.atomic(query_targets)){
    query_targets = list(query_targets)
  }
  # Count intersection of each tf with one query
  get_predicted_regulators = function(current_targets){
    X = load_grn_by_subnetwork(grn_name, subnetwork_name)
    regulators_predicted = X %>% 
      subset(toupper(target) %in% toupper(current_targets), select = "regulator", drop = T) %>% 
      table 
    if(normalize){
      totals = X[["regulator"]] %>% table
      do_one = function(x, y, z) fisher.test(matrix(c(29998, x, y, z), nrow = 2))$p
      regulators_predicted =
        mapply( do_one,
                x = totals[names(regulators_predicted)],
                y = length(current_targets),
                z = regulators_predicted,
                SIMPLIFY = T ) %>% log10 %>% multiply_by(-1)
                                    
    }
    if(do_sort){
      regulators_predicted %<>% sort(decreasing = T)
    }
    regulators_predicted
  }
  # process all queries
  lapply(query_targets, get_predicted_regulators)
}

# Check setup
check_networks_location = function(){
  tryCatch(
    invisible(load_grn_metadata()), 
    error = function(e) {
      message("Cannot find network metadata. Set options('GRN_PATH') to point to the folder containing 'published_networks.csv'.\nExample:\n    options(GRN_PATH = 'path/to/networks')")
    }
  )
}
check_networks_location()
