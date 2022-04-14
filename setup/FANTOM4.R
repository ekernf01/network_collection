dir.create("../not_ready/FANTOM4/networks", recursive = T, showWarnings = F)
curl::curl_download("https://swissregulon.unibas.ch/data/fantom4/REFSEQ_sites_Z.gff.gz", "../not_ready/FANTOM4/networks/MARA.txt.gz")
network = read.delim("../not_ready/FANTOM4/networks/MARA.txt.gz", comment.char = "#", header = F)
head(network)
