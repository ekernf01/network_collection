# This is my attempt to download a dataset of TF-target interactions that the authors of Mogrify describe as MARA.
# The important point here is to try to give a fair evaluation of this dataset in terms of the value it 
# contributes to Mogrify's cell engineering predictions by its ability to reflect causal structure
# of transcriptional control. If we overexpress some TF's nominated by Mogrify, to what extent do
# the effects propagate along the MARA network as opposed to other paths?
#
# Rackham, O. J., Firas, J., Fang, H., Oates, M. E., Holmes, M. L., Knaupp, A. S., ... & Gough, J. (2016). A 
# predictive computational framework for direct reprogramming between human cell types. Nature genetics, 48(3), 331-335.
# 
# Confusion arises on a few points.
# MARA is not a dataset; it is a computational method. It was first deployed on FANTOM4 data to yield
# a dataset of TF-promoter interactions. 
# 
# Center, R. O. S., & FANTOM Consortium. (2009). The transcriptional network that controls growth arrest 
# and differentiation in a human myeloid leukemia cell line. Nature genetics, 41(5), 553.
# 
# I think this is that dataset but I am not sure. 
# 
# You'll need these packages.
# BiocManager::install(c("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg18.knownGene"))

dir.create("not_ready/FANTOM4/networks", recursive = T, showWarnings = F)
curl::curl_download("https://swissregulon.unibas.ch/data/fantom4/REFSEQ_sites_Z.gff.gz", "not_ready/FANTOM4/networks/MARA_motif_to_promoter.gff.gz")
fantom4_network = rtracklayer::readGFFAsGRanges(gzfile("not_ready/FANTOM4/networks/MARA_motif_to_promoter.gff.gz"))

# Most regulators look reasonable, except the arabidopsis gene Agamous. 
# Clearly that's not what's binding those motifs. 
fantom4_network$ALIAS %>% table %>% View

# Determine a gene symbol or set of genes symbols for each motif.
# I use the Entrez id instead of the provided alias, but they are very similar. 
library(org.Hs.eg.db)
converter = 
  select(org.Hs.eg.db,
         keys = fantom4_network$ENTREZ_ID %>% unique %>% sapply(strsplit, ",") %>% unlist %>% unique, 
         columns = c("ENTREZID","SYMBOL"))
fantom4_network %<>% merge(converter, by.x = "ENTREZ_ID", by.y = "ENTREZID")
fantom4_network %<>% dplyr::rename(regulator_symbol = SYMBOL)
fantom4_network %<>% as("GRanges")
# To see how similar, sort this co-occurrence table by frequency.
# table(fantom4_network$regulator_symbol, fantom4_network$ALIAS) %>% View

# Determine a gene symbol for each target promoter.
library(TxDb.Hsapiens.UCSC.hg18.knownGene) # Promoters from swissregulon are in hg18 coords
txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
hg18_transcripts = genes(txdb)
hits = GenomicRanges::nearest(fantom4_network, hg18_transcripts)
fantom4_network$target_gene_id = hg18_transcripts$gene_id[hits]
converter = 
  select(org.Hs.eg.db,
         keys = fantom4_network$target_gene_id %>% unique, 
         columns = c("ENTREZID","SYMBOL"))
fantom4_network %<>% merge(converter, by.x = "target_gene_id", by.y = "ENTREZID")
fantom4_network %<>% dplyr::rename(target_symbol = SYMBOL)

write_grn_by_subnetwork(
  grn_name = "MARA_FANTOM4", 
  subnetwork_name = "THP-1.csv.gz", 
  grn_df = 
    fantom4_network[c("regulator_symbol", "target_symbol")] %>% 
    set_colnames(c("regulator", "target"))
)

# I am also not sure how to match promoters to genes or how to match motifs to TF's.
# MARA author Erik van Nimwegen says (answering someone else, on Twitter) this is hard.
# https://twitter.com/NimwegenLab/status/1514697519560859659
# The Mogrify authors give no description of their process. 

