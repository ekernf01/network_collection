from celloracle import motif_analysis as ma
from celloracle.motif_analysis.process_bed_file import df_to_list_peakstr as tolist
import sys
import pandas as pd
from gimmemotifs.motif import default_motifs
import genomepy
import pereggrn_networks

# Genome information
ref_genome = "hg38"
if not ma.is_genome_installed(ref_genome=ref_genome, genomes_dir=None):
    genomepy.install_genome(name=ref_genome, provider="UCSC", genomes_dir=None)

# Peaks
peaks = pd.read_csv(sys.argv[1], sep = "\t", header=None)[[0, 1, 2]]
peaks.columns = ["chr", "start", "end"]
peaks["start"] = peaks["start"].astype("int")
peaks["end"] = peaks["end"].astype("int")
peaks["seqname"] = peaks["chr"] + "_" + peaks.start.astype("str") + "_" + peaks.end.astype("str")
# By default, CellOracle includes only promoters.
# https://github.com/morris-lab/CellOracle/blob/master/docs/notebooks/01_ATAC-seq_data_processing/option2_Bulk_ATAC-seq_data/01_preprocess_Bulk_ATAC_seq_peak_data.ipynb
# https://github.com/morris-lab/CellOracle/blob/master/celloracle/motif_analysis/tss_annotation.py
# We will widen the peaks by 10k radius because we want to also include nearby enhancers.
peaks["start"] = peaks["start"] - 10000
peaks["end"] = peaks["end"] + 10000
peaks["start"] = [max(0, s) for s in peaks["start"]]
peaks = ma.get_tss_info(peak_str_list=tolist(peaks), ref_genome=ref_genome)
# Undo the widening before doing the motif scan.
peaks["start"] = peaks["start"] + 10000
peaks["end"] = peaks["end"] - 10000

# fix formatting
peaks = pd.DataFrame({"peak_id": tolist(peaks), "gene_short_name": peaks.gene_short_name.values})
peaks = peaks.reset_index(drop=True)
peaks = ma.check_peak_format(peaks, ref_genome, genomes_dir=None)

# motif scan & filtering
motifs = default_motifs()
tfi = ma.TFinfo(peak_data_frame=peaks,
                ref_genome=ref_genome,
                genomes_dir=None)
tfi.scan(motifs=motifs)
tfi.scanned_df.head()
tfi.reset_filtering()
tfi.filter_motifs_by_score(threshold=10)
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

# convert format
pereggrn_networks.pivotNetworkWideToLong(tfi.to_dataframe()).to_parquet(sys.argv[2])