# Here we acquire the GRN underlying the tool IRENE. 
# 
# Jung, S., Appleton, E., Ali, M., Church, G. M., & Del Sol, A. (2021). 
# A computer-guided design tool to increase the efficiency of cellular conversions. Nature communications, 12(1), 1-12.
# 
#
# Cistromedb appears to be used as a fundamental component of IRENE, but without being cited. 
# (They instead cite chip-atlas.)
# The source code for IRENE is at https://github.com/saschajung/IRENE, and in commit 17562742c, the 
# README describes a file 
#
#     Cistrome_ChIPseq_sorted.bed: All human TF ChIP-seq data taken from Cistrome browser.
#
# So, this is the file we use here.
# 
# Unlike Cistrome, IRENE uses GeneHancer to pair enhancers with promoters. 
# 
# The Cistrome download is done using the URL described in the metadata.
# It can't be automated according to their terms of service. 
cistrome = read.delim("???.bed")
# The GeneHancer download is obtained directly from IRENE:
curl::curl_download("???")