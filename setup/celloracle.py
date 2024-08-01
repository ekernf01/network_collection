# This script should work in any environment where celloracle is installed. 
import os
import celloracle as co
import pandas as pd
os.environ["GRN_PATH"] = "../network_collection/networks"
def fix_format(x, net, subnet):
    x.drop(["peak_id"], axis = 1, inplace=True)
    x = pd.melt(x, id_vars = "gene_short_name")
    x.set_axis(["target", "regulator", "weight"], axis=1, inplace=True)
    x = x.loc[x['weight'] == 1]
    x["weight"] = -1
    x = x[["regulator", "target", "weight"]]
    print(x.head())
    write_to = os.path.join(os.environ["GRN_PATH"], net, "networks", subnet)
    try:
        os.makedirs(os.path.dirname(write_to))
    except:
        pass
    x.to_parquet(write_to)

fix_format(x = co.data.load_human_promoter_base_GRN(), net = "celloracle_human", subnet = "human_promoters.parquet")
fix_format(x = co.data.load_mouse_promoter_base_GRN(), net = "celloracle_mouse", subnet = "all_promoters.parquet")
fix_format(x = co.data.load_zebrafish_promoter_base_GRN(), net = "celloracle_zebrafish", subnet = "all_promoters.parquet")
