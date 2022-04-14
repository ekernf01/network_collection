# This assumes os.environ["GRN_PATH"] points to the networks folder adjacent to the README. 
import celloracle as co
import pandas as pd
x = co.data.load_human_promoter_base_GRN()
x.drop(["peak_id"], axis = 1, inplace=True)
x = pd.melt(x, id_vars = "gene_short_name")
x.set_axis(["target", "regulator", "is_present"], axis=1, inplace=True)
x = x.loc[x['is_present'] == 1]
x.is_present *= -1
write_to = os.path.join(os.environ["GRN_PATH"], "celloracle_human", "networks", "human_promoters.csv")
try:
    os.makedirs(os.path.dirname(write_to))
except:
    pass
x.to_csv(write_to, header = False, index = False)
os.system("gzip " + write_to)
