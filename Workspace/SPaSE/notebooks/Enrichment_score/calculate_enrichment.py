import pandas as pd
import numpy as np
import tqdm
import glob
import sys
from enrichment_scores import score_abid_pval_and_enrichment
import os

option_idx =  int(sys.argv[1])

options_df = pd.read_csv("options_ot_updated.csv", index_col = 0)
fno_idx = options_df.loc[option_idx, "fno"]
pno_idx = options_df.loc[option_idx, "pno"]

# load files
files = np.sort(list(set(glob.glob("./adata*.csv")) \
             - set(glob.glob("./*spatial*"))))

# load programs
sheet_url = "https://docs.google.com/spreadsheets/d/1N5Hlo6ASSnrpV2GoZHSKPhRX8j_oFOYSFhcgyEgrlpo/edit#gid=0"
sheet_url = sheet_url.replace('/edit#gid=', '/export?format=csv&gid=')
programs_df = pd.read_excel("ot_gene_lists.xlsx", sheet_name="programs")
programs = np.sort(programs_df["Program"].unique())

fname = files[fno_idx]
pname = programs[pno_idx]   

print(fname.split("/")[-1], pname)

df = pd.read_csv(fname, index_col=0)

nbins = 50
obs_avg = df.mean(axis=0)
n_items = int(np.round(len(obs_avg) / (nbins - 1)))
obs_cut = (obs_avg.rank(method='min') // n_items).astype("int")

results_pval = pd.DataFrame(index = df.index)
results_enrichment = pd.DataFrame(index = df.index)


program_df = programs_df[programs_df["Program"] == pname]
program_genes = list(program_df["Gene"].values)

pval, enrichment = score_abid_pval_and_enrichment(df, program_genes, 50, 1000, obs_cut=obs_cut)
results_pval[pname] = pval
results_enrichment[pname] = enrichment

os.makedirs(os.path.dirname(f"./results_pval/{fno_idx}_{pno_idx}.csv"), exist_ok=True)
os.makedirs(os.path.dirname(f"./results_enrichment/{fno_idx}_{pno_idx}.csv"), exist_ok=True)

results_pval.to_csv(f"./results_pval/{fno_idx}_{pno_idx}.csv")
results_enrichment.to_csv(f"./results_enrichment/{fno_idx}_{pno_idx}.csv")
