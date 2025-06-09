from map_count_ags import *
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import linregress
import ast
import seaborn as sns
from pygam import LinearGAM, s, f
import scipy.stats as stats
import sys



levels=["reptype","class"]


counts_df_dict={}
res_stats_df_dict={}

#set the maximum number of chromosomes an ambivalence group may encompass. All others are not included.
max_ag_chrs = 5

#command line input of ag annotated reads file from cask_make_tally.sh (example: /path/to/my/readsdir/cen_kmatch/k25.minusBKG/myinputfastq.fastq.gz.ANNOTATED.ci2_and_unq.ag.txt)
i = sys.argv[1]


print("Starting:", i)
out_counts_file = i.split(".fastq.gz.ANNOTATED.ci2_and_unq.ag.txt")[0] + ".ag.counts.txt"


for l in levels:
	print("Handling ags for:", i, l)
	out_counts_file = i.split(".fastq.gz.ANNOTATED.ci2_and_unq.ag.txt")[0] + "." + l + ".ag.counts.txt"
	resolved_ag_df, res_stats_df  = load_ag_dat(i, l, save=True)

	#count reads for each  ag:
	ag_counts_df = count_rep_reads(resolved_ag_df)
	ag_counts_df.to_csv(out_counts_file, sep="\t")

	print("Done counting reads. Saved counts as ", out_counts_file)


#RUN AS:
#python run_agmap_count.py /path/to/my_input_fastq_file.fastq.gz.ANNOTATED.ci2_and_unq.ag.txt


	