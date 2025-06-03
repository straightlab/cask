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
import os


def map_ag(ag_filepath, level, df, ci2_col, unq_col):
	""" 
	Maps ambiv group ID to repeat and class types. 
	Inputs: 
		ci2 and unq pickle file from CASK
		read level ambiv group assignments as a pandas dataframe (from loag_ag_dat function)
	Outputs:
		read level ambiv group assignments including repeat and class level information 
	"""
	df = df.copy()
	

	ci2_pickle = ag_filepath.split(".ci2_and_unq.ag.txt")[0] + ".ci2.aglut.pickle"
	unq_pickle = ag_filepath.split(".ci2_and_unq.ag.txt")[0] + ".unq.aglut.pickle"

	if level == "class":
		ag_dict_key = 'ag_to_chr_LUD_class'		
	else:
		ag_dict_key = 'ag_to_chr_LUD'
		
	with open(ci2_pickle, 'rb') as ci2p:
		ci2_ag_dict=pickle.load(ci2p)
		ci2agdict = ci2_ag_dict[ag_dict_key]
	with open(unq_pickle, 'rb') as unqp:
		unq_ag_dict=pickle.load(unqp)
		unqagdict = unq_ag_dict[ag_dict_key]

	#map the names of the repeat types to the ambiv groups
	df["ag_reptypes_ci2"]=df[ci2_col].map(ci2agdict)
	df["ag_reptypes_unq"]=df[unq_col].map(unqagdict)

	return df


def resolve_reptypes(row):
	""" 
	Compares repeat assignments according to ci2 k-mers (k-mers that were present at least 2 times in the full repeat database) and unq k-mers (only found once in the entire repeat database) to resolve ambivalence groups where possible.
	Example:
		A read that is assigned to an ambiv group encompassing {repeat_1, repeat_2, repeat_3} based on k-mers that were found more than once in the entire repeat database, and assigned {repeat_1} based on k-mers that were only found once (in this case they were found in repeat_1). As long as there are at least 2 unq k-mers found in the read, the read will be resolved to {repeat_1} assignment.
	Inputs: 
		Row from repeat/class level mapped ag data from load_ag_dat function (each row is one read)
	Outputs:
		1. Resolved repeat type list
		2. Type of resolution 
	Resolution Types:
		ci2z_unqres				No ci2 k-mers found (ag ID = 0) or repeat types assigned based on ci2 k-mers were conflicting (ag ID = -1). Resolved repeat types is based only on the unq k-mer assignment

		unqz_ci2res				No unq k-mers found (ag ID = 0) or repeat types assigned based on unq k-mers were conflicting (ag ID = -1). Resolved repeat types is based only on the ci2 k-mer assignment

		low_kcount_nores		No assignment (ag ID = 0) or conlficting assignment (ag ID = -1) from either unq or ci2 k-mers and too few k-mers from the other set (unq or ci2). No resolution could occur and returns NA

		same 					Assignment based on unq and ci2 k-mers is identical. Returns the assignment from unq k-mer assignment

		issubset				Types assigned based on unq k-mers is a subset of those assigned based on ci2 k-mers. Returns the assignment from unq k-mer assignment (as long as there are at least 2 unq k-mers present in the read). Example: ci2 = {repeat_1, repeat_2, repeat_3}, unq = {repeat_1}, resolved = {repeat_1}
		
		ci2_wins				Types assigned based on unq k-mers is not a subset of those assigned based on ci2 k-mers. There are more ci2 k-mers in the read than unq k-mers, so resolved repeat types = ci2 repeat types. Example: ci2 = {repeat_1, repeat_2, repeat_3},ci2_kcount = 8; unq = {repeat_6}, unq_kcount=2; resolved = {repeat_1, repeat_2, repeat_3}.

		unq_wins				Same as ci2_wins, but there were more unq k-mers than ci2 k-mers. Outputs repeat types based on unq k-mers

		conflicting				There are equal number of k-mers from unq and ci2 and the repeat types are not a subset of eachother. No resolution can occur. Outputs NA

	Note:
		Resolutions only occur if there are at least 2 k-mers on the side being used to resolve repeat types. This is to prevent the presence of a single k-mer (which could be the result of a sequencing error) from dominating the read assignment.
		CASK ambiv group ID = 0 --> No repeat k-mers were found in the read
		CASK ambiv group ID = -1 --> The possible repeat types for the read were conflicting. Example: 4 kmers found in the read are found in the repeat database (k1, k2, k3, k4). k1 and k2 are found in {repeat_1} k3 is found in {repeat_1, repeat_3, repeat_5}, and k4 is found in {repeat_5}. Since k1 and k2 can only be found in repeat_1 and k4 can only be found in repeat_5, CASK can not discern which repeat it comes from. This could be due to sequencing errors, reads overlapping multiple repeats, or genome differences (SNPs).

	"""

	ci2_reptypes_list, unq_reptypes_list = row['ag_reptypes_ci2'], row['ag_reptypes_unq']
	ci2kcount, unqkcount = row["ci2kcount_inread"], row['unqkcount_inread']

	if ci2_reptypes_list in [['0'], ['-1']] and unq_reptypes_list not in [['0'], ['-1']]:
		if unqkcount > 1:
			return unq_reptypes_list, "ci2z_unqres"
		else:
			return ["NA"], "low_kcount_nores"

	elif unq_reptypes_list in [['0'], ['-1']] and ci2_reptypes_list not in [['0'], ['-1']]:
		if ci2kcount > 1:
			return ci2_reptypes_list, "unqz_ci2res"
		else:
			return ["NA"], "low_kcount_nores"
	elif ci2_reptypes_list in [['0'], ['-1']] and unq_reptypes_list in [['0'], ['-1']]:
		return ["NA"], "both_z"
	else:
		ci2_set, unq_set = set(ci2_reptypes_list), set(unq_reptypes_list)
		if ci2_set == unq_set:
			return unq_reptypes_list, "same"
		elif unq_set.issubset(ci2_set):
			if unqkcount > 1:
				return unq_reptypes_list, "subset"
			elif ci2kcount > 1:
				return ci2_reptypes_list, "low_kcount_nores"
			else:
				return ["NA"], "low_kcount_nores"
		else:
			if ci2kcount > unqkcount and ci2kcount > 1:
				return ci2_reptypes_list, "ci2_wins"
			elif unqkcount > ci2kcount and unqkcount > 1:
				return unq_reptypes_list, "unq_wins"
			else:
				return ["NA"], "conflicting"




def load_ag_dat(ag_file, level, save=False):
	"""
	Loads read level ambiv group assignment file from CASK, maps the corresponding repeat types based on their IDs using the function map_ag, and resolves possible repeat types using the function resolve_reptypes.

	Inputs:
		1. Read level ambicv group assignment file_path from CASK (i.e. /path/to/dna.fastq.gz.ANNOTATED.ci2_and_unq.ag.txt)
		2. class or repeat type (level)
	Outputs:
		Dataframe with each row as an ag mapped, resolved read. This is also saved as a parquet file in the same directory as the CASK files (Input 1)
		Resolution stats file (how many reads were resolved based on each resolution type). Saved as .txt file in the same directory as above. 
	"""

	#cask_parent_dir = "/".join(ag_file.split("/")[0:-1]) + "/"
	out_file = ag_file.split(".ag.txt")[0] + "_" + level + ".resolved_agmapped_reads.parquet"
	#out_file = cask_parent_dir + ag_side + "_" + level + ".resolved_agmapped_reads.parquet"
	
	out_stats_file = ag_file.split(".ag.txt")[0] + "_" + level + ".resolved_agstats.txt" 
	#out_stats_file = cask_parent_dir + ag_side + "_" + "resolved_ag_stats_" + level + ".txt"

	if os.path.exists(out_file):
		print("Loading existing ag_dat:", out_file)
		#out_df = pd.read_csv(out_file, sep="\t", header=0, index_col=0)
		#out_df["ag_reptypes"] = out_df["ag_reptypes"].apply(ast.literal_eval)
		out_df = pd.read_parquet(out_file)
		out_df["ag_reptypes"] = out_df["ag_reptypes_str"].str.split(',')
		res_stats_df = pd.read_csv(out_stats_file, sep="\t", header=0, index_col=0)
	else:

		if level == "class":
			ci2_col = "class_aG_ci2"
			unq_col = "class_aG_unq"
		else:
			ci2_col = "reptype_ag_ci2"
			unq_col = "reptype_ag_unq"

		ag_df = pd.read_csv(ag_file, sep="\t", usecols=[0,1,2,3,4,5,6], names=['readid','reptype_ag_ci2','class_aG_ci2','ci2kcount_inread','reptype_ag_unq','class_aG_unq','unqkcount_inread'])
		level_df = ag_df[['readid',ci2_col,'ci2kcount_inread', unq_col,'unqkcount_inread']]
		mapped_ag_df = map_ag(ag_file, level, level_df, ci2_col, unq_col)
		print("Finished mapping ags")

		# Resolve reptypes based on both ci2 and unq
		mapped_ag_df[['ag_reptypes','res_type']] = mapped_ag_df.apply(resolve_reptypes, axis=1, result_type='expand')
		res_stats_df = mapped_ag_df.groupby(['res_type']).size()
		print("Finished resolving reptypes")
		
		
		mapped_ag_df['ag_reptypes_str'] = (mapped_ag_df['ag_reptypes'].transform(lambda x: ",".join(map(str,x))))

		#filter conflicts and zeros
		mapped_ag_df=mapped_ag_df[mapped_ag_df['ag_reptypes_str']!= "NA"]

		res_stats_df["Total_resolved"]=len(mapped_ag_df)

		out_df = mapped_ag_df.drop(columns=['res_type'])

		
		
		if save:
			out_df = out_df.drop(columns=['ag_reptypes'])
			out_df.to_parquet(out_file)
			#out_df['ag_reptypes']=out_df['ag_reptypes'].apply(lambda x: str(x))
			print("Saved counts as: ", out_file)

			
			res_stats_df.to_csv(out_stats_file, sep="\t", header=True)
			print("Saved stats as: ", out_stats_file)

	return out_df, res_stats_df




def count_rep_reads(res_agmapped_df):
	"""
	Counts the number of reads assigned to each ambiv group
	Inputs:
		output from load_ag_dat
	Output:
		Dataframe with the ambiv group  (a comma sep string of the repeat types possible for that amibv group) and the total number of reads assigned to that ambiv group
		{repeat_1,repeat_2: 4, repeat_5: 2, ...}
	"""
	ag_counts_df = res_agmapped_df.groupby('ag_reptypes_str').size().reset_index(name='count')

	return ag_counts_df


def count_rep_gbin_reads(res_agmapped_df, gbin_reads_file):
	"""
	Counts the number of reads assigned to each ambiv group as well as the number of reads aligned in each genomic bin. Reads that have a CASK assignment are excluded from gbin counts
	Inputs:
		output from load_ag_dat
		reads intersected with genome bins bed file
	Output:
		Dataframe with the ambiv group  (a comma sep string of the repeat types possible for that amibv group) and the total number of reads assigned to that ambiv group
		{repeat_1,repeat_2: 4, repeat_5: 2, ...}
	"""
	cask_readids = res_agmapped_df['readid'].tolist()
	gbin_df = pd.read_csv(gbin_reads_file, sep="\t", usecols=[3,9], names=['readid','gbin'])
	#print("gbin_df:", gbin_df)

	filt_gbin_df = gbin_df[~gbin_df['readid'].isin(cask_readids)]

	filt_gbin_df.rename(columns={'gbin':'ag_reptypes_str'}, inplace=True)
	#print("filt_gbin_df:", filt_gbin_df)


	gbin_counts_df = filt_gbin_df.groupby('ag_reptypes_str').size().reset_index(name='count')
	ag_counts_df = res_agmapped_df.groupby('ag_reptypes_str').size().reset_index(name='count')

	full_counts_df = pd.concat([ag_counts_df,gbin_counts_df])

	return full_counts_df

