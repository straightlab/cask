import pandas as pd
import numpy as np
import sys

#create a genome bin bedfile:
#run as: python generate_gbins.py my_chrom_sizes_file output_file_prefix
#saves a file for each genome bin size in 'bin_size_dict' as output_file_prefix + bin_size_key + "_genome_bins.bed"



chrom_sizes_file = sys.argv[1]
#chrom_sizes_file="/oak/stanford/groups/astraigh/T2T/CURRENT/T2T-CHM13v2.0-hs1_jan2022/hs1.chromosomes.bed"
output_filename_prefix = sys.argv[2]
#output_dir="/home/groups/astraigh/kelsey/chm13_genome/v2/hs1_"

bin_size_dict={"1kb":1000,"10kb":10000,"100kb":100000}


chrom_sizes_df=pd.read_csv(chrom_sizes_file, sep='\t', header=None, usecols=[0,1,2], names=["chromosome","start_coord","stop_coord"])

bin_size_df=pd.DataFrame.from_dict(bin_size_dict, orient="index", columns=["bin_size"])




for i in chrom_sizes_df['chromosome']:
    my_chr_stop=chrom_sizes_df.loc[chrom_sizes_df["chromosome"]==i,"stop_coord"].values[0]
    bin_size_df[i]=my_chr_stop/bin_size_df["bin_size"] + 1


for k in bin_size_dict.keys():
    
    a=bin_size_dict[k]
    

    my_anno_df=pd.DataFrame(columns=["chromosome","start","stop"])
    for c in chrom_sizes_df['chromosome']:
        chr_bins = bin_size_df.loc[bin_size_df["bin_size"]==a, c].values[0]
        #chr_bin_size = bin_size_df.loc[bin_size_df["cen_annotation_name"]==a, 'size'].values[0]
        start = np.arange(0, chr_bins * a, a)[:-1] + 1
        stop = np.arange(a, chr_bins * a + 1, a)
        stop[-1] = chrom_sizes_df.loc[chrom_sizes_df["chromosome"]==c, 'stop_coord'].values[0]

        if len(start) != len(stop):
            print("no matchy", a, c)
            continue
        chr_bins_df = pd.DataFrame({'chromosome': c, 'start': start, 'stop': stop})
        chr_bins_df['start'] = chr_bins_df['start'].astype(str).str.split('.').str[0]
        chr_bins_df['stop'] = chr_bins_df['stop'].astype(str).str.split('.').str[0]
    
        #create a bin name
        chr_bins_df.reset_index(inplace=True)
        chr_bins_df['index']=chr_bins_df['index']+1
        chr_bins_df['bin_name']=chr_bins_df['chromosome'].astype(str)+"_bin"+chr_bins_df['index'].astype(str)
        chr_bins_df.drop(columns=['index'], inplace=True)
        my_anno_df = pd.concat([my_anno_df, chr_bins_df], ignore_index=True)

    
    file_name = output_filename_prefix + k + "_genome_bins.bed"

    my_anno_df.to_csv(file_name, sep="\t", index=False, header=False)

    #print(my_anno_df)
    