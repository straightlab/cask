import os
import json
from cask import kdb_parser as kdb_parser
from subprocess import check_output
import pickle
import sys
from multiprocessing import cpu_count, Pool

def wc(filename):
    return int(check_output(["wc", "-l", filename]).split()[0])

nfiles=len(sys.argv)-2

#os.chdir(sys.argv[1])
#rootdir=os.getcwd()
rootdir=os.path.dirname(sys.argv[1])
typelist=os.path.join(rootdir,'repeats.withClass.txt')

typelist=sys.argv[1]
klut=sys.argv[2]
#klut=os.path.join(rootdir, klut_relpath)

dbpath=sys.argv[3]
dbname=sys.argv[4]

fqlist_file=sys.argv[5]

nprocess=cpu_count()
if len(sys.argv)==7:
    nprocess=int(sys.argv[6])

print("Parsing %s"%klut)
print("Using types in %s"%typelist)
nmax=wc(klut)+1


p=kdb_parser.Annotated_KDB(klut, nmax, typelist, with_class=True)

filein=open(fqlist_file,"r")
fqroots=[]
for line in filein:
    f0=line.strip()
    fdir=os.path.dirname(f0)
    fname=os.path.basename(f0)
    fqroot=os.path.join(fdir,dbpath,fname+".ANNOTATED."+dbname)
    fqroots+=[fqroot]

print("Ready to heavy lift! Files to be processed: ")
print(fqroots)

agg=kdb_parser.AmbivalenceGroupGenerator(p)

for f in fqroots:
#     heavy_lift(f)
# def heavy_lift(f):
    fq=os.path.join(f+".fastq.gz")
    print("Computing ag for : %s"%fq)
    outname=f+".ag.txt"
    agg.compute_ambivalence_groups(outname,0,fq)
    print('AG computed : %s'%outname)
    print("Preparint the LUT to save")
    ix_to_chr_LUD={v:k for k, v in agg.kdb.typeid_LUD.items()}
    ag_to_chr_LUD={v:[ix_to_chr_LUD[kk] for kk in k] for k, v in agg.ag_dict.items()}

    ag_to_chr_LUD[-1]=["-1"]
    ag_to_chr_LUD[0]=["0"]
    for k, v in ix_to_chr_LUD.items():
        ag_to_chr_LUD[k]=[v]

    ix_to_chr_LUD_class={v:k for k, v in agg.kdb.classid_LUD.items()}

    ag_to_chr_LUD_class={v:[ix_to_chr_LUD_class[kk] for kk in k] for k, v in agg.aG_dict.items()}
    ag_to_chr_LUD_class[-1]=["-1"]
    ag_to_chr_LUD_class[0]=["0"]
    for k, v in ix_to_chr_LUD_class.items():
        ag_to_chr_LUD_class[k]=[v]

    mydata2pickle={'ag_to_chr_LUD':ag_to_chr_LUD,'ag_to_chr_LUD_class':ag_to_chr_LUD_class}
    parserout=f+".aglut.pickle"
    with open(parserout, 'wb') as f1:
        pickle.dump(mydata2pickle, f1)
    print('LUT saved : %s'%parserout)

# pool = Pool(processes=nprocess)
# res = pool.map(heavy_lift, fqroots)
# pool.close()

print("ALL DONE")
