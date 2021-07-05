import csv
import gzip
import json
import pickle


class Annotated_KDB: 

    def __init__(self, klut, nmax, typelist, with_class=False): 
        
        if nmax<0: #preloaded
            print("Loading precomputed types LUD and kmers LUT")
            with open(klut) as f:
                json_data=f.read()
                x=json.loads(json_data)
                self.kmer_LUT=x[0]
                self.kmer_count_LUT=x[1]
                self.typeid_LUD=x[2]
                if with_class and (len(x)>3):
                    self.classid_LUD=x[3]
                    self.class_LUT=x[4]
                else:
                    if with_class:
                        print("Watch out, no class level LUT found in json file")
                    self.classid_LUD=None
                    self.class_LUT=None
        else:
            print("Computing types LUD...")
            typeid_LUD, classid_LUD, class_LUT = Annotated_KDB.load_typelist(typelist, with_class=with_class)
            self.typeid_LUD=typeid_LUD #{type: ID}
            self.classid_LUD=classid_LUD
            self.class_LUT=class_LUT #[1,10,3,...] the class to which each type belongs to project the types onto class id (can be used to regroup different types together, higher level annotation)
            print("Computing kmers LUT...")
            self.load_kmerlut(klut, nmax)
        print("Done.")

    @staticmethod
    def load_typelist(file, with_class=False):
        with open(file,'r') as f:
            myreader=csv.reader(f, delimiter='\t')
            # myreader.__next__() #get rid of header
            typeid_LUD={row[0]:(ix+1) for ix, row in enumerate(myreader)} #0 is unknown!
        if with_class:
            n=len(typeid_LUD)
            class_LUT=[0]*(n+1)
            classid_LUD={} #
            nmax=0
            with open(file,'r') as f:
                myreader=csv.reader(f, delimiter='\t')
                for ix, row in enumerate(myreader):
                    # print(ix)
                    # print(row)
                    x=classid_LUD.get(row[1],0)
                    if x==0: #new class
                        nmax+=1
                        classid_LUD[row[1]]=nmax
                        x=nmax
                    class_LUT[ix+1]=x
            return typeid_LUD, classid_LUD, class_LUT
        else:
            return typeid_LUD, None, None

    def load_kmerlut(self, kmers_annotations, nmax):
        self.kmer_counts=[0]*nmax
        self.kmer_LUT=[[] for _ in range(nmax)] #[[kid0_typeid0,kid0_typeid1,...],[kid1_typeid0,kid1_typeid1,...]]
        self.kmer_count_LUT=[[] for _ in range(nmax)] #[[kid0_COUNTtypeid0,kid0_COUNTtypeid1,...],[kid1_COUNTtypeid0,kid1_COUNTtypeid1,...]]
        nmax=0
#         1       38      AluJb,5,AluJo,1,AluJr,1,AluSc,1,AluSc5,1,AluSg,2,AluSg7,2,AluSp,1,AluSq,2,AluSq2,2,AluSx,2,AluSx1,3,AluSx4,1,AluSz,2,AluY,3,AluYc,1,AluYj4,1,(A)n,1,(T)n,5,(TT)n,1
# 2       28      AluJb,4,
        with open(kmers_annotations,'r') as f:
            myreader=csv.reader(f, delimiter='\t')
            
            for _, data in enumerate(myreader):
                ix_kmer=int(data[0])

                self.kmer_counts[ix_kmer]=data[1]
                klist=data[2].split(",")
                nk=int(len(klist)/2)
                for i in range(nk):
                    kix=self.typeid_LUD.get(klist[2*i],0)
                    if kix>0:
                        self.kmer_LUT[ix_kmer]+=[kix]
                        self.kmer_count_LUT[ix_kmer]+=[float(klist[2*i+1])]
        
    def get_ktypes_union(self,kid_list):
        s=set([])
        for kid in kid_list:
            if kid>0:
                s.update(self.kmer_LUT[kid])
        sl=list(s)
        sl.sort()
        ktypes=tuple(sl)
        if not self.class_LUT is None:
            s2=set([])
            for kt in ktypes:
                s2.add(self.class_LUT[kt])
            sl2=list(s2)
            sl2.sort()
            kclass=tuple(sl2)
        else:
            kclass=()
        return ktypes, kclass

    def get_ktypes_intersect_alllevels(self,kid_list):
        
        s=set(self.kmer_LUT[kid_list[0]])
        for i in range(1,len(kid_list)):
            s.intersection_update(self.kmer_LUT[kid_list[i]])
        sl=list(s)
        sl.sort()
        ktypes=tuple(sl)

        if (not self.class_LUT is None):
            # print(ktypes[0])
            
            s2=set([self.class_LUT[kt] for kt in self.kmer_LUT[kid_list[0]]])
            for i in range(1,len(kid_list)):
                s2.intersection_update([self.class_LUT[kt] for kt in self.kmer_LUT[kid_list[i]]])
            sl2=list(s2)
            sl2.sort()
            kclass=tuple(sl2)
        else:
            kclass=()
        return ktypes, kclass

    def get_ktypes_intersect(self,kid_list):
        s=set(self.kmer_LUT[kid_list[0]])
        for i in range(1,len(kid_list)):
            s.intersection_update(self.kmer_LUT[kid_list[i]])
        sl=list(s)
        sl.sort()
        ktypes=tuple(sl)

        if (len(sl)>0) and (not self.class_LUT is None):
            # print(ktypes[0])
            
            s2=set([])
            for kt in ktypes:
                s2.add(self.class_LUT[kt])
            sl2=list(s2)
            sl2.sort()
            kclass=tuple(sl2)
        else:
            kclass=()
        return ktypes, kclass
            
    def save(self, kdbfile):
        with open(kdbfile, 'w') as f:
            if not self.class_LUT is None:
                json.dump([self.kmer_LUT, self.kmer_count_LUT, self.typeid_LUD, self.classid_LUD, self.class_LUT], f, indent=4)
            else:
                json.dump([self.kmer_LUT, self.kmer_count_LUT, self.typeid_LUD], f, indent=4)


class AmbivalenceGroupGenerator:
    def __init__(self, kdb, agfile=None): #kdb is a kmer database
        self.kdb=kdb
        # if agfile is None:
        self.ag_dict={}

        self.nag=-1 #neg are incompatible ones
        self.naG=-1
        self.aG_dict={}
        # else:


    def _compute_ambivalence_group(self, kid_list):
        if len(kid_list)==0:
            return 0
        ktypes, kclass=self.kdb.get_ktypes_intersect_alllevels(kid_list)
        
        if len(ktypes)==0: #null intersection, means that there are conflicting assignments
            ag=-1
        elif len(ktypes)==1: #non ambivalent at the kmer level
            ag=ktypes[0]
        else:
            ag=self.ag_dict.get(ktypes,0)
            if ag==0:
                self.nag-=1
                ag=self.nag
                self.ag_dict[ktypes]=ag
        
        aG=0
        if len(kclass)==0:
            aG=-1
        elif len(kclass)==1: #non ambivalent at the class level
            aG=kclass[0]
        else:
            aG=self.aG_dict.get(kclass,0)
            if aG==0:
                self.naG-=1
                aG=self.naG
                self.aG_dict[kclass]=aG
    
        return ag, aG

                    
    def compute_ambivalence_group_streams(self, sync_streams, out_stream, nmax=0): #streams are syncronized
        n=0
        while ((nmax==0) or (n<nmax)):
            try:
                record_id, kdata=sync_streams.__next__()
            except StopIteration:
                break
            n+=1
            out_stream.write(record_id)
            for kid_list in kdata:
                ag, aG=self._compute_ambivalence_group(kid_list) 
                out_stream.write("\t")
                out_stream.write(str(ag))
                out_stream.write("\t")
                out_stream.write(str(aG))
                out_stream.write("\t")
                out_stream.write(str(len(kid_list)))
            out_stream.write("\n")
        return n

    # def save_ag(filename):

    
    def compute_ambivalence_groups(self, outfilename, nmax, *args):
        sync_streams=KFQReader(*args)
        out_stream=open(outfilename, 'w')
        self.compute_ambivalence_group_streams(sync_streams,out_stream, nmax=nmax)
        out_stream.close()
        sync_streams.close()

    def save(self, outjson):
        with open(outjson, 'w') as f:
            json.dump([{",".join([str(kk) for kk in k]):v for k,v in self.ag_dict.items()}, {",".join([str(kk) for kk in k]):v for k,v in self.aG_dict.items()}, self.nag, self.naG], f, indent=4)
            
    def save_pickle(self, out):
        with open(out, 'wb') as handle:
            pickle.dump({'ag_dict':self.ag_dict,'aG_dict':self.aG_dict,'nag':self.nag, 'naG':self.naG}, handle, protocol=pickle.HIGHEST_PROTOCOL)

    def load(self, fname):
        with open(fname, 'r') as handle:
            ags=json.load(handle)
        self.ag_dict={tuple(int(x) for x in k.split(',')):v for k, v in ags[0].items()}
        self.aG_dict={tuple(int(x) for x in k.split(',')):v for k, v in ags[1].items()}
        self.nag=ags[2]
        self.naG=ags[3]


class KFQReader():
    def __init__(self, *args):
        self.readers=[gzip.open(arg, 'rt') for arg in args]
        self.N=len(args)

    def __iter__(self):
        return self

    def __next__(self):
        L=self.readers[0].readline()
        self.readers[0].readline()
        self.readers[0].readline()
        self.readers[0].readline()
        
        if len(L)<2:
            # self.close()
            raise StopIteration

        L_split=L.split("\t")
        record_id=L_split[0][1:]
        h=[int((x.split("=")[0]).split(",")[0]) for x in L_split[1:]]
        kdata=[h]


        return record_id, kdata

    def close(self):
        for reader in self.readers:
            reader.close()
   