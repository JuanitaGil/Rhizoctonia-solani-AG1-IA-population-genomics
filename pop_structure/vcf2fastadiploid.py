import numpy as np
import itertools
import os
import argparse
import gzip
parser = argparse.ArgumentParser()
parser.add_argument("--vcf")
parser.add_argument("--out")
args = parser.parse_args()
IUPAC={'A/G':'R',
        'G/A':'R',
        'C/T':'Y',
        'T/C':'Y',
       'G/C':'S',
    'C/G':'S',
       'A/T':'W',
    'T/A':'W',
       'G/T':'K',
    'T/G':'K',
       'A/C':'M',
    'C/A':'M',
    'N/N':'N',
    'A/A':'A',
    'C/C':'C',
    'G/G':'G',
    'T/T':'T'}

############
index2isolate={}
c=0
NLR_nucls={}
with gzip.open(args.vcf, "rt") as VCF,open(args.out,'w') as OUT:
    for line in VCF:
        items=line.strip().split()
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                for i in range(9,len(items)):
                    index2isolate[i]=items[i]
                    NLR_nucls[items[i]]=[]
        else:
            pos=items[1]
            alleles={'0':items[3],'.':'N'}
            if ',' in items[4]:
                d=0
                for allele in items[4].split(','):
                    d+=1
                    alleles[str(d)]=allele
            else:
                alleles['1']=items[4]
            c+=1
            if c%10000==0:
                print(c)
            for i in range(9,len(items)):
                allele1=alleles[items[i].split(':')[0].split('/')[0]]
                allele2=alleles[items[i].split(':')[0].split('/')[1]]
                NLR_nucls[index2isolate[i]].append(IUPAC[allele1+'/'+allele2])
    for accession in NLR_nucls:
        OUT.write('>'+accession+'\n'+''.join(NLR_nucls[accession])+'\n')
