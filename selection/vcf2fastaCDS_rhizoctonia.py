import itertools
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import gzip
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument("--listvcf")
parser.add_argument("--gff")
parser.add_argument("--out")
parser.add_argument("--path")
parser.add_argument("--ref")
args = parser.parse_args()

IUPAC={'A|G':'R',
        'G|A':'R',
        'C|T':'Y',
        'T|C':'Y',
       'G|C':'S',
    'C|G':'S',
       'A|T':'W',
    'T|A':'W',
       'G|T':'K',
    'T|G':'K',
       'A|C':'M',
    'C|A':'M',
    'N|N':'N',
    'A|A':'A',
    'C|C':'C',
    'G|G':'G',
    'T|T':'T'}

with open(args.ref,'r') as IN:
    chromosomes = SeqIO.to_dict(SeqIO.parse(IN, "fasta"))

genes={}
with open(args.listvcf,'r') as IN:
    for line in IN:
        items=line.strip().split()
        gene=items[0]
        dir_=items[1]
        chr_=items[2].split(":")[0]
        if gene not in genes:
            genes[gene]={'positions':{},'chr':chr_,'dir':dir_}
        for coords in items[2].split(','):
            print(coords)
            position1=int(coords.split(':')[1].split('-')[0])
            position2=int(coords.split(':')[1].split('-')[1])
            genes[gene]['positions'][position1]=[position1,position2]
#read and parse REF NLR sequences and VCF
for gene in genes:
    genotypes={}
    sequences={}
    index2isolate={}
    print(gene,genes[gene])
    for file in os.listdir(args.path):
        if file.startswith(gene):
            myfile=file
    with gzip.open(args.path+myfile, "rt") as VCF:
        for line in VCF:
            items=line.strip().split()
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    for i in range(9,len(items)):
                        index2isolate[i]=items[i]
                        sequences[items[i]]=[]
            else:
                pos=int(items[1])
                genotypes[pos]={}
                alleles={'0':items[3],'.':'N'}
                if ',' in items[4]:
                    d=0
                    for allele in items[4].split(','):
                        d+=1
                        alleles[str(d)]=allele
                else:
                    alleles['1']=items[4]
#                print(alleles)
                for i in range(9,len(items)):
                    nucl1=alleles[items[i].split(':')[0].replace('|','')[0]]
                    nucl2=alleles[items[i].split(':')[0].replace('|','')[1]]
                    if len(nucl1)!=1:
                        if len(nucl1)==len(nucl2):
                            seq=[]
                            for j in range(len(nucl1)):
                                seq.append(IUPAC[nucl1[j]+'|'+nucl2[j]])
                            genotypes[pos][index2isolate[i]]=''.join(seq)
                        else:
                            genotypes[pos][index2isolate[i]]=len(nucl1)*'N'
                    else:
                        if len(nucl2)!=1:
                            genotypes[pos][index2isolate[i]]='N'
                        else:
                            genotypes[pos][index2isolate[i]]=IUPAC[nucl1+'|'+nucl2]
    chromosome=chromosomes[genes[gene]['chr']]
    for position in sorted(list(genes[gene]['positions'].keys())):
        interval=genes[gene]['positions'][position]
        for i in range(interval[0],interval[1]+1):
            if i in genotypes:
                for ind in genotypes[i]:
                    sequences[ind].append(genotypes[i][ind])
            else:
                for ind in genotypes[i]:
                    sequences[ind].append('N')
    records=[]
    for ind in sequences:
#        print(gene,dirs[gene])
        if genes[gene]['dir']=='+':
            records.append(SeqRecord(Seq(''.join(sequences[ind])), id = ind, description = ''))
        else:
            records.append(SeqRecord(Seq(''.join(sequences[ind])).reverse_complement(), id = ind, description = ''))
#    print(records)
    SeqIO.write(records,args.out+gene+'.fasta','fasta')
