import itertools
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--path")
parser.add_argument("--input")
parser.add_argument("--outdir")
args = parser.parse_args()

gene_info={}
with open(args.input,'r') as IN:
    for line in IN:
        items=line.strip().split()
        gene_info[items[0].split('.')[0]]=items[1:]

c=0
genes={}
for file in os.listdir(args.path):
    gene=file.split('.')[0]
    c+=1
    if c%400 not in genes:
        genes[c%400]=[]
    genes[c%400].append(gene+'\t'+'\t'.join(gene_info[gene]))
for batch in genes:
    with open(args.outdir+'genes_batch'+str(batch)+'.txt','w') as OUT:
        OUT.write('\n'.join(genes[batch]))
