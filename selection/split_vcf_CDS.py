import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--input")
parser.add_argument("--scripts")
parser.add_argument("--vcf")
parser.add_argument("--outdir")
parser.add_argument("--batchsize")
args = parser.parse_args()

genes={}
with open(args.input,'r') as GENEs:
    for line in GENEs:
        items=line.strip().split()
#        print(items)
        genes[items[0]]=items[-1]

c=0
commands={}
for gene in genes:
    c+=1
    win=c%int(args.batchsize)
#    print(win)
    if win not in commands:
        commands[win]=[]
    commands[win].append('bcftools view -r '+genes[gene].replace('Chr','')+' '+args.vcf+' | bgzip > '+args.outdir+gene+'.vcf.gz;\n')
#print(commands)
for win in commands:
    with open(args.scripts+'batch'+str(win)+'.sh','w') as SH:
#        SH.write('#!/bin/bash\n#SBATCH -J fas'+str(win)+'\n')
#        SH.write('#SBATCH -o fas'+str(win)+'.out\n')
#        SH.write('#SBATCH -e fas'+str(win)+'.err\n')
#        SH.write('module load htslib/1.9\nmodule load bcftools/1.9\n')
        SH.write(''.join(commands[win]))
#    os.system('sbatch ./sh_scripts/batch'+str(win)+'.sh')
