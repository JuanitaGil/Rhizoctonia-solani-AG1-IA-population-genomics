import egglib
import sys
import os
import argparse
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument("--path")
parser.add_argument("--min_size")
parser.add_argument("--min_length")
parser.add_argument("--max_missing")
parser.add_argument("--pathout")
args = parser.parse_args()
c=0
with open(args.pathout,'w') as OUT:
    OUT.write('Name\tlseff\tnseff\tS\tHe\tthetaW\tPi\tD\tnumNS\tnumS\tPiN\tPiS\tPiN/PiS\n')
    c+=1
    output=[]
    aln = egglib.io.from_fasta(args.path, alphabet=egglib.alphabets.DNA, labels=True)
    #struct = egglib.get_structure(aln,lvl_clust=0,lvl_pop=1,lvl_indiv=2)
    cs = egglib.stats.ComputeStats()
    cs.add_stats('lseff','nseff','S','He','thetaW','Pi','D')
    stat_order=['lseff','nseff','S','He','thetaW','Pi','D']
    stats = cs.process_align(aln)#,max_missing=30)
    output.append(args.path.split('/')[-1].split('.')[0])
    output.append(stats['lseff'])
#    print(stats)
    if stats['lseff']>int(args.min_length) and stats['nseff']>int(args.min_size):
        print(c,args.path,'processed')
        output.append(stats['nseff'])
        output.append(stats['S'])
        output.append(stats['He'])
        output.append(stats['thetaW']/stats['lseff'])
        output.append(stats['Pi']/stats['lseff'])
        output.append(stats['D'])
        aln.to_codons()
        cd = egglib.stats.CodingDiversity(aln,max_missing=float(args.max_missing))
        sitesNS=cd.sites_NS
        sitesS=cd.sites_S
        numS=cd.num_sites_S
        numNS=cd.num_sites_NS
        codoneff=cd.num_codons_eff
        cs = egglib.stats.ComputeStats()
        cs.add_stats('nseff','S','Pi','lseff')
        statsS = cs.process_sites(sitesS)#,max_missing=float(args.max_missing))
        statsNS = cs.process_sites(sitesNS)#,max_missing=float(args.max_missing))
        PiS=statsS['Pi']
        PiNS=statsNS['Pi']
        print(statsNS['S'],statsS['S'])
        if numNS!=None and numNS!=0 and PiNS!=None:
            pins=PiNS/numNS
        else:
            pins=None
        if numS!=None and numS!=0 and PiS!=None:
            pis=PiS/numS
        else:
            pis=None
        if pins!=0 and pins!=None and pis!=None:
            pinspis=pins/pis
        else:
            pinspis=None
        output.append(numNS)
        output.append(numS)
        output.append(pins)
        output.append(pis)
        output.append(pinspis)
    else:
        print(c,args.path,'not enough sites left')
        output.append('None')
        output.append('None')
        output.append('None')
        output.append('None')
        output.append('None')
        output.append('None')
        output.append('None')
        output.append('None')
        output.append('None')
        output.append('None')
        output.append('None')
    OUT.write('\t'.join([str(x) for x in output])+'\n')
