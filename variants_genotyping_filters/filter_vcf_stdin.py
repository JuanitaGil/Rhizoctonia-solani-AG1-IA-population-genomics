import fileinput

for line in fileinput.input():
    line=line.strip()
    if line.startswith('#'):
        print(line)
    else:
        items=line.split('\t')
        if items[5]=='.':
            print(line)
        else:
            QUAL=float(items[5])
            MQ=float(items[7].split(';')[-1].split('=')[1])
            if QUAL>=5 and MQ>=10:
                print(line)
            else:
                print('\t'.join(items[0:8])+'\tGT\t'+'\t'.join(145*['./.']))
