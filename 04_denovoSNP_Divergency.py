import argparse
from itertools import combinations

###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--index',metavar='File',dest='index',help='Input index file',type=open,required=True)
parser.add_argument('-m','--matrix',metavar='File',dest='matrix',help='Input matrix file',type=open,required=True)
parser.add_argument('-o','--Output',metavar='File',dest='output',help='Output file',type=argparse.FileType('w'),required=True)
args = parser.parse_args()
###### arguments ######



##############################################################################
# This script will call denovo mutation rate for each haplotype of a whole block
# 2018-06-10
# This script is improved from the two aspects:
# 1. denovo SNP is divided by BLOCK LENGTH;
# 2. callFreq has the threshold now
##############################################################################

def bigRemove(n,l):
    while n in l:
        l.remove(n)
    return l

def pair2divergency(pair):
    if pair in [(1,1), (2,2), (3,3)]:
        return 0
    elif pair in [(1,2), (2,1), (2,3), (3,2)]:
        return 0.5
    elif pair in [(1,3), (3,1)]:
        return 1


def make_pair(c):
    c = bigRemove(0,c)
    cpair = list(combinations(c, 2))
    pair = []
    for i in cpair:
        pair.append(pair2divergency(i))
    return sum(pair)/len(cpair)

#########################################

def calRate(ifile, mfile, wfile):
    bbin1 = {}
    bbin3 = {}
    for line in ifile:
        line = line.strip().split()

        # if Haplotype I is empty
        if line[2] == '0':
            bbin1[(int(line[0]), int(line[1]))] = [0]
            bbin3[(int(line[0]), int(line[1]))] = [int(i) for i in line[3].split(':')]
        # if Haplotype III is empty
        elif line[3] == '0':
            bbin1[(int(line[0]), int(line[1]))] = [int(i) for i in line[2].split(':')]
            bbin3[(int(line[0]), int(line[1]))] = [0]
        # if both H I and H III are not empty
        else:
            bbin1[(int(line[0]), int(line[1]))] = [int(i) for i in line[2].split(':')]
            bbin3[(int(line[0]), int(line[1]))] = [int(i) for i in line[3].split(':')]

####################################################################
    # INItialize a new index dictionary
    freqbbin1Dict = {}
    for i in bbin1.keys():
        freqbbin1Dict[i] = {}
        freqbbin1Dict[i]['cc'] = []
        freqbbin1Dict[i]['realLen'] = []

    freqbbin3Dict = {}
    for i in bbin3.keys():
        freqbbin3Dict[i] = {}
        freqbbin3Dict[i]['cc'] = []
        freqbbin3Dict[i]['realLen'] = []

    # LOOKING THROUGH MATRIX FILE
    for line in mfile:
        line = line.strip().split()
        for i in bbin1.keys():
            if i[0]*1000000 < int(line[0]) < i[1]*1000000:
                if len(bbin1[i]) > 1:
                    c = [int(line[i]) for i in bbin1[i]]
                    cc = make_pair(c)
                    if cc/((i[1]-i[0]+1)*1000000) > 0.00067573:
                        freqbbin1Dict[i]['cc'].append(int(line[0]))
                else:
                    freqbbin1Dict[i]['cc'] = [0]
                    freqbbin1Dict[i]['realLen'] = [1]

        for j in bbin3.keys():
            if j[0]*1000000 < int(line[0]) < j[1]*1000000:
                if len(bbin3[j]) > 1:
                    c = [int(line[i]) for i in bbin3[j]]
                    cc = make_pair(c)
                    if cc/((j[1]-j[0]+1)*1000000) > 0.00067573:
                        freqbbin3Dict[j]['cc'].append(int(line[0]))
                else:
                    freqbbin3Dict[j]['cc'] = [0]
                    freqbbin3Dict[j]['realLen'] = [1]

####################################################################

    for i in bbin1.keys():
        wfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
                    .format(i[0],i[1],(i[1]-i[0]+1)
                            ,len(bbin1[i]),'%.8f'%(len(bbin1[i])/((i[1]-i[0]+1)*1000000))
                            ,len(bbin3[i]),'%.8f'%(len(bbin3[i])/((i[1]-i[0]+1)*1000000))))
    wfile.close()

calRate(ifile=args.index,mfile=args.matrix,wfile=args.output)
