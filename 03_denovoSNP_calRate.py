import argparse

###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--index',metavar='File',dest='index',help='Input index file',type=open,required=True)
parser.add_argument('-m','--matrix',metavar='File',dest='matrix',help='Input matrix file',type=open,required=True)
parser.add_argument('-t','--threshold',metavar='Float',dest='threshold',help='mutation rate below the threshold will be selected',type=float,default=1.0)
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

def callFreq(c,threshold):
    lenc = len(c)
    count0 = c.count(0)
    count1 = c.count(1)
    count2 = c.count(2)
    count3 = c.count(3)
    realLen = (lenc-count0)*2
    refFreq = 2*count1+count2
    altFreq = 2*count3+count2
    cc = min(refFreq,altFreq)
    if cc <= threshold:
        ccc = cc
    else:
        ccc = 0
    return ccc,realLen

#########################################

def calRate(ifile, mfile, wfile, threshold):
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
                    cc, realLen = callFreq(c,threshold)
                    freqbbin1Dict[i]['cc'].append(cc)
                    freqbbin1Dict[i]['realLen'].append(realLen)
                else:
                    freqbbin1Dict[i]['cc'] = [0]
                    freqbbin1Dict[i]['realLen'] = [1]

        for j in bbin3.keys():
            if j[0]*1000000 < int(line[0]) < j[1]*1000000:
                if len(bbin3[j]) > 1:
                    c = [int(line[i]) for i in bbin3[j]]
                    cc,realLen = callFreq(c,threshold)
                    freqbbin3Dict[j]['cc'].append(cc)
                    freqbbin3Dict[j]['realLen'].append(realLen)
                else:
                    freqbbin3Dict[j]['cc'] = [0]
                    freqbbin3Dict[j]['realLen'] = [1]

####################################################################

    denovoMutationRate1 = {}
    for key in freqbbin1Dict.keys():
        denovoMutationRate1[key] = []
        #ssall = sum(freqbbin1Dict[key]['realLen'])
        ssdenovo = sum(freqbbin1Dict[key]['cc'])
        denovoMutationRate1[key] = ssdenovo/(key[1]-key[0]+1)

    denovoMutationRate3 = {}
    for key in freqbbin3Dict.keys():
        denovoMutationRate3[key] = []
        #ssall = sum(freqbbin3Dict[key]['realLen'])
        ssdenovo = sum(freqbbin3Dict[key]['cc'])
        denovoMutationRate3[key] = ssdenovo/(key[1]-key[0]+1)

####################################################################

    for i in denovoMutationRate1.keys():
        wfile.write('{}\t{}\t{}\t{}\n'.format(i,(i[1]-i[0]+1),'%.5f'%(denovoMutationRate1[i]),'%.5f'%(denovoMutationRate3[i])))
    wfile.close()

calRate(ifile=args.index,mfile=args.matrix,wfile=args.output,threshold=args.threshold)
