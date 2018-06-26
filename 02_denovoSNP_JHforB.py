import argparse
import numpy as np
from collections import defaultdict

###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input',metavar='File',dest='input',help='Input file',type=open,required=True)
parser.add_argument('-pf','--position_file',metavar='File',dest='position',help='Input position file',type=open,required=True)
parser.add_argument('-s','--size',metavar='Int',dest='size',help='size',type=int,default=10)
parser.add_argument('-o','--Output',metavar='File',dest='output',help='Output file',type=argparse.FileType('w'),required=True)
args = parser.parse_args()
###### arguments ######

# This script will judge the haplotype for each block

##################################################

def loc_block_bins(posfile,size):
    binList = []
    for line in posfile:
        line = line.strip().split()
        binList.append(int(line[0]))

    bbin = []
    for i in range(0,len(binList)-1):
        if abs(binList[i] - binList[i+1]) >= size:
            bbin.append((binList[i],binList[i+1]))

    # CHECK CONSISTENCY
    new = []; abort = []
    for j in range(0,len(bbin)-1):
        if bbin[j][1] == bbin[j+1][0]:
            k1 = bbin[j][0]; k2 = bbin[j+1][1]
            abort.append(bbin[j]); abort.append(bbin[j+1])
            new.append((k1,k2))

    for k in abort:
        bbin.remove(k)

    for l in new:
        bbin.append(l)

    return bbin

######################################################

def drop_into_num(mean_ab):
    norm_lv_block = []
    col_num = mean_ab.shape[1]
    for i in range(0,col_num):
        if 0.0 <= mean_ab[0,i] < 0.5:
            norm_lv_block.append(0)
        elif 0.5 <= mean_ab[0,i] < 1.5:
            norm_lv_block.append(1)
        elif 1.5 <= mean_ab[0,i] < 2.5:
            norm_lv_block.append(2)
        elif 2.5 <= mean_ab[0,i]:
            norm_lv_block.append(3)
    return norm_lv_block


######################################################

def find_deletion(source,deletion):
    deletion_index=[]
    s_index = 0;e_index = len(source)
    while(s_index < e_index):
        try:
            temp = source.index(deletion,s_index,e_index)
            deletion_index.append(temp+1)
            s_index = temp + 1
        except ValueError:
            break
    # 1 BASE
    return deletion_index

######################################################

def list2str(a):
    cc = [str(i)+':' for i in a]
    d = ''
    for i in cc:
        d += i
    ddd = d[0:-1]
    return ddd

######################################################



def denovoSNP_jhforb(ifile,posfile,size,wfile):
    """
    197523  2       2       2       2       1       2       2       1       2       1       1       2       1       1       1       1
    198048  1       1       2       2       1       1       2       2       1       1       2       1       1       1       1       2
    198163  1       1       2       2       1       1       1       2       1       1       2       1       1       1       1       2
    198249  1       1       2       2       1       2       2       2       2       1       2       2       1       1       1       2
    263400  1       1       1       1       1       1       1       1       0       1       1       1       1       1       1       1
    """

    bbin = loc_block_bins(posfile,size)
    # Initialization check block
    ab = defaultdict(list)
    for i in bbin:
        ab[i[0]] = []

    for line in ifile:
        line = line.strip().split()
        pos = int(line[0])
        #######LANDRACE AND VARIETY#######
        lv1 = line[1:46]; lv2 = line[47:49]; lv3 = line[60:65]; lv4 = line[68:73]; lv5 = line[83::]
        lv = lv1+lv2+lv3+lv4+lv5
        numlv = np.array([int(i) for i in lv])

        #######LANDRACE AND VARIETY#######

        for key in ab.keys():
            if key*1000000 < int(pos) < (key+size)*1000000:
                ab[key].append(numlv)

    mean_ab = {}
    haploBlock = {}
    for k in ab.keys():
        catall = np.matrix(ab[k])
        mean_ab[k]=drop_into_num(np.mean(catall,0))
    for l in mean_ab.keys():
        haploBlock[l] = []
        ddd1 = list2str(find_deletion(mean_ab[l],1)); ddd3 = list2str(find_deletion(mean_ab[l],3))
        haploBlock[l].append(ddd1); haploBlock[l].append(ddd3)
    print(haploBlock)

    ### Replace key with full region ###
    for m in haploBlock.keys():
        for n in bbin:
            if int(m) == int(n[0]):
                wfile.write('{}\t{}\t{}\t{}\n'.format(n[0],n[1],haploBlock[m][0],haploBlock[m][1]))


denovoSNP_jhforb(ifile=args.input,posfile=args.position,size=args.size,wfile=args.output)
