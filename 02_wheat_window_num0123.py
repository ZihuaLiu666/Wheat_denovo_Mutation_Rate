import argparse
from collections import defaultdict
import math
import numpy as np


###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input',metavar='File',dest='input',help='Input file',type=open,required=True)
parser.add_argument('-chr','--chromosome',metavar='Str',dest='chr',help='chromosome',type=str,default='chr1A')
parser.add_argument('-l','--length',metavar='File',dest='length',help='length file',type=open,required=True)
parser.add_argument('-w','--window',metavar='Int',dest='window',help='window',type=int,default=5000)
parser.add_argument('-s','--step',metavar='Int',dest='step',help='step',type=int,default=5000)
parser.add_argument('-o','--Output',metavar='File',dest='output',help='Output file',type=argparse.FileType('w'),required=True)
args = parser.parse_args()
###### arguments ######

def look_for_chr_length(length_file,chr):
    for line in length_file:
        line = line.strip().split()
        if line[0] == chr:
            ll = int(line[1])
            return ll


def drop_into_num(lv_block):
    norm_lv_block = []
    col_num = lv_block.shape[1]
    for i in range(0,col_num):
        if 0.0 <= lv_block[:,i] < 0.5:
            norm_lv_block.append(0)
        elif 0.5 <= lv_block[:,i] < 1.5:
            norm_lv_block.append(1)
        elif 1.5 <= lv_block[:,i] < 2.5:
            norm_lv_block.append(2)
        elif 2.5 <= lv_block[:,i]:
            norm_lv_block.append(3)
    return norm_lv_block

def make_bin_and_write(length_file,chr,ifile,window,step,wfile):

    ll = look_for_chr_length(length_file,chr)

    """
    77488   0       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1
    77617   1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1
    77641   1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1
    77662   1       1       3       1       1       1       1       1       3       3       3       1       3       3       1       1       3       1       3       1       1       1
    77663   1       1       3       1       1       1       1       1       3       3       3       1       3       3       1       1       3       1       3       1       1       1
    77856   1       1       3       1       1       1       1       1       3       3       3       1       3       3       1       1       3       1       3       1       1       1
    77875   1       1       3       1       1       0       1       1       3       3       3       1       3       3       1       1       3       1       3       1       1       1
    """

    ## COLLECTING DATA
    ab = {}
    for lline in ifile:
        lline = lline.strip().split()
        ab[int(lline[0])] = np.array([int(i) for i in lline[1::]])


    aab = defaultdict(list)
    i = 1
    total_bin = int(math.ceil(((ll-window)/step)+1))

    while i <= total_bin:

        start_bin = 1+(i-1)*step
        end_bin = window + (i-1)*step

        for keys in ab.keys():
            if keys >= start_bin and keys <= end_bin:
                aab[i].append(ab[keys])
        #print('Working... on the bin {}/{}'.format(i,total_bin))


        ## DEAL WITH THE EMPTY MATRIX
        if aab[i] == []:
            aab[i] = [1]*61
            str_aab = [str(i) for i in aab[i]]
            wfile.write('{}\t{}\n'.format(i,'\t'.join(str_aab)))

        else:

            ######################DDD##################
            # EXTRACT ONLY LANDRACE AND VARIETY COLUMNS
            """
            num1 = np.array(aab[i])
            num2 = np.matrix(num1)
            aab[i] = np.mean(num2,0)
            # EXTRACT ONLY LANDRACE AND VARIETY COLUMNS
            lv = aab[i][:,:]
            norm_lv_block = drop_into_num(lv)
            """


            ######################AB###################

            num1 = np.array(aab[i])
            num2 = np.matrix(num1)
            aab[i] = np.mean(num2,0)

            # EXTRACT ONLY LANDRACE AND VARIETY COLUMNS
            
            lv1 = aab[i][:,0:45]
            lv2 = aab[i][:,46:48]
            lv3 = aab[i][:,59:64]
            lv4 = aab[i][:,67:72]
            lv5 = aab[i][:,82::]
            lv_block = np.hstack((lv1,lv2,lv3,lv4,lv5))

            norm_lv_block = drop_into_num(lv_block)


            str_lv_block = [str(i) for i in norm_lv_block]
            wfile.write('{}\t{}\n'.format(i,'\t'.join(str_lv_block)))

        ## BIN ORDER...
        i = i + 1

make_bin_and_write(length_file=args.length,chr=args.chr,ifile=args.input,window=args.window,step=args.step,wfile=args.output)
