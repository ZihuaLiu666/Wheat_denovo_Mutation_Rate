import argparse
from collections import defaultdict
import math
import numpy as np


###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input',metavar='File',dest='input',help='Input file',type=open,required=True)
parser.add_argument('-chr','--chromosome',metavar='Str',dest='chr',help='chromosome',type=str,default='chr1A')
parser.add_argument('-l','--length',metavar='File',dest='length',help='length file',type=open,required=True)
parser.add_argument('-w','--window',metavar='Int',dest='window',help='window',type=int,default=1000000)
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

def make_bin_and_write(length_file,chr,ifile,window,wfile):

    """
    77488   0       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1
    77617   1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1
    77641   1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1
    77662   1       1       3       1       1       1       1       1       3       3       3       1       3       3       1       1       3       1       3       1       1       1
    77663   1       1       3       1       1       1       1       1       3       3       3       1       3       3       1       1       3       1       3       1       1       1
    77856   1       1       3       1       1       1       1       1       3       3       3       1       3       3       1       1       3       1       3       1       1       1
    77875   1       1       3       1       1       0       1       1       3       3       3       1       3       3       1       1       3       1       3       1       1       1
    """

    ## CREAT BIN DIC
    ll = look_for_chr_length(length_file,chr)

    step = window

    i = 1
    total_bin = int(math.ceil(((ll-window)/step)+1))
    bbin = defaultdict(list)
    while i <= total_bin:
        start_bin = 1+(i-1)*step
        end_bin = window + (i-1)*step
        bbin[(start_bin,end_bin)] = []
        i += 1

    ## COLLECTING DATA
    #aab = defaultdict(list)
    for lline in ifile:
        lline = lline.strip().split()

        snp_pos = int(lline[0])
        snp_pos_info = np.array([int(i) for i in lline[1::]])
        for i in bbin.keys():
            if i[0] <= snp_pos <= i[1]:
                bbin[i].append(snp_pos_info)


    for j in bbin.keys():
        ## DEAL WITH THE EMPTY MATRIX
        if bbin[j] == []:
            bbin[j] = [1]*61
            str_bbin = [str(j) for j in bbin[j]]
            wfile.write('{}\t{}\n'.format(int(j[1]/window),'\t'.join(str_bbin)))

        else:
            if chr in ['chr1D','chr3D','chr3D','chr4D','chr5D','chr6D','chr7D']:
                ######################DDD##################
                # EXTRACT ONLY LANDRACE AND VARIETY COLUMNS
                num1 = np.array(bbin[j])
                num2 = np.matrix(num1)
                bbin[j] = np.mean(num2,0)

                lv = bbin[j][:,:]
                norm_lv_block = drop_into_num(lv)

                str_lv_block = [str(j) for j in norm_lv_block]
                wfile.write('{}\t{}\n'.format(int(j[1]/window),'\t'.join(str_lv_block)))
            else:
                ###################ABABABAB##################
                num1 = np.array(bbin[j])
                num2 = np.matrix(num1)
                bbin[j] = np.mean(num2,0)

                # EXTRACT ONLY LANDRACE AND VARIETY COLUMNS
                lv1 = bbin[j][:,0:45]
                lv2 = bbin[j][:,46:48]
                lv3 = bbin[j][:,59:64]
                lv4 = bbin[j][:,67:72]
                lv5 = bbin[j][:,82::]
                lv_block = np.hstack((lv1,lv2,lv3,lv4,lv5))

                norm_lv_block = drop_into_num(lv_block)

                str_lv_block = [str(i) for i in norm_lv_block]
                wfile.write('{}\t{}\n'.format(int(j[1]/window),'\t'.join(str_lv_block)))

make_bin_and_write(length_file=args.length,chr=args.chr,ifile=args.input,window=args.window,wfile=args.output)
