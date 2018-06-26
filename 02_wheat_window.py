import argparse
from collections import defaultdict
import math
import numpy as np


###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input',metavar='File',dest='input',help='Input file',type=open,required=True)
parser.add_argument('-chr','--chromosome',metavar='Str',dest='chr',help='chromosome',type=str,default='chr1A')
parser.add_argument('-l','--length',metavar='File',dest='length',help='length file',type=open,required=True)
parser.add_argument('-w','--window',metavar='Int',dest='window',help='window',type=int,default=10000)
parser.add_argument('-s','--step',metavar='Int',dest='step',help='step',type=int,default=10000)
parser.add_argument('-o','--Output',metavar='File',dest='output',help='Output file',type=argparse.FileType('w'),required=True)
args = parser.parse_args()
###### arguments ######

def look_for_chr_length(length_file,chr):
    for line in length_file:
        line = line.strip().split()
        if line[0] == chr:
            ll = int(line[1])
            return ll


def cal_each_bin(a):
    new_a = []
    col_num = a.shape[1]
    for i in range(0,col_num):
        ref = 0
        col_a = a[:,i]

        ## DELETION ./.
        del_col_a = [i == './.' for i in col_a].count(True)

        col_a_pop_size = 2*(len(col_a) - del_col_a)
        for j in col_a:
            if j == '0/0':
                ref += 2
            elif j == '0/1':
                ref += 1
            else:
                continue
        cal_col_a = float('%.4f'%(ref/(col_a_pop_size+.0000000000000001)))
        new_a.append(cal_col_a)
    return new_a



def make_bin_and_write(length_file,chr,ifile,window,step,wfile):

    ll = look_for_chr_length(length_file,chr)

    """
    1145847 ./.     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     1/1     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/
    1158042 0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     ./.     1/1     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/
    1158958 0/0     1/1     0/0     0/0     0/0     0/0     0/0     0/0     ./.     0/0     0/0     0/0     0/0     1/1     1/1     0/0     0/0     0/0     1/1     0/0     0/0     0/
    1159689 0/0     1/1     0/1     0/0     0/0     0/0     0/0     ./.     0/0     0/0     0/0     0/0     0/0     0/0     1/1     0/0     0/0     0/0     1/1     0/0     0/0     0/
    1159868 0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     1/1     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/
    """

    ## COLLECTING DATA
    ab = {}
    for lline in ifile:
        lline = lline.strip().split()
        ab[int(lline[0])] = np.array([i for i in lline[1::]])


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
            aab[i] = [1.]*61
            str_aab = [str(i) for i in aab[i]]
            wfile.write('{}\t{}\n'.format(i,'\t'.join(str_aab)))

        else:
            # EXTRACT ONLY LANDRACE AND VARIETY COLUMNS
            num1 = np.array(aab[i])

            new_num1 = cal_each_bin(num1)

            # EXTRACT ONLY LANDRACE AND VARIETY COLUMNS
            lv1 = new_num1[0:45]
            lv2 = new_num1[46:48]
            lv3 = new_num1[59:64]
            lv4 = new_num1[67:72]
            lv5 = new_num1[82::]

            lv_block = lv1 + lv2 + lv3 + lv4 + lv5
            str_lv_block = [str(i) for i in lv_block]
            wfile.write('{}\t{}\n'.format(i,'\t'.join(str_lv_block)))

        ## BIN ORDER...
        i = i + 1

make_bin_and_write(length_file=args.length,chr=args.chr,ifile=args.input,window=args.window,step=args.step,wfile=args.output)
