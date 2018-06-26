import argparse
import math
from collections import defaultdict

###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input',metavar='File',dest='input',help='Input file',type=open,required=True)
parser.add_argument('-w','--window',metavar='Int',dest='window',help='window',type=int,default=1000000)
parser.add_argument('-chr','--chromosome',metavar='Str',dest='chr',help='chromosome',type=str,default='chr4A')
parser.add_argument('-l','--length',metavar='File',dest='length',help='length file',type=open,required=True)
parser.add_argument('-o','--Output',metavar='File',dest='output',help='Output file',type=argparse.FileType('w'),required=True)
args = parser.parse_args()
###### arguments ######

def look_for_chr_length(length_file,chr):
    for line in length_file:
        line = line.strip().split()
        if line[0] == chr:
            ll = int(line[1])
            return ll

def SNP_density_vcf(ifile,window,length_file,chr,wfile):
    """
    82493   0       2       2       2       2       1       1       2       2       2       2       2       1       1       2       2       2       2       1       1       2       2       2       1       2
    82529   0       2       2       2       2       1       1       2       2       2       2       2       1       1       2       2       2       2       1       1       2       2       2       1       2
    82561   0       2       2       2       2       1       1       2       2       2       2       2       1       1       1       2       2       2       1       1       2       2       2       1       2
    82617   2       2       2       2       2       1       1       1       2       2       1       2       1       1       1       1       2       2       2       1       1       2       2       1       2
    82675   0       2       2       2       2       1       1       2       2       2       2       2       2       1       2       2       2       2       2       1       2       2       2       1       2
    82684   0       2       2       2       2       1       1       2       2       2       2       2       2       1       2       2       2       2       2       1       2       2       2       1       2
    82814   2       2       2       2       2       1       1       1       2       1       2       2       1       1       2       2       2       2       2       1       2       2       2       1       2
    """

    ll = look_for_chr_length(length_file,chr)

    i = 1
    step = window
    total_bin = int(math.ceil(((ll-window)/step)+1))
    bbin = defaultdict(list)

    while i <= total_bin:
        start_bin = 1+(i-1)*step
        end_bin = window + (i-1)*step
        bbin[(start_bin,end_bin)] = []
        i += 1

    # read ifile
    for line in ifile:
        line = line.strip().split()
        snp_pos = int(line[0])

        for i in bbin.keys():
            if i[0] <= snp_pos <= i[1]:
                bbin[i].append(snp_pos)

    # count snp_pos of each bin
    for j in bbin.keys():
        wfile.write('{}\t{}\t{}\n'.format(j[0],j[1],len(bbin[j])))
    wfile.close()


SNP_density_vcf(ifile=args.input,window=args.window,length_file=args.length,chr=args.chr,wfile=args.output)
