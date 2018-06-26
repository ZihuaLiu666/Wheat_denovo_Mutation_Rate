import argparse

###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input',metavar='File',dest='input',help='Input file',type=open,required=True)
parser.add_argument('-ts','--tag_snp',metavar='File',dest='ts',help='Input tag file',type=open,required=True)
parser.add_argument('-o','--Output',metavar='File',dest='output',help='Output file',type=argparse.FileType('w'),required=True)
args = parser.parse_args()
###### arguments ######

def collect_tag_data(tgfile):
    tglist = []
    for line in tgfile:
        line = line.strip().split()
        tglist.append(int(line[3]))
    return tglist


def cal_maf(tgfile,ifile,wfile):
    """
    82493   0       2       2       2       2       1       1       2       2       2       2       2       1       1       2       2
    82529   0       2       2       2       2       1       1       2       2       2       2       2       1       1       2       2
    82561   0       2       2       2       2       1       1       2       2       2       2       2       1       1       1       2
    82617   2       2       2       2       2       1       1       1       2       2       1       2       1       1       1       1
    82675   0       2       2       2       2       1       1       2       2       2       2       2       2       1       2       2
    82684   0       2       2       2       2       1       1       2       2       2       2       2       2       1       2       2
    """
    tglist = collect_tag_data(tgfile)

    for line in ifile:
        line = line.strip().split()
        if int(line[0]) in tglist:
            wfile.write('{}\t{}\n'.format(int(line[0]),'\t'.join(line[1::])))

cal_maf(tgfile=args.ts,ifile=args.input,wfile=args.output)
