import argparse

###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input',metavar='File',dest='input',help='Input file',type=open,required=True)
parser.add_argument('-s','--size',metavar='Int',dest='size',help='size',type=int,default=10)
args = parser.parse_args()
###### arguments ######
def loc_block_bins(ifile,size):
    binList = []
    for line in ifile:
        line = line.strip().split()
        binList.append(int(line[0]))

    bbin = []
    for i in range(0,len(binList)-1):
        if abs(binList[i] - binList[i+1]) >= size:
            bbin.append((binList[i],binList[i+1]))

    return bbin

loc_block_bins(ifile=args.input,size=args.size)
