import argparse

###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input',metavar='File',dest='input',help='Input file',type=open,required=True)
parser.add_argument('-f','--fourfold',metavar='File',dest='fourfold',help='Input fourfold file',type=open,required=True)
parser.add_argument('-o','--Output',metavar='File',dest='output',help='Output file',type=argparse.FileType('w'),required=True)
args = parser.parse_args()
###### arguments ######
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


def make_pair(wildlist, landlist):
    wildlist = bigRemove(0,wildlist)
    landlist = bigRemove(0,landlist)
    pair = []
    for i in wildlist:
        for j in landlist:
            pair.append(pair2divergency((i,j)))
    return sum(pair)/(len(wildlist)*len(landlist))

# SNP data set is filtered
def collect_data(ffile,ifile,wfile):

    ######################################
    lrList = []
    for line in ffile:
        line = line.strip().split()
        ffpos = int(line[1])
        lrList.append(ffpos)
    # four fold snp position
    ######################################

    divList = []
    for line in ifile:
        line  = line.strip().split()
        snpos = int(line[0])
        if snpos in lrList:
            wildlist = [int(i) for i in (line[49:57] + line[58:60] + line[73:83])]
            landlist = [int(j) for j in (line[1:46] + line[47:49] + line[60:65] + line[68:73] + line[83::])]
            divList.append(make_pair(wildlist, landlist))
    tot = sum(divList)
    wfile.write('{}\t{}'.format(tot, len(lrList)))

collect_data(ffile=args.fourfold,ifile=args.input,wfile=args.output)
