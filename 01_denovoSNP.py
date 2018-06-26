import argparse
import gzip

###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input',metavar='File',dest='input',help='Input file',type=open,required=True)
parser.add_argument('-s','--size',metavar='Int',dest='size',help='size',type=int,default=10)
parser.add_argument('-g','--gzfile',metavar='File',dest='gzfile',help='Input gzfile file',required=True)
parser.add_argument('-chr','--chromosome',metavar='Str',dest='chr',help='chromosome',type=str,default='chr4A')
parser.add_argument('-wdmaf','--wild_durum_maf',metavar='Int',dest='wdmaf',help='wild_durum_maf',type=float,default=0.01)
parser.add_argument('-lvmaf','--land_var_maf',metavar='Int',dest='lvmaf',help='land_var_maf',type=float,default=0.05)
parser.add_argument('-o','--Output',metavar='File',dest='output',help='Output file',type=argparse.FileType('w'),required=True)
args = parser.parse_args()
###### arguments ######


###################LOCATE THE BLOCK BINS#####################

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

###################LOCATE THE BLOCK BINS#####################

def trans_GT(aa):
    matrix = []
    for i in aa[2::]:
        gt = i.split(':')[0]
        if gt == '0/0':
            matrix.append('1')
        elif gt == '0/1':
            matrix.append('2')
        elif gt == '1/1':
            matrix.append('3')
        elif gt == './.':
            matrix.append('0')
    return [aa[1]] + matrix

def cal_maf(line,wdmaf,lvmaf):

    ## WILD AND DURUM POPULATIONS
    wd1 = [line[56].split(':')[0]]
    wd2 = [i.split(':')[0] for i in line[59:70]]
    wd3 = [j.split(':')[0] for j in line[75:78]]
    wd4 = [k.split(':')[0] for k in line[83:93]]
    wd_pop = wd1 + wd2 + wd3 + wd4

    ## LANDRACE AND VARIETY POPULATIONS
    lv1 = [i1.split(':')[0] for i1 in line[9:44]]
    lv2 = [i2.split(':')[0] for i2 in line[45:53]]
    lv3 = [i3.split(':')[0] for i3 in line[54:56]]
    lv4 = [i4.split(':')[0] for i4 in line[57:59]]
    lv5 = [i5.split(':')[0] for i5 in line[70:75]]
    lv6 = [i6.split(':')[0] for i6 in line[78:83]]
    lv7 = [i7.split(':')[0] for i7 in line[93::]]
    lv_pop = lv1 + lv2 + lv3 + lv4 + lv5 + lv6 + lv7


    ## INITIALIZATION
    wd_pop_size = 2*(len(wd_pop) - wd_pop.count('./.'))
    wd_ref = 0
    for i in wd_pop:
        if i == '0/0':
            wd_ref += 2
        elif i == '0/1':
            wd_ref += 1
        else:
            continue
    cal_wd = wd_ref/wd_pop_size


    lv_pop_size = 2*(len(lv_pop) - lv_pop.count('./.'))
    lv_ref = 0
    for j in lv_pop:
        if j == '0/0':
            lv_ref += 2
        elif j == '0/1':
            lv_ref += 1
        else:
            continue
    cal_lv = lv_ref/lv_pop_size

    # wild maf = 0.01 (ABORT) AND landrace variety maf = 0.05 (KEEP)
    if (cal_wd <= wdmaf or cal_wd >= 1-wdmaf) and 1-lvmaf > cal_lv > lvmaf:
        return line

def wheat_maf_filter(ifile,size,igzfile,wdmaf,lvmaf,wfile,chr):
    bbin = loc_block_bins(ifile,size)

    with gzip.open(igzfile,'rb') as f_in:
        for line in f_in:
            line = line.decode()
            line = line.strip().split()
            if line[0][0:2] == '##':
                continue
            elif line[0][0:2] == '#C':
                continue
            else:
                if line[0] == chr:
                    for ii in bbin:
                        if ii[0] * 1000000 < int(line[1]) < ii[1] * 1000000:
                            new_line = cal_maf(line,wdmaf,lvmaf)

                            if new_line != None:
                                aa = new_line[0:2] + new_line[9:44] + new_line[45:53] + new_line[54:56] + new_line[57:59]\
                                     + new_line[70:75] + new_line[78:83] + new_line[93::]

                                matrix = trans_GT(aa)
                                wfile.write('{}\n'.format('\t'.join(matrix)))
                            else:
                                continue


wheat_maf_filter(ifile=args.input,size=args.size,igzfile=args.gzfile,
                 wdmaf=args.wdmaf,lvmaf=args.lvmaf,wfile=args.output,chr=args.chr)
