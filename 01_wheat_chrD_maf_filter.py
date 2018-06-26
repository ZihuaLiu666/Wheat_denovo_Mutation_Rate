import argparse
import gzip

###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input',metavar='File',dest='input',help='Input file',required=True)
parser.add_argument('-chr','--chromosome',metavar='Str',dest='chr',help='chromosome',type=str,default='chr1D')
parser.add_argument('-lvmaf','--land_var_maf',metavar='float',dest='lvmaf',help='land_var_maf',type=float,default=0.00)
parser.add_argument('-o','--Output',metavar='File',dest='output',help='Output file',type=argparse.FileType('w'),required=True)
args = parser.parse_args()
###### arguments ######

def cal_maf(line,lvmaf):


    ## LANDRACE AND VARIETY POPULATIONS
    lv1 = [i1.split(':')[0] for i1 in line[9:44]]
    lv2 = [i2.split(':')[0] for i2 in line[45]]
    lv3 = [i3.split(':')[0] for i3 in line[47:54]]
    lv4 = [i4.split(':')[0] for i4 in line[55:64]]
    lv5 = [i5.split(':')[0] for i5 in line[68::]]

    lv_pop = lv1 + lv2 + lv3 + lv4 + lv5


    ## INITIALIZATION

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

    if 1-lvmaf >= cal_lv >= lvmaf:
        return line

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



def wheat_maf_filter(ifile,lvmaf,wfile,chr):
    with gzip.open(ifile,'rb') as f_in:
        for line in f_in:
            line = line.decode()
            line = line.strip().split()
            if line[0][0:2] == '##':
                continue
            elif line[0][0:2] == '#C':
                continue
            else:
                if line[0] == chr:
                    new_line = cal_maf(line,lvmaf)

                    if new_line != None:
                        aa = new_line[0:2] + new_line[9:44] + [new_line[45]] + line[47:54] + new_line[55:64] + new_line[68::]
                        matrix = trans_GT(aa)
                        wfile.write('{}\n'.format('\t'.join(matrix)))
                    else:
                        continue


wheat_maf_filter(ifile=args.input,lvmaf=args.lvmaf,wfile=args.output,chr=args.chr)
