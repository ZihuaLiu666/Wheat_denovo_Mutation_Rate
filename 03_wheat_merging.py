import argparse
import numpy as np
import pandas as pd


###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input',metavar='File',dest='input',help='Input file',required=True)
parser.add_argument('-r','--correlation_coefficient',metavar='Float',dest='coeff',help='correlation_coefficient',type=float,default=0.1)
parser.add_argument('-o','--Output',metavar='File',dest='output',help='Output file',type=argparse.FileType('w'),required=True)
args = parser.parse_args()
###### arguments ######


def merge(ifile,r_tresh,wfile):
    df = pd.read_table('{}'.format(ifile),header = None)
    i=0
    matrix_ln = df.shape[0]
    matrix = np.matrix(df)[:,1::]

    while i <= matrix_ln-2:
        r = abs(np.corrcoef(matrix[0,:],matrix[1,:]))

        if  np.isnan(r[1,0]) or (r[1,0] > r_tresh):
            fst_cat = np.array((matrix[0,:] + matrix[1,:])/2)
            matrix = np.vstack((fst_cat,matrix[2::,:]))
            i += 1
        else:
            wfile.write('{}\n'.format(i+1))
            matrix = matrix[1::,:]
            i += 1

merge(ifile=args.input,r_tresh=args.coeff,wfile=args.output)
