import argparse
import numpy as np
import pandas as pd


###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input',metavar='File',dest='input',help='Input file',required=True)
parser.add_argument('-r1','--correlation_coefficient1',metavar='Float',dest='coeff1',help='correlation_coefficient1',type=float,default=0.7)
parser.add_argument('-r2','--correlation_coefficient2',metavar='Float',dest='coeff2',help='correlation_coefficient2',type=float,default=0.4)
parser.add_argument('-o','--Output',metavar='File',dest='output',help='Output file',type=argparse.FileType('w'),required=True)
args = parser.parse_args()
###### arguments ######


def LD_detector(ifile,r1_tresh,r2_tresh,wfile):
    df = pd.read_table('{}'.format(ifile),header = None)
    i=0
    matrix_ln = df.shape[0]
    matrix = np.matrix(df)[:,1::]
    matrix = np.hstack((matrix[:,0:45],matrix[:,46:48],matrix[:,59:64],matrix[:,67:72],matrix[:,82:86]))

    while i <= matrix_ln-3:
        r1 = abs(np.corrcoef(matrix[0,:],matrix[1,:]))

        if  np.isnan(r1[1,0]) or (r1[1,0] > r1_tresh):
            r2 = abs(np.corrcoef(matrix[0,:],matrix[2,:]))
            r3 = abs(np.corrcoef(matrix[1,:],matrix[2,:]))

            if (np.isnan(r2[1,0]) or (r2[1,0] > r2_tresh)) and (np.isnan(r3[1,0]) or (r3[1,0] > r2_tresh)):
                matrix = np.vstack((matrix[0,:],matrix[2::,:]))
                i += 1

            else:
                wfile.write('{}\n'.format(i+1))
                matrix = matrix[1::,:]
                i += 1
        else:
            wfile.write('{}\n'.format(i+1))
            matrix = matrix[1::,:]
            i += 1

LD_detector(ifile=args.input,r1_tresh=args.coeff1,r2_tresh=args.coeff2,wfile=args.output)
