import pandas as pd
import matplotlib.pyplot as plt
import argparse
import gzip

###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input',metavar='File',dest='input',help='Input file',required=True)
args = parser.parse_args()
###### arguments ######

def denovo_plot(ifile):
    df = pd.read_table(ifile)
    x = df.ix[:,0]
    y = df.ix[:,1]
    size = df.ix[:,2]*10
    ccolor = df.ix[:,3]

    plt.scatter(x,y,s=size,c=ccolor,cmap="coolwarm",alpha=0.6)
    plt.show()

denovo_plot(ifile=args.input)
