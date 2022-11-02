#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
import sys, getopt
import os.path
import argparse
import scipy.stats
import statsmodels.stats.multitest

dir = os.path.dirname(os.path.abspath(__file__))
version_py = os.path.join(dir, "_version.py")
exec(open(version_py).read())

def DLR_ICF_comparison(inputpath,caseid,ctrlid,bin,outpath,outfile):

    case = pd.read_csv(inputpath+"/"+caseid,sep="\t",header=None)
    case.rename(columns = {0:'chr', 1:'start', 2:'end'}, inplace = True)
    ctrl = pd.read_csv(inputpath+"/"+ctrlid,sep="\t",header=None)
    ctrl.rename(columns = {0:'chr', 1:'start', 2:'end'}, inplace = True)
    result = pd.merge(case,ctrl,on=['chr','start','end'])
    result.rename(columns = {3:'casescore', 4:'ctrlscore'}, inplace = True)
    result.fillna(0,inplace=True)
    result['diff'] = result['casescore'] - result['ctrlscore']
    result['zscore'] = scipy.stats.zscore(result['diff'])
    result['p'] = scipy.stats.norm.sf(abs(result['zscore']))*2
    result['fdr'] = statsmodels.stats.multitest.fdrcorrection(result['p'])[1]
    #if sigid:
    #    topID = result.nlargest(round(result.shape[0] * n), ['zscore'])
    #    bottomID = result.nsmallest(round(result.shape[0] * n), ['zscore'])
    #    result = pd.concat([topID,bottomID])
    result.to_csv(outpath+'/'+outfile+'.bed',sep="\t",header=False,index=False)
    result[['chr','start','end','zscore']].to_csv(outpath+'/'+outfile+'_zscore.bedgraph',sep="\t",header=False,index=False)
    result[['chr','start','end','diff']].to_csv(outpath+'/'+outfile+'_difference.bedgraph',sep="\t",header=False,index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputpath', dest='inputpath',
                        required=True,
                        help='path of input file')
    parser.add_argument('-t', '--treatment', dest='treatment',
                        required=True,
                        help='name of treatment file')
    parser.add_argument('-c', '--control', dest='control',
                        required=True,
                        help='name of control file')
    parser.add_argument('-r', '--resolution', dest='resolution',
                        required=True,type=int,
                        help='resolution of contact matrix')
    parser.add_argument('-O', '--outpath', dest='outpath',
                        required=True,type=float,
                        help='path of output file')
    parser.add_argument('-o', '--outfile', dest='outfile',
                        required=True,
                        help='name of output file')
    parser.add_argument("-V", "--version", action="version",version="DLR_ICF_comparison {}".format(__version__)\
                      ,help="Print version and exit")
    args = parser.parse_args()
    print('###Parameters:')
    print(args)
    print('###Parameters')
    DLR_ICF_comparison(args.inputpath,args.treatment,args.control,args.resolution,args.outpath,args.outfile,args.version)

if __name__ == '__main__':
    main()
