#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
import sys, getopt
import os.path
import argparse
import cooler
import h5py
import math

dir = os.path.dirname(__file__)
version_py = os.path.join(dir, "_version.py")
exec(open(version_py).read())

def annotation(balanced,inputpath,filename,bin,outpath,chrsize,PC,outfile):
    #balanced = int(balanced)
    bin = int(bin)
    if balanced:
        ### load balanced contact matrix
        contact = cooler.Cooler(inputpath+'/'+filename+'_'+str(bin)+'.mcool::resolutions/'+str(bin))
    else:
        ### load ICE/KR normalized contact matrix #SCA-Veh_iced_100000.cool
        contact = cooler.Cooler(inputpath+'/'+filename+'_'+str(bin)+'.cool')
    bins = contact.bins()[:]
    pix = contact.pixels()[:]
    input = cooler.annotate(pix, bins)

    ##### ICF

    input['bin1_id'] = (input['start1'] / bin).astype('int')
    input['bin2_id'] = (input['start2'] / bin).astype('int')
    input = input[input['bin1_id'] != input['bin2_id']]  ## remove diagol values
    input = input[['chrom1','start1','end1','chrom2','start2','end2','count']]
    compartment = pd.read_csv(PC,sep="\t",header=None)
    compartment.rename(columns = {0:'chrom1', 1:'start1',2:'end1',3:'PCscore'}, inplace = True)
    compartment['PCid'] = np.where(compartment['PCscore'] <= 0, 'B','A')
    compartment = compartment[['chrom1','start1','end1','PCid']]

    input = pd.merge(input,compartment,how='left',on=['chrom1','start1','end1'])
    input.rename(columns = {'PCid':'left_PCid'}, inplace = True)

    compartment.rename(columns = {'chrom1':'chrom2','start1':'start2','end1':'end2','PCid':'right_PCid'}, inplace = True)
    input = pd.merge(input,compartment,how='left',on=['chrom2','start2','end2'])

    input['type'] = np.where(input['chrom1'] != input['chrom2'],'inter','intra')

    ICFmatrix1 = input.copy()
    ICFmatrix2 = input[['chrom2','start2','end2','chrom1','start1','end1','count','right_PCid','left_PCid','type']]
    ICFmatrix2.rename(columns = {'chrom2':'chrom1', 'start2':'start1', 'end2':'end1','chrom1':'chrom2', 'start1':'start2', 'end1':'end2','right_PCid':'left_PCid','left_PCid':'right_PCid'}, inplace = True)

    result = pd.concat([ICFmatrix1,ICFmatrix2])
 
    chrfile = pd.read_csv(chrsize,sep="\t",header=None)
    chrfile.rename(columns = {0:'chrom1', 1:'size'}, inplace = True)

    result = pd.merge(result,chrfile,on=['chrom1'])

    for i in range(result.shape[0]):
        if(result.iloc[i,2] > result.iloc[i,10]):
            result.iloc[i,2] = result.iloc[i,10]
    chrfile.rename(columns = {0:'chrom2', 1:'size'}, inplace = True)
    for i in range(result.shape[0]):
        if(result.iloc[i,5] > result.iloc[i,10]):
            result.iloc[i,5] = result.iloc[i,10]
    result = result[['chrom1','start1','end1','chrom2','start2','end2','count','right_PCid','left_PCid','type']]
    result.to_csv(outpath+'/'+outfile+'_ICF_AB.bed',sep="\t",header=True,index=False)
    del ICFmatrix,ICFmatrix1,ICFmatrix2,result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--balanced', dest='balanced',
                        default=False,action='store_false',
                        help='type of contact matrix')
    parser.add_argument('-I', '--inputpath', dest='inputpath',
                        required=True,
                        help='path of input file')
    parser.add_argument('-f', '--filename', dest='filename',
                        required=True,
                        help='name of input file')
    parser.add_argument('-r', '--resolution', dest='resolution',
                        required=True, type=int,
                        help='resolution of contact matrix')
    parser.add_argument('-O', '--outpath', dest='outpath',
                        required=True,
                        help='path of output file')
    parser.add_argument('-c', '--chrsize', dest='chrsize',
                        required=True,
                        help='chromosome size file')
    parser.add_argument('-p', '--compartment', dest='compartment',
                        required=True,
                        help='compartment file')
    parser.add_argument('-o', '--outfile', dest='outfile',
                        required=True,
                        help='name of output file')
    parser.add_argument("-V", "--version", action="version",version="DLR_ICF_comparison {}".format(__version__)\
                      ,help="Print version and exit")
    args = parser.parse_args()
    print('###Parameters:')
    print(args)
    print('###Parameters')
    annotation(args.balanced,args.inputpath,args.filename,args.resolution,args.outpath,args.chrsize,args.compartment,args.outfile,args.version)

if __name__ == '__main__':
    main()
