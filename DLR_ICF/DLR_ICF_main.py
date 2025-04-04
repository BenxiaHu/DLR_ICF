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
import matplotlib.pyplot as plt

dir = os.path.dirname(__file__)
version_py = os.path.join(dir, "_version.py")
exec(open(version_py).read())

def annotation(format,inputpath,filename,rangeid,bin,outpath,chrsize,outfile):
    # c.matrix(balance=False, as_pixels=True, join=True)[1000:1005, 1000:1005] check it
    bin = int(bin)
    contact = cooler.Cooler(inputpath+'/'+filename+'.cool')
    if format == "balance":
        ### load cooler balanced contact matrix
        ### AD_rep1.cool
        bins = contact.bins()[:]
        pix = contact.pixels()[:]
        weights, metadata = cooler.balance_cooler(contact, rescale_marginals=False)
        pix['weight1'] = weights[pix['bin1_id']]
        pix['weight2'] = weights[pix['bin2_id']]
        input = cooler.annotate(pix, bins)
        del weights,metadata,pix, bins
        # the "weights" as used by cooler are the reciprocals of the "biases"
        # w_i = 1 / beta_i
        input['count'] = input['count'] * input['weight1'] * input['weight2']
        input = input.dropna(subset=['count'])
        input = input[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2','count']]
    elif format == "ICE":
        ### load HiC-Pro ICE normalized contact matrix
        input = contact.matrix(balance=False, as_pixels=True, join=True)[:]
        input = input.dropna(subset=['count'])
        input = input[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2','count']]
    else:
        print("input format is wrong")
    del contact
    #### analyze DLR
    N = rangeid / bin

    input.loc[:, 'bin1_id'] = (input['start1'] / bin).astype('int')
    input.loc[:, 'bin2_id'] = (input['start2'] / bin).astype('int')

    ## remove chrY
    input = input[input['chrom1'] != 'chrY']
    input = input[input['chrom2'] != 'chrY']

    DLRmatrix = input[input['chrom1'] == input['chrom2']]
    DLRmatrix = DLRmatrix[DLRmatrix['bin1_id'] != DLRmatrix['bin2_id']]

    DLRmatrix = DLRmatrix[['chrom1','bin1_id','bin2_id','count']]
    DLRmatrix2 = DLRmatrix[['chrom1','bin2_id','bin1_id','count']]
    DLRmatrix2.rename(columns = {'bin2_id':'bin1_id', 'bin1_id':'bin2_id'}, inplace = True)
    result = pd.concat([DLRmatrix,DLRmatrix2])
    del DLRmatrix2,DLRmatrix
    result['distance'] = result['bin2_id'] - result['bin1_id']
    result = result[abs(result['distance']) >1]
    result['type'] = np.where(abs(result['distance']) <= N, 'local','distal')
    result = result[['chrom1','bin1_id','type','count']]
    result = result.groupby(['chrom1','bin1_id','type'],observed=True).agg('sum')
    result = result.reset_index(level=['chrom1','bin1_id', 'type'])

    #Reshape from long to wide
    result = pd.pivot(result, index=['chrom1','bin1_id'], columns = 'type',values = 'count')
    result = result.reset_index(level=['chrom1','bin1_id'])
    result = result[(result['distal'] > 0)  & (result['local'] > 0)]

    result['DLR_ratio'] = np.log2((result['distal'])/(result['local']))
    result['DLR_ratio'] = np.round(result['DLR_ratio'], decimals=4)
    result['distal'] = np.round(result['distal'], decimals=4)
    result['local'] = np.round(result['local'], decimals=4)
    result['start'] = result['bin1_id'] * bin
    result['end'] = result['start'] + 1 * bin
    result = result[['chrom1','start','end','distal','local','DLR_ratio']]
    print(result)
    chrfile = pd.read_csv(chrsize,sep="\t",header=None)
    chrfile.rename(columns = {0:'chrom1', 1:'size'}, inplace = True)
    result = pd.merge(result,chrfile,on=['chrom1'])

    for i in range(result.shape[0]):
        if(result.iloc[i,2] > result.iloc[i,6]):
            result.iloc[i,2] = result.iloc[i,6]
    result = result[['chrom1','start', 'end','distal','local','DLR_ratio']]
    result.to_csv(outpath+'/'+outfile+'_DLR_ratio.bedgraph',sep="\t",header=False,index=False)
    
    del result

    ##### ICF
    ICFmatrix = input.copy()
    del input
    ICFmatrix['type'] = np.where(ICFmatrix['chrom1'] != ICFmatrix['chrom2'],'inter','intra')
    ICFmatrix1 = ICFmatrix[['chrom1','start1','end1','count','type']]
    ICFmatrix2 = ICFmatrix[['chrom2','start2','end2','count','type']]
    del ICFmatrix
    ICFmatrix1 = ICFmatrix1.rename(columns={'chrom1': 'chrom', 'start1': 'start', 'end1': 'end'})
    ICFmatrix2 = ICFmatrix2.rename(columns={'chrom2': 'chrom', 'start2': 'start', 'end2': 'end'})
    result = pd.concat([ICFmatrix1,ICFmatrix2])
    del ICFmatrix1,ICFmatrix2

    result['count'] = np.round(result['count'], decimals=4)
    result = result.groupby(['chrom','start','end','type'],observed=True).agg('sum')
    result = result.reset_index(level=['chrom','start', 'end','type'])
    result = pd.pivot(result, index=['chrom','start', 'end'], columns = 'type',values = 'count')
    result = result[(result['inter'] > 0)  & (result['intra'] > 0)]
    result = result.reset_index(level=['chrom','start', 'end'])
    result['ICF_ratio'] = np.log2((result['inter'])/(result['inter']+result['intra']))
    result['ICF_ratio'] = np.round(result['ICF_ratio'], decimals=4)
    result['inter'] = np.round(result['inter'], decimals=4)
    result['intra'] = np.round(result['intra'], decimals=4)
    chrfile = pd.read_csv(chrsize,sep="\t",header=None)
    chrfile.rename(columns = {0:'chrom', 1:'size'}, inplace = True)
    result = pd.merge(result,chrfile,on=['chrom'])

    for i in range(result.shape[0]):
        if(result.iloc[i,2] > result.iloc[i,6]):
            result.iloc[i,2] = result.iloc[i,6]
    result = result[['chrom','start', 'end','inter','intra','ICF_ratio']]
    result.to_csv(outpath+'/'+outfile+'_ICF_ratio.bedgraph',sep="\t",header=False,index=False)
    del result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-F', '--format', type=str, default='balance', 
                       choices=['balance', 'ICE'], help='Format of .cool file.')
    parser.add_argument('-I', '--inputpath', dest='inputpath',
                        required=True,
                        help='path of input file')
    parser.add_argument('-f', '--filename', dest='filename',
                        required=True,
                        help='name of input file')
    parser.add_argument('-d', '--distance', dest='distance',
                        required=True, type=int,
                        help='the distance of distal chromation interactions')
    parser.add_argument('-r', '--resolution', dest='resolution',
                        required=True, type=int,
                        help='resolution of contact matrix')
    parser.add_argument('-O', '--outpath', dest='outpath',
                        required=True,
                        help='path of output file')
    parser.add_argument('-c', '--chrsize', dest='chrsize',
                        required=True,
                        help='chromosome size file')
    parser.add_argument('-o', '--outfile', dest='outfile',
                        required=True,
                        help='name of output file')
    parser.add_argument("-V", "--version", action="version",version="DLR_ICF_main {}".format(__version__)\
                      ,help="Print version and exit")
    args = parser.parse_args()
    print('###Parameters:')
    print(args)
    print('###Parameters')
    annotation(args.format,args.inputpath,args.filename,args.distance,args.resolution,args.outpath,args.chrsize,args.outfile)

if __name__ == '__main__':
    main()
