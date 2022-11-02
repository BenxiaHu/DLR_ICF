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

dir = os.path.dirname(os.path.abspath(__file__))
version_py = os.path.join(dir, "_version.py")
exec(open(version_py).read())

def annotation(balanced,inputpath,filename,rangeid,bin,outpath,chrsize,outfile):
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

    #### analyze DLR
    N = rangeid / bin
    input['bin1_id'] = (input['start1'] / bin).astype('int')
    input['bin2_id'] = (input['start2'] / bin).astype('int')

    DLRmatrix = input[input['chrom1'] == input['chrom2']]
    DLRmatrix = DLRmatrix[DLRmatrix['bin1_id'] != DLRmatrix['bin2_id']]   

    DLRmatrix = DLRmatrix[['chrom1','bin1_id','bin2_id','count']]
    DLRmatrix2 = DLRmatrix[['chrom1','bin2_id','bin1_id','count']]
    DLRmatrix2.rename(columns = {'bin2_id':'bin1_id', 'bin1_id':'bin2_id'}, inplace = True)
    result = pd.concat([DLRmatrix,DLRmatrix2])
    result['distance'] = result['bin2_id'] - result['bin1_id']
    result['type'] = np.where(abs(result['distance']) <= N, 'local','distal')
    result = result[['chrom1','bin1_id','type','count']] # 151470258
    result = result.groupby(['chrom1','bin1_id','type'],observed=True).agg('sum') # 1188240
    result = result.reset_index(level=['chrom1','bin1_id', 'type'])

    #Reshape from long to wide
    result = pd.pivot(result, index=['chrom1','bin1_id'], columns = 'type',values = 'count') # 594120
    result = result.reset_index(level=['chrom1','bin1_id'])
    result = result[(result['distal'] > 0)  & (result['local'] > 0)] #280189

    result['DLR_ratio'] = np.log2((result['distal'])/(result['local']))
    #result['DLR_ratio'] = np.round(result['DLR_ratio'], decimals=4)
    result['start'] = result['bin1_id'] * bin
    result['end'] = (result['bin1_id'] + 1) * bin
    result = result[['chrom1','start','end','DLR_ratio']]
    chrfile = pd.read_csv(chrsize,sep="\t",header=None)
    chrfile.rename(columns = {0:'chrom1', 1:'size'}, inplace = True)
    result = pd.merge(result,chrfile,on=['chrom1'])

    for i in range(result.shape[0]):
        if(result.iloc[i,2] > result.iloc[i,4]):
            result.iloc[i,2] = result.iloc[i,4]
    result = result[['chrom1','start', 'end','DLR_ratio']]
    result.to_csv(outpath+'/'+outfile+'_DLR_ratio.bedgraph',sep="\t",header=False,index=False)
    
    del DLRmatrix2,DLRmatrix,result

    ##### ICF
    ICFmatrix = input
    ICFmatrix['type'] = np.where(input['chrom1'] != input['chrom2'],'inter','intra')
    ICFmatrix1 = ICFmatrix[['chrom1','start1','end1','count','type']]  
    ICFmatrix2 = ICFmatrix[['chrom2','start2','end2','count','type']]
    ICFmatrix1.rename(columns = {'chrom1':'chrom', 'start1':'start', 'end1':'end'}, inplace = True)
    ICFmatrix2.rename(columns = {'chrom2':'chrom', 'start2':'start', 'end2':'end'}, inplace = True)
    result = pd.concat([ICFmatrix1,ICFmatrix2])
    chrfile = pd.read_csv(chrsize,sep="\t",header=None)
    chrfile.rename(columns = {0:'chrom', 1:'size'}, inplace = True)

    #for i in range(chrfile.shape[0]):
    #ICF = result[result['chrom'] == chrfile.iloc[i,0]]
    
    result = result.groupby(['chrom','start','end','type'],observed=True).agg('sum')
    result = result.reset_index(level=['chrom','start', 'end','type'])
    result = pd.pivot(result, index=['chrom','start', 'end'], columns = 'type',values = 'count')
    result = result[(result['inter'] > 0)  & (result['intra'] > 0)]
    result = result.reset_index(level=['chrom','start', 'end'])
    result = result[(result['inter'] > 0)  & (result['intra'] > 0)]
    result['ICF_ratio'] = np.log2((result['inter'])/(result['inter']+result['intra']))
    #result['ICF_ratio'] = np.round(result['ICF_ratio'], decimals=4)

    result = pd.merge(result,chrfile,on=['chrom'])

    for i in range(result.shape[0]):
        if(result.iloc[i,2] > result.iloc[i,6]):
            result.iloc[i,2] = result.iloc[i,6]
    result = result[['chrom','start', 'end','ICF_ratio']]
    result.to_csv(outpath+'/'+outfile+'_ICF_ratio.bedgraph',sep="\t",header=False,index=False)
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
    parser.add_argument("-V", "--version", action="version",version="DLR_ICF_comparison {}".format(__version__)\
                      ,help="Print version and exit")
    args = parser.parse_args()
    print('###Parameters:')
    print(args)
    print('###Parameters')
    annotation(args.balanced,args.inputpath,args.filename,args.distance,args.resolution,args.outpath,args.chrsize,args.outfile,args.version)

if __name__ == '__main__':
    main()
