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

def annotation(balanced,inputpath,filename,rangeid,bin,outpath,chrsize,PC,outfile):
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
    input = input[input['bin1_id'] != input['bin2_id']]

    compartment = pd.read_csv(PC,sep="\t",header=None)
    compartment.rename(columns = {0:'chrom1', 1:'start1',2:'end1',3:'PCscore'}, inplace = True)
    compartment['PCid'] = np.where(compartment['PCscore'] <= 0, 'B','A')
    compartment = compartment[['chrom1','start1','end1','PCid']]

    input = pd.merge(input,compartment,how='left',on=['chrom1','start1','end1'])
    input.rename(columns = {'PCid':'left_PCid'}, inplace = True)

    compartment.rename(columns = {'chrom1':'chrom2','start1':'start2','end1':'end2','PCid':'right_PCid'}, inplace = True)
    input = pd.merge(input,compartment,how='left',on=['chrom2','start2','end2'])

    #### analyze DLR
    DLRmatrix = input[input['chrom1'] == input['chrom2']]
    DLRmatrix1 = DLRmatrix[['chrom1','bin1_id','bin2_id','count','left_PCid','right_PCid']]
    DLRmatrix2 = DLRmatrix[['chrom1','bin2_id','bin1_id','count','right_PCid','left_PCid']]
    DLRmatrix1['PCid'] = DLRmatrix1['left_PCid'] + DLRmatrix1['right_PCid']
    DLRmatrix2['PCid'] = DLRmatrix1['right_PCid'] + DLRmatrix1['left_PCid']

    DLRmatrix1.rename(columns = {'chrom1':'chrom'}, inplace = True)
    DLRmatrix2.rename(columns = {'chrom1':'chrom','bin2_id':'bin1_id', 'bin1_id':'bin2_id'}, inplace = True)

    #result = pd.merge(result,compartment,how='left',on=['chrom','bin1_id'])

    result = pd.concat([DLRmatrix1,DLRmatrix2])
    result = result[['chrom','bin1_id','bin2_id','count','PCid']]
    result['distance'] = result['bin2_id'] - result['bin1_id']
    result['type'] = np.where(abs(result['distance']) <= N, 'local','distal')

    result['typeid'] = result['type'] + result['PCid']
    result = result[['chrom','bin1_id','typeid','count']] # 151470258

    result = result.groupby(['chrom','bin1_id','typeid'],observed=True).agg('sum') # 1188240
    result = result.reset_index(level=['chrom','bin1_id', 'typeid'])

    #Reshape from long to wide
    result = pd.pivot(result, index=['chrom','bin1_id'], columns = 'typeid',values = 'count') # 594120
    result = result.reset_index(level=['chrom','bin1_id'])

    result1 = result[['chrom','bin1_id', 'distalAA','localAA']]
    result2 = result[['chrom','bin1_id', 'distalAB','localAA']]
    result3 = result[['chrom','bin1_id', 'distalAA','localAB']]
    result4 = result[['chrom','bin1_id', 'distalAB','localAB']]

    result5 = result[['chrom','bin1_id', 'distalBA','localBA']]
    result6 = result[['chrom','bin1_id', 'distalBB','localBA']]
    result7 = result[['chrom','bin1_id', 'distalBA','localBB']]
    result8 = result[['chrom','bin1_id', 'distalBB','localBB']]

    result1 = result1[(result1['distalAA'] > 0)  & (result1['localAA'] > 0)]
    result2 = result2[(result2['distalAB'] > 0)  & (result2['localAA'] > 0)]
    result3 = result3[(result3['distalAA'] > 0)  & (result3['localAB'] > 0)]
    result4 = result4[(result4['distalAB'] > 0)  & (result4['localAB'] > 0)]
    result5 = result5[(result5['distalBA'] > 0)  & (result5['localBA'] > 0)]
    result6 = result6[(result6['distalBB'] > 0)  & (result6['localBA'] > 0)]
    result7 = result7[(result7['distalBA'] > 0)  & (result7['localBB'] > 0)]
    result8 = result8[(result8['distalBB'] > 0)  & (result8['localBB'] > 0)]

    result1['DLR_ratio'] = np.log2((result1['distalAA'])/(result1['localAA']))
    result1['type'] = 'distalAA' + "_"+ 'localAA'
    result1 = result1[['chrom','bin1_id','DLR_ratio','type']]
    result2['DLR_ratio'] = np.log2((result2['distalAB'])/(result2['localAA']))
    result2['type'] = 'distalAB' + "_"+ 'localAA'
    result2 = result2[['chrom','bin1_id','DLR_ratio','type']]
    result3['DLR_ratio'] = np.log2((result3['distalAA'])/(result3['localAB']))
    result3['type'] = 'distalAA' + "_"+ 'localAB'
    result3 = result3[['chrom','bin1_id','DLR_ratio','type']]
    result4['DLR_ratio'] = np.log2((result4['distalAB'])/(result4['localAB']))
    result4['type'] = 'distalAB' + "_"+ 'localAB'
    result4 = result4[['chrom','bin1_id','DLR_ratio','type']]
    result5['DLR_ratio'] = np.log2((result5['distalBA'])/(result5['localBA']))
    result5['type'] = 'distalBA' + "_"+ 'localBA'
    result5 = result5[['chrom','bin1_id','DLR_ratio','type']]
    result6['DLR_ratio'] = np.log2((result6['distalBB'])/(result6['localBA']))
    result6['type'] = 'distalBB' + "_"+ 'localBA'
    result6 = result6[['chrom','bin1_id','DLR_ratio','type']]
    result7['DLR_ratio'] = np.log2((result7['distalBA'])/(result7['localBB']))
    result7['type'] = 'distalBA' + "_"+ 'localBB'
    result7 = result7[['chrom','bin1_id','DLR_ratio','type']]
    result8['DLR_ratio'] = np.log2((result8['distalBB'])/(result8['localBB']))
    result8['type'] = 'distalBB' + "_"+ 'localBB'
    result8 = result8[['chrom','bin1_id','DLR_ratio','type']]

    result = pd.concat([result1,result2,result3,result4,result5,result6,result7,result8])
    result['start'] = result['bin1_id'] * bin
    result['end'] = (result['bin1_id'] + 1) * bin
    result = result[['chrom','start','end','type','DLR_ratio']]
    chrfile = pd.read_csv(chrsize,sep="\t",header=None)
    chrfile.rename(columns = {0:'chrom', 1:'size'}, inplace = True)
    result = pd.merge(result,chrfile,on=['chrom'])

    for i in range(result.shape[0]):
        if(result.iloc[i,2] > result.iloc[i,5]):
            result.iloc[i,2] = result.iloc[i,5]
    result = result[['chrom','start', 'end','type','DLR_ratio']]
    result.to_csv(outpath+'/'+outfile+'_DLR_ratio.bedgraph',sep="\t",header=False,index=False)
    
    del DLRmatrix2,DLRmatrix,result

    ##### ICF
    ICFmatrix = input.copy()
    ICFmatrix['type'] = np.where(input['chrom1'] != input['chrom2'],'inter','intra')
    ICFmatrix['typeid'] = ICFmatrix['left_PCid']+ ICFmatrix['right_PCid']
    ICFmatrix1 = ICFmatrix[['chrom1','start1','end1','count','type','typeid']]
    ICFmatrix2 = ICFmatrix[['chrom2','start2','end2','count','type','typeid']]
    ICFmatrix1.rename(columns = {'chrom1':'chrom', 'start1':'start', 'end1':'end'}, inplace = True)
    ICFmatrix2.rename(columns = {'chrom2':'chrom', 'start2':'start', 'end2':'end'}, inplace = True)

    result = pd.concat([ICFmatrix1,ICFmatrix2])

    chrfile = pd.read_csv(chrsize,sep="\t",header=None)
    chrfile.rename(columns = {0:'chrom', 1:'size'}, inplace = True)

    result = result.groupby(['chrom','start','end','type','typeid'],observed=True).agg('sum')
    result = result.reset_index(level=['chrom','start', 'end','type','typeid'])
    result['type'] = result['type']+ result['typeid']
    result = result[['chrom','start', 'end','type','count']]
    result = pd.pivot(result, index=['chrom','start', 'end'], columns = 'type',values = 'count')
    result = result.reset_index(level=['chrom','start', 'end'])
    result1 = result[['chrom','start', 'end','interAA','intraAA']]
    result2 = result[['chrom','start', 'end','interAB','intraAB']]
    result3 = result[['chrom','start', 'end','interBA','intraBA']]
    result4 = result[['chrom','start', 'end','interBB','intraBB']]

    result1 = result1[(result1['interAA'] > 0)  & (result1['intraAA'] > 0)]
    result1['ICF_ratio'] = np.log2((result1['interAA'])/(result1['interAA']+result1['intraAA']))
    result1['type'] = 'interAA_intraAA'
    result1 = result1[['chrom','start', 'end','ICF_ratio','type']]

    result2 = result2[(result2['interAB'] > 0)  & (result2['intraAB'] > 0)]
    result2['ICF_ratio'] = np.log2((result2['interAB'])/(result2['interAB']+result2['intraAB']))
    result2['type'] = 'interAB_intraAB'
    result2 = result2[['chrom','start', 'end','ICF_ratio','type']]

    result3 = result3[(result3['interBA'] > 0)  & (result3['intraBA'] > 0)]
    result3['ICF_ratio'] = np.log2((result3['interBA'])/(result3['interBA']+result3['intraBA']))
    result3['type'] = 'interBA_intraBA'
    result3 = result3[['chrom','start', 'end','ICF_ratio','type']]

    result4 = result4[(result4['interBB'] > 0)  & (result4['intraBB'] > 0)]
    result4['ICF_ratio'] = np.log2((result4['interBB'])/(result4['interBB']+result4['intraBB']))
    result4['type'] = 'interBB_intraBB'
    result4 = result4[['chrom','start', 'end','ICF_ratio','type']]

    result = pd.concat([result1,result2,result3,result4])

    #result['ICF_ratio'] = np.round(result['ICF_ratio'], decimals=4)

    result = pd.merge(result,chrfile,on=['chrom'])

    for i in range(result.shape[0]):
        if(result.iloc[i,2] > result.iloc[i,5]):
            result.iloc[i,2] = result.iloc[i,5]
    result = result[['chrom','start', 'end','type','ICF_ratio']]
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
    annotation(args.balanced,args.inputpath,args.filename,args.distance,args.resolution,args.outpath,args.chrsize,args.compartment,args.outfile,args.version)

if __name__ == '__main__':
    main()
