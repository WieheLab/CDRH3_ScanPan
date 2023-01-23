#!/usr/bin/python


#libraries
import xlrd, re
import glob
import os
import string
import sys,math
import numpy,scipy

def CalcSE(posFile,position):
    aaDict=dict()
    seq_counts=0
    with open(posFile,'r') as f:
        for line in f:
            if ">" in line:
                continue
            AA=line.rstrip()[int(position)]
            if AA in aaDict:
                aaDict[AA]+=1
            else:
                aaDict[AA]=1
            seq_counts+=1

    SE=0
    for key in aaDict.keys():
        #print("{}\t{}\t{}".format(key,aaDict[key],float(aaDict[key])/float(seq_counts)))
        Prob=float(aaDict[key])/float(seq_counts)
        SE+=Prob*math.log2(Prob)
    return -SE

def main(posFile):
    seq_len=0
    with open(posFile,'r') as f:
        for line in f:
            if ">" in line:
                continue
            else:
                seq_len=len(line.rstrip())
                break
    for i in range(seq_len):
        SE=CalcSE(posFile,i)
        print("{}\t{:F}".format(i+1,SE))
if __name__=="__main__":
    main(sys.argv[1])
