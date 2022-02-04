#!/usr/bin/python3

#libraries
import xlrd, re
import glob
import os
import string
import sys
import numpy,scipy
import argparse
import gzip

from matplotlib import rcParams
rcParams['font.family']='monospace'
import matplotlib.pyplot as plt

MIN_PYTHON=(3,5)
if sys.version_info < MIN_PYTHON:
    sys.exit("Python %s.%s or later is required.\n" % MIN_PYTHON)


class SeqDistCalc():
    def __init__(self,filename,matrixfile):
        count=0
        matrix=self.read_Matrix(matrixfile)
        if filename[-2:]=="gz":
            self.distfile=filename.replace(".gz",".dist")
            self.processGZfile(filename,matrix,self.distfile)
        else:
            self.distfile=filename+".dist"
            self.processTXTfile(filename,matrix,self.distfile)
        return
    def read_Matrix(self,matrix_file):
        scoring_matrix=dict()
        with open(matrix_file,"rb") as f:
            count=0
            while True:
                line=f.readline()
                if count==0:
                    count+=1
                    continue
                if not line:
                    break
                line=line.decode("utf-8").rstrip().split()
                scoring_matrix[line[0]]=line[1:]
                count+=1
        return scoring_matrix
    def score_Seq(self,seq,matrix):
        all_dist=[]
        #for cutoff in [x/4.0 for x in range(-12,9)]:
        for cutoff in [x/5.0 for x in range(-15,1)]:
            dist=20
            for i,nt in enumerate(seq):
                if float(matrix[nt][i])>float(cutoff):
                    dist-=1
            all_dist.append(dist)
        distance="\t".join(["{}".format(i) for i in all_dist])
        return distance
    def processGZfile(self,filename,matrix,outname):
        with open(outname,'w') as out:
            with gzip.open(filename, 'rb') as f:
                while True:
                    line=f.readline()
                    if not line:
                        break
                    seq,Dvalue=line.rstrip().split()
                    score=self.score_Seq(seq.decode("utf-8")[1:21],matrix)
                    out.write("{}\t{}\n".format(seq.decode("utf-8"),score))

    def processTXTfile(self,filename,matrix,outname):
        with open(outname,'w') as out:
            with open(filename, "r") as f:
                while True:
                    line=f.readline()
                    if not line:
                        break
                    seq=line.rstrip()
                    if len(seq)==20:
                        score=self.score_Seq(seq,matrix)
                        out.write("{}\t{}\n".format(seq,score))

class PlotHist():
    def __init__(self,filename,outname,scaled_value):
        self.threshold_list=[x/5.0 for x in range(-15,1)]
        distanceHash=dict()
        for T in self.threshold_list:
            distanceHash[T]=numpy.zeros([1,21])
        with open(filename, "r") as infile:
            while True:
                line=infile.readline()
                if not line:
                    break
                dataline=line.rstrip().split()
                for i,value in enumerate(dataline[1:]):
                    distanceHash[self.threshold_list[i]][0,int(value)]+=1
                    
        heatmap=numpy.zeros([len(self.threshold_list),21])
        scaledheatmap=numpy.zeros([len(self.threshold_list),21])
        for i,T in enumerate(self.threshold_list):
            fullCount=numpy.sum(distanceHash[T])
            #print("standard scaling {}".format(fullCount))
            for pos,n in enumerate(distanceHash[T][0,:]):
                heatmap[i,pos]=n/float(fullCount)
                scaledheatmap[i,pos]=n/float(scaled_value)
        self.plot_heatmap(outname,heatmap)
        writeTable(outname.replace(".pdf",".table.txt"),heatmap, self.threshold_list)
        self.plot_cum_heatmap(outname.replace(".pdf",".cumulative.pdf"),heatmap)

        self.plot_heatmap(outname.replace(".pdf",".scaled.pdf"),scaledheatmap)
        writeTable(outname.replace(".pdf",".scaled.table.txt"),scaledheatmap, self.threshold_list)   
        self.plot_cum_heatmap(outname.replace(".pdf",".scaled.cumulative.pdf"),scaledheatmap)

        self.plot_legend(outname.replace(".pdf","legend.pdf"))

    def GenerateBins(self,data_list):
        bins=numpy.zeros([1,21])
        data=numpy.array(data_list)
        for i in range(0,21):
            bins[0,i]=len(data[data==i])/float(len(data_list))
        return bins
    def color_lookup(self,v):
        c=[247,252,240]
        if   v<=1/1000000000.0:
            #c=[247,252,240]
            c=[235,255,230]
        elif v<=1/100000000.0:
            c=[224,243,219]
        elif v<=1/10000000.0:
            c=[204,235,197]
        elif v<=1/1000000.0:
            c=[168,221,181]
        elif v<=1/100000.0:
            c=[123,204,196]
        elif v<=1/10000.0:
            c=[78,179,211]
        elif v<=1/1000.0:
            c=[43,140,190]
        elif v<=1/100.0:
            c=[8,104,172]
        else:
            c=[8,64,129]

        return [i/255.0 for i in c]
    def plot_heatmap(self,outpdf,data):
        plt.clf()
        fig,ax=plt.subplots(figsize=(10,10))
        ax.axis('tight')

        ax.set_ylim(-0.5,data.shape[0]-0.5)
        ax.set_xlim(-0.5,data.shape[1]-0.5)
        #ax.imshow(data)
        for i,vec_thresh in enumerate(data):
            for n,data in enumerate(vec_thresh):
                if data==0:
                    color=[1,1,1]
                else:
                    color=self.color_lookup(data)
                rect=plt.Rectangle((n-0.5,i-0.5),1,1,fc=color,ec="black")
                ax.add_patch(rect)
        ax.set_yticks(range(0,len(self.threshold_list)))
        ax.set_yticklabels(["{0:.2f}".format(x) for x in self.threshold_list])
        ax.set_xticks(range(0,21))
        ax.grid(b=None)
        plt.grid(b=None)
        plt.xlabel("Amino Acid Distance")
        plt.ylabel("Matrix Threshold Value")
        plt.savefig(outpdf,bbox_inches='tight')
        plt.close("all")
    def plot_cum_heatmap(self,outpdf,data):
        plt.clf()
        fig,ax=plt.subplots(figsize=(10,10))
        ax.axis('tight')
        ax.set_ylim(-0.5,data.shape[0]-0.5)
        ax.set_xlim(-0.5,data.shape[1]-0.5)
        for i,vec_thresh in enumerate(data):
            cum_data=0
            for n,data in enumerate(vec_thresh):
                cum_data+=data
                if cum_data==0:
                    color=[1,1,1]
                else:
                    color=self.color_lookup(cum_data)
                rect=plt.Rectangle((n-0.5,i-0.5),1,1,fc=color,ec="black")
                ax.add_patch(rect)
        ax.set_yticks(range(0,len(self.threshold_list)))
        ax.set_yticklabels(["{0:.2f}".format(x) for x in self.threshold_list])
        ax.set_xticks(range(0,21))
        ax.set_xticklabels(["$\leq${}".format(x) for x in range(0,21)])
        ax.grid(b=None)
        plt.grid(b=None)
        plt.xlabel("Amino Acid Distance (cumulative)")
        plt.ylabel("Matrix Threshold Value")
        plt.savefig(outpdf,bbox_inches='tight')
        plt.close("all")

    def plot_legend(self,outpdf):
        fig,ax=plt.subplots(figsize=(3,5))
        ax.axis('off')
        ax.set_ylim(-0.5,10.5)
        ax.set_xlim(-0.5,6.5)
        values=[10,100,1000,10000,100000,1000000,10000000,100000000,1000000000,0]
        for i,data in enumerate(values):
            if data==0:
                color=[1,1,1]
                ax.text(1.5,i+0.5,"Not Detected".format(data))
            else:
                color=self.color_lookup(0.9/data)
                if i==8:
                    ax.text(1.5,i+0.5,"$\leq$ 1 in {:,d}".format(data))
                else:
                    ax.text(1.5,i+0.5,"> 1 in {:,d}".format(values[i+1]))
            rect=plt.Rectangle((0.0,i),1,1,fc=color,ec="black")
            ax.add_patch(rect)
        plt.savefig(outpdf,bbox_inches='tight')
        plt.close("all")

def writeTable(outname,data,Vect):#need to reverse the data
    thresholdVect=Vect[::-1]
    with open(outname,'w') as outfile:
        for i,vec_thresh in enumerate(data[::-1]):
            outfile.write("{}\t".format(thresholdVect[i]))
            for n,values in enumerate(vec_thresh):
                outfile.write("{:0.4e}\t".format(values))
            outfile.write("\n")
        outfile.write("\n")
        outfile.write(" ")
        for i in range(len(data[0])):
            outfile.write("\t{}".format(i))
        outfile.write("\n")
                        
def parse_args():
    # Construct an argument parser
    all_args = argparse.ArgumentParser()
    
    # Add arguments to the parser
    all_args.add_argument("-m","--matrix",required=True,help="matrix data")
    all_args.add_argument("-g","--gene",required=False,help="gene")
    all_args.add_argument("-s","--sequences",required=True,help="set of sequences")
    all_args.add_argument("-c","--cutoff",required=False,default=-0.2,type=float,help="cut off for a hit")
    all_args.add_argument("-o","--pdfname",required=False,default=[],help="name of the final pdf file")
    all_args.add_argument("--scale",required=False,default=85149053,type=int,help="value to scale the data by")
    #scaled_value=85149053
    args = vars(all_args.parse_args())
    return args

def main(args):
    
    matrixfile=args["matrix"]
    sequencefile=args["sequences"]
    #cutoffvalue=args["cutoff"]
    SDC=SeqDistCalc(sequencefile,matrixfile)

    if args["pdfname"]:
        fn,f_ext=os.path.splitext(args["pdfname"])
        if not f_ext:
            fn=fn+".pdf"
        elif ".pdf" not in f_ext:
            fn=args["pdfname"]+".pdf"
        PlotHist(SDC.distfile,fn,args["scale"])
    else:
        PlotHist(SDC.distfile,SDC.distfile.replace(".txt","").replace(".dist",".pdf"),args["scale"])

if __name__=="__main__":
    args=parse_args()
    main(args)
    
#Program for paper:\
#Input:\
#Matrix\
#Gene\
#DB\
#Cutoff\
#\
#Output:\
#List | dist | histogram\
