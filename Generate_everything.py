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
    def __init__(self):
        return
    def calcDist(self,filename,matrixfile,single_cutoff=""):
        count=0
        matrix=self.read_Matrix(matrixfile)
        self.seqlen=len(list(matrix.values())[0])
        k=[k for k in matrix.keys()][0]
        self.seqsize=len(matrix[k[0]])
        self.seqs=dict()
        self.cutoff_point=single_cutoff
        self.distance_cutoff=dict()
        if filename[-2:]=="gz":
            distfile=filename.replace(".gz",".dist")
            self.processGZfile(filename,matrix,distfile)
        else:
            distfile=filename+".dist"
            self.processTXTfile(filename,matrix,distfile)
        return distfile
    def add_templates(self,seqTemplates):
        self.templates=seqTemplates

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
            dist=self.seqlen
            for i,nt in enumerate(seq):
                if float(matrix[nt][i])>float(cutoff):
                    dist-=1
            all_dist.append(dist)
        distance="\t".join(["{}".format(i) for i in all_dist])
        dist= self.seqlen
        score=0
        for i,nt in enumerate(seq):
            score+=float(matrix[nt][i])
            if float(matrix[nt][i])>float(self.cutoff_point):
                dist-=1
        if dist in self.distance_cutoff:
            self.distance_cutoff[dist]+=1
        else:
            self.distance_cutoff[dist]=1
            
        if dist<=numpy.min([9,self.seqlen/2]):
            if dist in self.seqs:
                self.seqs[dist].append((seq,score))
            else:
                self.seqs[dist]=[(seq,score)]
            
        return distance
    def print(self,filename,maxNumber=100):
        fn,f_ext=os.path.splitext(filename)
        outseqs=[]
        number=0
        for key in self.seqs.keys():
            for seq in self.seqs[key]:
                if number<maxNumber:
                    outseqs.append((seq[0],seq[1]))
                    number+=1
            outseqs.sort(key=lambda y: y[1])
        with open(fn+f_ext,"w") as ofile:
            for i,seq in enumerate(outseqs):
                ofile.write(">seq_{} score={:0.4f}\n{}\n".format(i+1,seq[1],seq[0]))                    
        return True
    def printN(self,filenamebase):
        fn,f_ext=os.path.splitext(filenamebase)
        for key in self.seqs.keys():
            outseqs=[]
            if key>4:
                continue
            for seq in self.seqs[key]:
                outseqs.append((seq[0],seq[1]))
            outseqs.sort(key=lambda y: y[1])
            if len(outseqs)<1:
                continue
            with open(fn+".N{}".format(key)+".fasta","w") as ofile:
                for seq in outseqs:
                    ofile.write(">score={:0.4f}\n{}\n".format(seq[1],seq[0]))                    
        return True
    def processGZfile(self,filename,matrix,outname):
        with open(outname,'w') as out:
            with gzip.open(filename, 'rb') as f:
                while True:
                    line=f.readline()
                    if not line:
                        break
                    tmpseq=line.rstrip().split()[0]
                    seq=tmpseq.decode("utf-8")
                    if len(seq)==self.seqlen:
                        matchTemp=True
                        if "templates" in dir(self):
                            for template in self.templates:
                                matchTemp=self.matchTemplate(seq,template)
                                if matchTemp:
                                    break
                        if not matchTemp:
                            continue
                        score=self.score_Seq(seq,matrix)
                        out.write("{}\t{}\n".format(seq,score))

    def processTXTfile(self,filename,matrix,outname):
        with open(outname,'w') as out:
            with open(filename, "r") as f:
                while True:
                    line=f.readline()
                    if not line:
                        break
                    seq=line.rstrip()
                    if len(seq)==self.seqlen:
                        matchTemp=True
                        if "templates" in dir(self):
                            for template in self.templates:
                                matchTemp=self.matchTemplate(seq,template)
                                if matchTemp:
                                    break
                        if not matchTemp:
                            continue
                        score=self.score_Seq(seq,matrix)
                        out.write("{}\t{}\n".format(seq,score))
    def matchTemplate(self,seq,template):
        for i,AA in enumerate(template):
            if "X"==AA:
                continue
            elif seq[i]==AA:
                continue
            else:
                return False
        return True

class PlotHist():
    def __init__(self):
        self.threshold_list=[x/5.0 for x in range(-15,1)]
        self.maxdist=20
    def run(self,filename,outname,scaled_value):
        distanceHash=dict()
        for T in self.threshold_list:
            distanceHash[T]=numpy.zeros([1,self.maxdist+1])
        with open(filename, "r") as infile:
            while True:
                line=infile.readline()
                if not line:
                    break
                dataline=line.rstrip().split()
                matchTemp=True
                if "templates" in dir(self):
                    seq=dataline[0]
                    for template in self.templates:
                        matchTemp=self.matchTemplate(seq,template)
                        if matchTemp:
                            break
                if not matchTemp:
                    continue
                for i,value in enumerate(dataline[1:]):
                    distanceHash[self.threshold_list[i]][0,int(value)]+=1
                    
        heatmap=numpy.zeros([len(self.threshold_list),self.maxdist+1])
        scaledheatmap=numpy.zeros([len(self.threshold_list),self.maxdist+1])
        for i,T in enumerate(self.threshold_list):
            fullCount=numpy.sum(distanceHash[T])
            for pos,n in enumerate(distanceHash[T][0,:]):
                heatmap[i,pos]=n/float(fullCount)
                scaledheatmap[i,pos]=n/float(scaled_value)
        self.plot_heatmap(outname,heatmap)
        writeTable(outname.replace(".pdf",".table.txt"),heatmap, self.threshold_list)
        self.plot_cum_heatmap(outname.replace(".pdf",".cumulative.pdf"),heatmap)

        self.plot_heatmap(outname.replace(".pdf",".scaled.pdf"),scaledheatmap)
        writeTable(outname.replace(".pdf",".scaled.table.txt"),scaledheatmap, self.threshold_list)   
        self.plot_cum_heatmap(outname.replace(".pdf",".scaled.cumulative.pdf"),scaledheatmap)

        self.plot_legend(outname.replace(".pdf",".legend.pdf"))

    def add_templates(self,seqTemplates):
        self.templates=seqTemplates
    def GenerateBins(self,data_list):
        bins=numpy.zeros([1,self.maxdist+1])
        data=numpy.array(data_list)
        for i in range(0,self.maxdist+1):
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
        ax.set_xticks(range(0,self.maxdist+1))
        ax.grid(visible=None)
        plt.grid(visible=None)
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
        ax.set_xticklabels(["$\leq${}".format(x) for x in range(0,self.maxdist+1)])
        ax.grid(visible=None)
        plt.grid(visible=None)
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
        
    def matchTemplate(self,seq,template):
        #print(template)
        for i,AA in enumerate(template):
            if "X"==AA:
                continue
            elif seq[i]==AA:
                continue
            else:
                return False
        return True

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
    all_args.add_argument("-m","--matrix",default="",required=False,help="matrix data")
    all_args.add_argument("-s","--sequences",default="",required=False,help="set of sequences")
    all_args.add_argument("-d","--distancefile",required=False,help="Distance file, generated from program or independently",default="")
        
    all_args.add_argument("-t","--seqtemplate",action="extend",nargs="+",type=str,required=False,help="gene")
    all_args.add_argument("-c","--cutoff",required=False,default=-0.2,type=float,help="cut off for a hit, (default:-0.2)")
    all_args.add_argument("-o","--pdfname",required=False,default=[],help="name of the final pdf file")
    all_args.add_argument("--scale",required=False,default=85149053,type=int,help="value to scale the data by")
    all_args.add_argument("--writeN",dest='writetop',required=False,action='store_true',help="set flag to write out with less then 4 distance seqs")
    all_args.add_argument("--writeout",nargs='?',required=False,default=-1,type=int,help="write out top N sequences")
    all_args.add_argument("--noseqout",dest='writetop',required=False,action='store_false',help="set flag to NOT write out less then 4 distance seqs")
    all_args.set_defaults(writetop=False)
    #scaled_value=85149053
    #need arg for no pdf printout
    args = vars(all_args.parse_args())
    return args

def main(args):
    
    matrixfile=args["matrix"]
    sequencefile=args["sequences"]

    SDC=SeqDistCalc()
    plothist=PlotHist()
    
    if args["seqtemplate"]:
        SDC.add_templates(args["seqtemplate"])
    if not args["writeout"]:
        args["writeout"]=100

    #need to deal with nothing getting through
    distfile=""
    if os.path.exists(args["distancefile"]):
        distfile=args["distancefile"]
        if args["seqtemplate"]:
            plothist.add_templates(args["seqtemplate"])
    elif os.path.exists(matrixfile) and os.path.exists(sequencefile):
        distfile=SDC.calcDist(sequencefile,matrixfile,args["cutoff"])

        if args["writetop"]:
            SDC.printN(sequencefile)
        if args["writeout"]>0:
            SDC.print(sequencefile+".fasta",args["writeout"])
        for dist in range(0,SDC.seqsize+1):
            if dist in SDC.distance_cutoff:
                print("{}\t{}".format(dist,SDC.distance_cutoff[dist]))
            else:
                print("{}\t{}".format(dist,0))

        if len(SDC.seqs.keys())<=0:
            print("No sequences pass templates or size")
            exit()

    if args["pdfname"]:
        fn,f_ext=os.path.splitext(args["pdfname"])
        if not f_ext:
            fn=fn+".pdf"
        elif ".pdf" not in f_ext:
            fn=args["pdfname"]+".pdf"
        else:
            fn=fn+f_ext
        plothist.run(distfile,fn,args["scale"])
    elif os.path.exists(distfile):
        plothist.run(distfile,distfile.replace(".txt","").replace(".dist",".pdf"),args["scale"])
    else:
        print("no distance file given or created")

if __name__=="__main__":
    args=parse_args()
    main(args)
