import numpy as np
import sys
import dihmm_ext
import bed_writer
import os
from os import listdir,mkdir
from os.path import isfile,join,exists
from optparse import OptionParser
import pandas as pd

def printlist(path,value,sep):
    # "print the list to text file"
    fid=open(path,'w')
    nline=[]
    for listvalue in value:
        if isinstance(listvalue[0],float):
            nline=sep.join(str(i) for i in listvalue)
        elif isinstance(listvalue[0],int):
            nline=sep.join(str(i) for i in listvalue)
        elif isinstance(listvalue[0],str):
            nline=sep.join(listvlaue)
        fid.write("%s\n" % nline)
    fid.close()

def main(argv):
    parser = OptionParser()
    parser.add_option("-t", "--table", action="store", type="string", dest="table", metavar="<file>", help="the file for currently annotated samples")
    parser.add_option("-n", "--name", action="store", type="string", dest="name", metavar="<file>", help="the name of the annotation pair to work on")
    parser.add_option("-i", "--indata", action="store", type="string", dest="indata", metavar="<file>", help="the path for prepared signal matrix in the target cell type")
    parser.add_option("-m", "--model", action="store", type="string", dest="model", metavar="<file>", help="the path for trained model")
    parser.add_option("-o", "--outfolder", action="store", type="string", dest="outfolder", metavar="<file>", help="the path for output folder to store the annotated chroms")
    
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 8:
        parser.print_help()
        sys.exit(1)
    
    domain_size = 20
    path1=opt.indata
    file1=[path1+f for f in listdir(path1)]
    allfile=file1

    cdf = pd.read_csv(opt.table, sep='\t', header=None)
    currentAnnotations = list(cdf[0])
    name = opt.name
    if name not in currentAnnotations:
        print("working on "+name+"...")
        os.mkdir(opt.outfolder)
        x=dihmm_ext.load_model(opt.model+'/',domain_size)
        p=x.emission_probabilities
        bt=x.bin_transition_probabilities
        dt=x.domain_transition_probabilities
        nb=x.n_bin_states
        nd=x.n_domain_states

        output=opt.outfolder+'/'
        for file in allfile:
            tmp=file.split('/')[-1].split('_')
            cellline=tmp[0]

            chrom=tmp[1]
            noutput=output+cellline+'/'
            if not os.path.exists(noutput):
                os.mkdir(noutput)
            a=dihmm_ext.annotate(x,[file])
            b=bed_writer.BedWriter(a[0],x)
            b.write_bed_files(noutput,cellline,chrom)
            anno=a[0].annotations
            bsd=a[0].bin_state_distributions
            dsd=a[0].domain_state_distributions
            printlist(noutput+chrom+'_bin_domain_states.txt',anno.tolist(),'\t')
    else:
        print(opt.name+" has already exsited...skip...")
        
if __name__ == "__main__":
    main(sys.argv)
