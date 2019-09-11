#!/usr/bin/env python
import pandas as pd;
import numpy as np;
import sys;
if __name__ == '__main__':
    Query1 = sys.argv[1]
    pfam_out = Query1.replace('.tsv','').replace('.txt','').replace('.csv','')
    pf = open("{}_formatted.tsv".format(pfam_out),"w+")
    ColTitles = "target_name    pfam_accession   tlen    query_name  accession   qlen    E-value     score   bias   #  of      c-Evalue    i-Evalue    score   bias    hmm_from    hmm_to  ali_from    ali_to  env_from    env_to  acc description_of_target"
    ColTitles=ColTitles.split()
    
    ColTitles=str(ColTitles)
    ColTitles=ColTitles.strip('[').strip(']').replace(',', '\t').replace("'","")

    pf.write(ColTitles + '\n')
    with open("{}".format(Query1),'r') as fp:
        for i in range(2):
            next(fp)
        for line in fp:

            if '#' in line:
                line = ''

            line = str(line)
            line = (line.split())
            line = str(line)
            line=line.strip('[').strip(']').replace(',','\t').replace("'","")
            if len(line) != 0:
            	pf.write(line + '\n')
    pf.close()
