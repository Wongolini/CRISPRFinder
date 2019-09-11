#!/usr/bin/env python
#import matplotlib.pyplot as plt
import sys;
import numpy as np
import os
from Bio import SeqIO;

if __name__ == '__main__':
    accessions = []
    Query1 = sys.argv[1]
 
    length_list = []
    for sec_lengths in SeqIO.parse(Query1, "fasta"):
        longboi = sec_lengths.seq
        length_list.append(len(longboi))


    accessions=[]
    for seq_rec in SeqIO.parse(Query1, "fasta"):
        accessions.append(len(seq_rec.seq))

    proteinseq = Query1.replace('.fasta', '').join('_proteins.faa')
    #faa = open("filtered_{}".format(proteinseq),"w+")	
    f = open("filtered_{}".format(Query1),"w+")
    filename = Query1.replace('.fasta' , '').replace('pfam_accessions', '')
    g = open("recollect_aa_homologs_{}.sh".format(filename),"w+")
    for check in SeqIO.parse(Query1, "fasta"):
        checker = check.seq
        checkid = check.id
        s = 'lcl|'
	cds = '_cds_'
        description = check.description
        if s in description:
	    #description = description.split()[1:]
	    description = description.split()
	    description = str(description).replace("'","")
	if cds in checkid:
	    checkid = checkid.replace('_cds_','@')
	    cds_index=checkid.index('@')
	    checkid = checkid[cds_index:]
	    checkid = checkid.replace('@','')
	    if '.' in checkid:
		dtrm = checkid.index('.')
            	checkid = checkid[:dtrm]

	   #have protein accessions that don't have accessions? 

        for i in range(len(checker)):
            if checker[i] not in ['A','T','C','G']:
                print(description)
		check.seq = ''
                check.id= ''
		

            else:
                pass;


        if ((len(check.id)) != 0 or len(check.seq)!=0):
            f.write('>' + str(checkid)+' '+ str(description) +'\n')
            f.write(str(checker)+'\n')
	    gstring = str(checkid)+','

            #faa.write('>' + str(checkid)+ ' ' + str(description) + '\n')
	    #faa.write(str(checker.translate())+'\n')
    	div = 500
	num_commas = gstring.count(',')
	while(num_commas%div!=0):
	    div-=1;
	num_commas=num_commas/div;
	gstring = gstring.split(',')
	glist = [[]]
	print(gstring)
	for k in range(1,num_commas):
	    glist.append(gstring[(k-1)*div:k*div])
	for kk in range(len(glist)):
	    good = str(glist[kk]).replace("'","")
	    g.write('efetch -db sequences -format fasta_cds_aa -id '+ str(good)+'\n')
	
	
            
    
    f.close()
    g.close()
    #faa.close()
