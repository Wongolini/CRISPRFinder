#!/usr/bin/env python
import sys;
if __name__ == '__main__':
	Query1 = sys.argv[1];
	with open (Query1,"r") as mefile:
		accessions = mefile.read()
	accessions = str(accessions)
	num_commas = accessions.count(',')
	gstring = accessions
		
	Acc_write = open("recollect_aa_homologs.sh","w+")
	Acc_write_nt = open("recollect_nt_homologs.sh","w+")
	div = 500

	if num_commas > 999:
            while(num_commas%div!=0):
            	div-=1;
            num_commas=num_commas/div;
            gstring = gstring.split(',')
            glist = []

            for k in range(1,num_commas):
                glist.append(gstring[(k-1)*div:k*div])
            for kk in range(1,len(glist)):
                good = str(glist[kk]).replace("'","").replace('[','').replace(']','').replace(' ','')
                Acc_write.write('efetch -db sequences -format fasta_cds_aa -id '+ str(good)+'\n')
		Acc_write_nt.write('efetch -db sequences -format fasta_cds_na -id ' +str(good) + '\n')

	Acc_write.close()
	Acc_write_nt.close()
	
		
