#!/usr/bin/env python
#Need Accessions_env.tsv and final_gene_nt_sequences.fastsa and combine data to fetch domain sequences within the genes
import sys;
from Bio import SeqIO;

if __name__ == '__main__':
    Query1 = sys.argv[1]
    Query2 = sys.argv[2]
    Query3 = sys.argv[3]
    domain,regions,fast='','',''
    if Query3.endswith(".fasta"):
        fast = Query3
    if Query1.endswith(".fasta"):
        fast = Query1
    if Query2.endswith(".fasta"):
        fast = Query2
    if (Query1.endswith(".tsv") or Query1.endswith(".txt") or Query1.endswith(".csv")):
        regions = Query1
    if (Query2.endswith(".tsv") or Query2.endswith(".txt") or Query2.endswith(".csv")):
        regions = Query2
    if (Query3.endswith(".tsv") or Query2.endswith(".txt") or Query2.endswith(".csv")):
        regions = Query3
    if ('PF' in Query1):
        domain = Query1
    if ('PF' in Query2):
        domain = Query2
    if ('PF' in Query3):
        domain = Query3
    keys = ['query_name', 'E-value', 'env_to', 'env_from']
    treasures = {'query_name':[],
                'E-value':[],
                'env_to':[],
                'env_from':[],
                'pfam_accession':[]
                }

    name_domain = domain.replace('.','_')
    for_out_fast = fast.replace(".fasta","")
    output = open("{}_with_pfam_domain_{}.fasta".format(for_out_fast,name_domain), "w+")
    output_control = open("{}_control.fasta".format(for_out_fast),"w+")
    output_control_aa = open("{}_control_aa.fasta".format(for_out_fast),"w+")
    pfam_accession_index,query_accession_index,Eval_index,env_to_index,env_from_index,pfam_accession_index=0,0,0,0,0,0;
    with open("{}".format(regions)) as fp:
        for line in fp:

            line = line.split()
            query_name_index = line.index("query_name")
            Eval_index = line.index("E-value")
            env_to_index = line.index("env_to")
            env_from_index = line.index("env_from")
            pfam_accession_index = line.index("pfam_accession")
            break;
        for lines in fp:
            lines = lines.split() #only append if it sees the pfam domain!!!!!
            
	    if domain in lines:
		print(lines[query_name_index])
		if '_cds_' in lines[query_name_index]:
		    query_name = lines[query_name_index]
		    cdsm = query_name.index('_cds_')
		    query_name = query_name[cds:]
		    query_name =query_name.replace('_cds_','')
                    if '.' in query_name:
                        strife = query_name.index('.')
                        query_name = query_name[:strife]
		    lines[query_name_index] = query_name
		if '_prot_' in lines[query_name_index]:
		    query_name = lines[query_name_index]
                    protm = query_name.index('_prot_')
                    query_name = query_name[protm:]
		    query_name =query_name.replace('_prot_','')
                    if '.' in query_name:
                        strife = query_name.index('.')
                        query_name = query_name[:strife]
                    lines[query_name_index] = query_name
		#print(lines[query_name_index])
                    #print(lines[query_name_index])
                treasures["query_name"].append(lines[query_name_index].replace('lcl|',''))

                treasures["E-value"].append(lines[Eval_index]);
                treasures["env_to"].append(lines[env_to_index]);
                treasures["env_from"].append(lines[env_from_index]);
                treasures["pfam_accession"].append(lines[pfam_accession_index])



    for sec in SeqIO.parse(fast, "fasta"):
        recid = sec.id #this requires the accessions in all files be the same!
	print(recid)
        rec_seq = sec.seq
        if recid  in treasures["query_name"]:
            fuck = (treasures["query_name"].index(recid))
#           print(treasures["query_name"][fuck])
        if 'lcl' in recid:
            stripper = recid.index('_cds_')
            recid = recid[stripper:]
            recid=recid.replace('_cds_','')
            if '.' in recid:
                dot = recid.index('.')
                recid = recid[:dot]
            for unc in ['!','@','#','$','%','^','&','*','(',')',"'",'/','_','-','+','=']:
                recid = recid.replace(unc,'')
#       print(recid)

        if recid in treasures["query_name"]:
            pos = treasures["query_name"].index(recid)
            seq_end = int(treasures["env_to"][pos])
            seq_start = int(treasures["env_from"][pos])


            if seq_start >1:
                seq_start = seq_start*3-2-1
                if seq_start > 20:
                    seq_start-=15
                while (seq_start%3 != 0):
                    seq_start-=1

            else:
                seq_start=seq_start*3-3
                while(seq_start%3!=0):
                    seq_start-=1
            if seq_end*3+15 <= len(rec_seq):
                seq_end = seq_end*3+15
                while((seq_end-seq_start)%3!=0):
                    if (seq_end<=len(rec_seq)):
                        seq_end+=1
                    else:
                        seq_end-=5

            else:
                seq_end = seq_end*3-1
                while(((seq_end-seq_start)%3!=0) and (seq_end%3!=0)):
                    if (seq_end>len(rec_seq)):
                        seq_end-=1
                    else:
                        seq_end-=1
                    #print(seq_end)
                    #print('*')


            #print((seq_end-seq_start)/3)
            #print(seq_start/3)
            #print(seq_start)
            #print(seq_end)
            #print('len ' + str(len(rec_seq)))
            #print(rec_seq[seq_start:seq_end])
            description = sec.description.split()[1:]

            description= str(description)
            description = description.replace("'","").replace(',',' ')

            if (((seq_end-seq_start)%3)==0):
                output.write(">" + str(recid) +' ' + str(description) + " region with pfam domain " + str(domain)+ " E-val: " + str(treasures["E-value"][pos]) + " sequence_length: " + str(len(rec_seq[seq_start:seq_end])) + ' range: (' + str(seq_start) +','+ str( seq_end) +')' +'\n')
                output.write(str(rec_seq[seq_start:seq_end])+'\n')
                output_control.write('>' + str(recid) + '\n')
                output_control.write(str(sec.seq) + '\n')
                output_control_aa.write('>' + str(recid)+' ' + str(sec.description) + 'pfam domain sequence translated \n')
                output_control_aa.write(str(rec_seq[seq_start:seq_end].translate())+'\n')
                #print('>' + recid+ ' aa seq with domain')
                #print(rec_seq[seq_start:seq_end].translate())
                #print('>' + recid + 'full protein')
                #print(rec_seq.translate())
                #if rec_seq[seq_start:seq_end] in rec_seq:
                    #print('true')
        #print(recid)
                    #print(sec.description)


                #else:
                 #   pass;
            else:
                print(str(recid) + 'not in')
                break;
    output.close()
    output_control.close()

