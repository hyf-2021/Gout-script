import sys
import os
import gzip
import pprint
from optparse import OptionParser
from collections import Counter
import argparse as ap

parser = ap.ArgumentParser()
parser.add_argument('-p', "--pe", dest="soap_pe_file", help="input SOAP PE file", metavar="PE")
parser.add_argument('-s', "--se", dest="soap_se_file", help="input SOAP SE file", metavar="SE")
parser.add_argument('-l', "--len", dest="gene_len_file", help="input gene length file", metavar="LEN")
parser.add_argument('-o', "--out", dest="out_file", help="output abundance result", metavar="ABUNDANCE")
options  = vars(parser.parse_args())


soap_pe_file =options['soap_pe_file']
soap_se_file =options['soap_se_file']
gene_len_file=options['gene_len_file']
out_ab_file  =options['out_file']

def parse_soap_pe_file(reads,genes,soap_file, all_ref_length):
	all_insert_size = []
	flag = 1	### flag = 1 means reading read1 now;
	f_in=gzip.open(soap_file,'rb')
	for line in f_in:
		arr=line.decode().strip().split('\t')
		read_id = arr[0]
		gene_id = arr[7]
		map_len = int(arr[5])
		site = int(arr[8])
		if arr[3] == "1":
			if flag:
				site1 = site
				L1 = map_len
				flag = 0
			else:
				site2 = site
				L2 = map_len
				insert_size = site2-site1+L2 if arr[6]=="-" else site1-site2+L1
				all_insert_size.append(insert_size)
				flag = 1

		if read_id not in reads:
			reads[read_id]={}
		reads[read_id][gene_id]=None #Initialize Co

		if gene_id not in genes:
			ref_length = all_ref_length[gene_id]
			genes[gene_id]=[None,None,None,None,None,None,{},None] #[GeneLength,U,M,Ab(U),Ab(M),Ab(G),{read_id1:1,read_id2:1},gene_num]
			genes[gene_id][0] = ref_length[0]
			genes[gene_id][7] = ref_length[1]
		if read_id not in genes[gene_id][6]:
			genes[gene_id][6][read_id]=1
		else:
			genes[gene_id][6][read_id]+=1

	f_in.close()
	all_insert_size = dict(Counter(all_insert_size))
	most_size = sorted(all_insert_size.items(), key=lambda all_insert_size:all_insert_size[1], reverse=True)[0]
	print("insertsize:"+str(most_size[0])+"\n" +"abundance:"+str(most_size[1])+"\n")
	return most_size[0]


def parse_soap_se_file(reads, genes, soap_file, insert_size, all_ref_length):
	f_in=gzip.open(soap_file,'rb')
	for line in f_in:
		arr=line.decode().strip().split('\t')
		read_id = arr[0] + "_" + arr[4]
		gene_id = arr[7]
		map_len = int(arr[5])
		site = int(arr[8])
		ref_length = all_ref_length[gene_id]
		if arr[6] == "+":
			if site + insert_size < ref_length[0]:continue
		else:
			if site + map_len > insert_size:continue

		if read_id not in reads:
			reads[read_id]={}
		reads[read_id][gene_id]=None #Initialize Co

		if gene_id not in genes:
			genes[gene_id]=[None,None,None,None,None,None,{},None] #[GeneLength,U,M,Ab(U),Ab(M),Ab(G),{read_id1:1,read_id2:1},gene_num]
			genes[gene_id][0] = ref_length[0]
			genes[gene_id][7] = ref_length[1]

		if read_id not in genes[gene_id][6]:
			genes[gene_id][6][read_id]=1
		else:
			genes[gene_id][6][read_id]+=1

	f_in.close()
	all_ref_length.clear()

def read_ref_length(gene_len_file):
	f_in=gzip.open(gene_len_file,'rb')
	tmp=f_in.readline()
	all_ref_length = {}
	for line in f_in:
		arr=line.decode().strip().split('\t')
		gene_num=arr[0]
		gene_id=arr[1]
		length = int(arr[2])
		all_ref_length[gene_id] = [length, gene_num]

	#if(genes.has_key(gene_id)):
	#	genes[gene_id][0]=length
	#		   genes[gene_id][7]=gene_num
	f_in.close()
	return all_ref_length

reads={}
genes={}
all_ref_length = read_ref_length(gene_len_file)
insert_size = parse_soap_pe_file(reads,genes,soap_pe_file, all_ref_length)
parse_soap_se_file(reads,genes,soap_se_file, insert_size, all_ref_length)

for gene_id in genes:
	reads_cluster=genes[gene_id][6]
	U=0
	M=0
	for read_id in reads_cluster:
		genes_cluster=reads[read_id]
		if len(genes_cluster)==1: #Uniq mapping
			U+=1
		elif len(genes_cluster)>1: #Multi mapping
			M+=1
		else:
			sys.exit('Mapping must be either Uniq or Multi!')
	genes[gene_id][1]=U
	genes[gene_id][2]=M

	length=genes[gene_id][0]
	if length:
		ab_U=U/float(length)
		genes[gene_id][3]=ab_U
	else:
		sys.exit('Gene length must NOT be 0 or None!')
for read_id in reads:
	genes_cluster=reads[read_id]
	gene_arr=genes_cluster.keys()
	
	len_arr =map(lambda x:genes[x][0],gene_arr)
	U_arr   =map(lambda x:genes[x][1],gene_arr)
	M_arr   =map(lambda x:genes[x][2],gene_arr)
	abU_arr =map(lambda x:genes[x][3],gene_arr)

	total_U   =sum(U_arr)
	total_len =sum(len_arr)
	total_ab_U=sum(abU_arr)

	for gene_id in gene_arr:
		ab_U=genes[gene_id][3]
		Co=ab_U/float(total_ab_U) if ab_U > 0 else 0
		reads[read_id][gene_id]=Co

total_ab_G=0.0
for gene_id in genes:
	reads_cluster=genes[gene_id][6]
	length=genes[gene_id][0]
	ab_U  =genes[gene_id][3]

	total_Co=0.0
	for read_id in reads_cluster:
		genes_cluster=reads[read_id]
		if len(genes_cluster)>1: #Multi mapping
			total_Co+=genes_cluster[gene_id]
	ab_M=total_Co/float(length)

	genes[gene_id][4]=ab_M
	ab_G=ab_U+ab_M
	genes[gene_id][5]=ab_G
	total_ab_G += ab_G

f_out=open(out_ab_file,'w')
f_out.write('#Gene_num\tGene_length\tUnique_align_reads\tMultiplt_align_reads\tUnique_align_copynumber(U)\tMultiplt_align_copynumber(M)\tTotal_copynumber(G)\tRelative_abundance(G)\n')
for gene_id in sorted(genes.keys()):
	if genes[gene_id][1] == 0:continue
	length=genes[gene_id][0]
	U	 =genes[gene_id][1]
	M	 =genes[gene_id][2]
	ab_U  =genes[gene_id][3]
	ab_M  =genes[gene_id][4]
	ab_G  =genes[gene_id][5]
	gene_num=genes[gene_id][7]
#	gene_num = gene_id
	re_ab_G = ab_G/total_ab_G

	f_out.write(gene_num+'\t'+str(length)+'\t'+str(U)+'\t'+str(M)+'\t'+str(ab_U)+'\t'+str(ab_M)+'\t'+str(ab_G)+'\t'+str(re_ab_G)+'\n')
f_out.close()

