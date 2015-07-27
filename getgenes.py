#This script uses the biomart server to search for genes in the interval
#of the genome defined by the user.
#The specie/database used is easly edited by the user
#After the search the method writer will use the search information
#to make 3 inputs to circos
#one with the exon/intron structure of the genes color coded by type of gene
#other 2 with the labels, to be placed in diferent directions
#Usage:
#python getgenes.py chromossome start_interval end_interval out_exon/intro out_label1 out_label2

from sys import argv
from biomart import BiomartServer


#Setup of the server
server = BiomartServer( "http://www.biomart.org/biomart" )

#avaiable databases
#print server.show_databases()
#avaiable datasets
#print server.show_datasets()

#setup the dataset
hg=server.datasets["hsapiens_gene_ensembl"]

#avaiable filters/atributes
#print hg.show_filters()
#print hg.show_attributes()

#biotypes and colors avaiable
color={'protein_coding':'chr5', 'lincRNA': 'chr8', 'miRNA':'chr13', 'processed_transcript':'chr14', 'snRNA':'chr1', 'snoRNA':'chr1'}

def search(chrr, start, end):
	"""makes the search against the biomart database, using as filter
	the cromossome, start and end coordinates, and returns lists of 
	gene byotype, exon start/end and the name of the gene. If the gene
	as several exons, it will retrive a list for each exon.
	Then, it organizes the information in 2 dictionaries:
	one with the coordinates: dic[genename]=[start1, end1, start2, end2...]
	and another with lab[genename]=biotype of gene
	returns the two dictionaries"""
	dic={}
	lab={}
	response = hg.search({
	  'filters': {
		 'chromosome_name': [chrr], 'start': [start], 'end':[end]
	  },
	  'attributes': ['gene_biotype',  'exon_chrom_start', 'exon_chrom_end', 'external_gene_name'
		  
	  ]
	})
	print("search done!")
	for line in response.iter_lines():
		d=line.split("\t")
		if d[3].split(".")[0] not in dic:
			dic[d[3].split(".")[0]]=[d[1], d[2]]
			lab[d[3].split(".")[0]]=d[0]
		else:
			dic[d[3].split(".")[0]].append(d[1])
			dic[d[3].split(".")[0]].append(d[2])
	return dic, lab

def writer (chrr, dictionary, lab, outfile_coords, outfile_label1, outfile_label2):
	"""organizes the two dictionaries from search, into usable files by circos
	divides the labels by protein_coding and non protein_coding. This will
	be color coded by the types of genes represented. if the biotype not in color, 
	it will not appear on the outputs"""
	outcord=open(outfile_coords, "w")
	outlab1=open(outfile_label1, "w")
	outlab2=open(outfile_label2, "w")
	for key, value in dictionary.items():
		if key!=" ":
			if lab[key] in color:
				if lab[key]=='protein_coding':
					outlab1.write("hs"+str(chrr)+" "+value[0]+" "+value[-1]+" "+key+"\n")
				else:
					outlab2.write("hs"+str(chrr)+" "+value[0]+" "+value[-1]+" "+key+" color="+color[lab[key]]+"\n")
				a=1
				while(a<len(value)):
					if a%2==1:
						outcord.write("hs"+str(chrr)+" "+value[a-1]+" "+value[a]+" fill_color="+color[lab[key]]+",r0=0.97r,r1=0.97r+50p\n")
						a+=1
					else:
						outcord.write("hs"+str(chrr)+" "+value[a-1]+" "+value[a]+" fill_color="+color[lab[key]]+",r0=0.985r,r1=0.985r+5p\n")
						a+=1
	outcord.close()
	outlab1.close()
	outlab2.close()
	
dic,lab =search(argv[1], argv[2], argv[3])
writer(argv[1], dic,lab, argv[4], argv[5], argv[6])
