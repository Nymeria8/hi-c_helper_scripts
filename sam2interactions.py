#This script uses the sam files from alignments of 3 diferent runs to construct a interaction matrix
#Except the intermidiate steps made by bash command line and samtools
#The data used to test was the data from the paper from dixon et al.
#The process should be made as:
#1 - Submit the sam file to filter the interactions outside the region of interest defined by the user
#that can variate from a small interval to a entire cromossome.
#2 - Writes the sam filtrated file
#the two steps above are made for the 3 runs simultaniosly
#3 - Merge the 3 runs to a single file plus the header using bash comand line: cat header ficheiro1 ficheiro2 ficheiro3 > final.sam
#4 - Convert sam to bam using samtools : samtools view -bS -o final.bam final.sam
#5 - sort: samtools sort final.bam final.sorted
#6 - remove duplicates: samtools rmdup final.sorted.bam final.rmdup.sorted.bam
#7 - sort by reads: samtools sort -n final.rmdup.sorted.bam final.rmdup.readsorted
#8 - convert bam to sam: samtools view -h -o final.rmdup.readsorted.sam final.rmdup.readsorted.bam
#9 - Parse the sam file to 2 diferent dictionaries, one with the interactions within the same cromossome, and other with the interactions between cromossomes.
#10 - Bins the interactions by intervals given by the user
#11 - Merges the informations from the interactions and the bins and writes the files, one with intra and other with inter cromossomal interactions
#12 - (optional) merge the intra and inter files to have all the interactions for that interval togheter: cat intra.matrix inter.matrix > all_interactions.matrix

#Since the 3 - 8 steps are made by bash comand line or samtools, this script can be runed in two diferent parts:

#The first part runs the step 1 and 2 described above
#The second part runs the steps 9 to 11
#To inform the script witch part of the script it shoud run, we give it the first argument as "1"
#if we want it to run the first part ou "2" if we want it to run the second part.
#Caution: the arguments acepted by the scripts are diferent when we run part one or part two of the script.

#PART 1:
#python sam2interactions.py 1 [chromossome of interst] [start of the interval] [end of the interval] [samfile1] [samfile2] [samfile3] [filtred1.sam] [filtered2.sam] [filtered3.sam]
#PART 2:
#python sam2interactions.py 2 [merged.sam] [chromossome of interest] [size of the windows or interval] [intracromossome_interactions.matrix] [intercromossome_interactions.matrix]


from sys import argv
import numpy as np
from collections import Counter

###########Global variables##########
#cromossome sizes by hg38
chrms={"chr1": 	248956422, "chr2":242193529, "chr3":198295559, "chr4":190214555, "chr5":181538259, "chr6":170805979, "chr7":159345973, "chr8":145138636, "chr9":138394717, "chr10":133797422, "chr11":135086622, "chr12":133275309, "chr13":114364328, "chr14":	107043718, "chr15":	101991189, "chr16":	90338345, "chr17":83257441, "chr18":80373285, "chr19":58617616, "chr20": 64444167, "chr21": 46709983, "chr22":50818468}
#####################################

#1
def samfilter(samfile, chrr, start, end):
	"""Recives a samfile, the cromossome name and the limits of an interval
	and filter the sam by that interval. Returns a list with the lines filtered"""
	sam=open(samfile)
	l=[]
	for i in sam:
		if i.startswith("@")==False:
			if i.split("\t")[2]==chrr or i.split("\t")[6]==chrr:
				if  int(i.split("\t")[3])> int(start) and int(i.split("\t")[3])< int(end):
					l.append(i)
				elif int(i.split("\t")[7])> int(start) and int(i.split("\t")[7])< int(end):
					l.append(i)
	sam.close()
	return l
	

#2
def samwriter(l, outfile):
	"""writes the sam list to a file."""
	out=open(outfile, "w")
	for i in l:
		out.write(i)
	out.close()

#3
#4
#5
#6
#7 
#8 
#9
def dividematrix(samfile, chrr):#divide o sam por intra relacoes e inter relacoes
	"""Slipts the samfile by inter and intra cromossome interactions. 
	Returns two dictionaries with the inter and intra cromossome interactions, as:
	intracromossome - intra[coordinate1 in interest cromossome]=[coordinate2, coordinate3, coordinate4]
	intercromossome - inter[other cromossome that not the interest one]=
	=[[coordinate1 in interest cromossome, coordinate1 in that not the interest one], [coordinate2 in interest cromossome, coordinate2 in that not the interest one]...] """
	sam=open(samfile)
	intrachm={}
	interchm={}
	temp=""
	for i in sam:
		if i.startswith("@")==False:
			if i.split("\t")[6]=="=" and i.split("\t")[0]!=temp:
				intrachm[i.split("\t")[3]]=i.split("\t")[7]
				temp=i.split("\t")[0]
			elif i.split("\t")[6]!="=" and i.split("\t")[6] in chrms and i.split("\t")[0]!=temp:
				if i.split("\t")[2]== chrr:
					if i.split("\t")[6] not in interchm:
						interchm[i.split("\t")[6]]=[[i.split("\t")[3],i.split("\t")[7]]]
						temp=i.split("\t")[0]
					else:
						interchm[i.split("\t")[6]].append([i.split("\t")[3],i.split("\t")[7]])
						temp=i.split("\t")[0]
				else:
					if i.split("\t")[2] not in interchm:
						interchm[i.split("\t")[2]]=[[i.split("\t")[7],i.split("\t")[3]]]
						temp=i.split("\t")[0]
					else:
						interchm[i.split("\t")[2]].append([i.split("\t")[7],i.split("\t")[3]])
						temp=i.split("\t")[0]					
	sam.close()
	return intrachm, interchm

#10
def bins (chrr, interval):
	"""makes the bins taking into acount the sizes of the cromossomes
	Reuturns a list of the average values in each bin
	Ex: We have a interval of 1-20, and we want bins of 5
	The list in this case would be:
	l =[2.5, 7.5, 12.5, 17.5]
	This list is used by numpy for bining"""
	windows=[]
	incr=0
	end=chrms[chrr]
	while (incr<end):
		if end-incr < int(interval):
			windows.append((incr+end)/2)
		else:
			windows.append((incr+incr+int(interval))/2)
		incr+=int(interval)
	return windows

def intra(intrachdic, chrr, interval):
	"""Uses the dictionary from make matrix, and the method bins to bin the interactions
	It returns a dictionary as 
	intra[bin of the interaction]=[correspondent bin1 of the interaction, correspondent bin2 of the interaction.... ]"""
	binss=np.array(bins(chrr,interval))
	dic_intervals={}
	for key, value in intrachdic.items():
		x=np.array([int(key),int(value)])
		inds=np.digitize(x,binss)
		if inds[0] not in dic_intervals:
			dic_intervals[inds[0]]=[inds[1]]
		else:
			dic_intervals[inds[0]].append(inds[1])
	return dic_intervals


def inter(interchdic, chrr, interval):
	"""Uses the dictionary from make matrix, and the method bins to bin the interactions
	It returns a dictionary as 
	inter[other cromossome that not the interest one]=
	=[[bin1 in interest cromossome, bin1 in that not the interest one], [bin2 in interest cromossome, bin2 in that not the interest one]...] """
	dic_intervals={}
	binss_prin=np.array(bins(chrr,interval))
	for key in interchdic:
		if key in chrms:
			binss=np.array(bins(key,interval))
			for el in interchdic[key]:
				x=np.array([int(el[0])])
				inds=np.digitize(x,binss_prin)
				x=np.array([int(el[1])])
				inds2=np.digitize(x,binss)
				if key not in dic_intervals:
					dic_intervals[key]=[[inds[0],inds2[0]]]
				else:
					dic_intervals[key].append([inds[0],inds2[0]])
	return dic_intervals


#11
def write_intra(dicintra, chrr, interval, outfile):
	"""Uses the intra output and writes the output matrix as:
	chrX	start	stop	chrX	start	stop	number of interactions"""
	binss=np.array(bins(chrr,int(interval)))
	out=open(outfile, "w")
	for key, value in dicintra.items():
		refs=Counter(value)
		origin= chrr.split("r")[-1]+"\t"+str(binss[key]-int(interval)/2)+"\t"+str(binss[key]+(int(interval)/2))
		for key2, value2 in refs.items():
			target= chrr.split("r")[-1]+"\t"+str(binss[key2]-(int(interval)/2))+"\t"+str(binss[key2]+(int(interval)/2))
			out.write(origin+"\t"+target+"\t"+str(value2)+"\n")
	out.close()


def write_inter(dicinter, chrr, interval, outfile):
	"""Uses the inter output and writes the output matrix as:
	chrY	start	stop	chrX	start	stop	number of interactions"""
	out=open(outfile, "w")
	binss_prin=np.array(bins(chrr,int(interval)))
	for key, value in dicinter.items():
		interdic={}
		binss=np.array(bins(key,int(interval)))
		for el in value:
			if el[0] not in interdic:
				interdic[el[0]]=[el[-1]]
			else:
				interdic[el[0]].append(el[-1])
		for key2, value2 in interdic.items():
			refs=Counter(value2)
			origin= chrr.split("r")[-1]+"\t"+str(binss_prin[key2]-(int(interval)/2))+"\t"+str(binss_prin[key2]+(int(interval)/2))
			for key3,value3 in refs.items():
				target= key.split("r")[-1]+"\t"+str(binss[key3]-(int(interval)/2))+"\t"+str(binss[key3]+(int(interval)/2))
				out.write(origin+"\t"+target+"\t"+str(value3)+"\n")
	out.close()
			
########################################################## RUNS THE ABOVE METHODS BY PARTS###############################################################################################
def define (step):
	"""Runs the script by the parts determinated by the user"""
	if step =="1":
		samwriter(samfilter(argv[5], argv[2], argv[3], argv[4]), argv[8])
		samwriter(samfilter(argv[6], argv[2], argv[3], argv[4]), argv[9])
		samwriter(samfilter(argv[7], argv[2], argv[3], argv[4]), argv[10])
	elif step == "2":
		intrac, interc =dividematrix(argv[2], argv[3])
		write_intra(intra(intrac, argv[3], argv[4]), argv[3], argv[4], argv[5])
		write_inter(inter(interc, argv[3], argv[4]), argv[3], argv[4], argv[6])

#CALLS THE SCRIPT
define(argv[1])
		
