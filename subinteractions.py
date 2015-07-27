#This script searches for the equal target zones of two diferent intervals
#i.e. it first uses a big matrix, and selects the interactions of the interval
#given by the user, to the rest. It does it for the two intervals given by
#the user
#Then it crosses the interactions of the two intervals and retrive the
#interactions that have targets common to the two intervals7
#it can be more specific, if you set the script to genes mode.
#in the genes mode, the script restricts even more the comon targets,
#to the targets that are genes, or contain genes. The genes can be 
#defined in a file as:
#hs(cromossome start end other_camps)

#if you want all the targets in common, use:
#python subinteractions.py all matrix chr1 start1 end1 chr2 start2 end2 outfile

#if you want the targets in genes, use:
##python subinteractions.py genes matrix chr1 start1 end1 chr2 start2 end2 genes_label_file outfile

from sys import argv


colors={"1":"thickness=1p,color=lgrey\n","2":"thickness=2p,color=green\n","3":"thickness=3p,color=blue\n","4":"thickness=4p,color=black\n","5":"thickness=5p,color=red\n"}


def make_table(infile, ch, start, end):
	"Returns a list with the interactions of the interval defined by the user"
	l=[]
	f=open(infile)
	for i in f:
		line=i.split("\t")
		intt=int(line[2])-int(line[1])
		if line[0]==ch and line[3]==ch:
			if line[1]>=start and line[2]<=end:
				if line[4]!=line[1] and int(line[1])+intt!=int(line[4]) and int(line[4])+intt!=int(line[1]):
					l.append(line[0]+" "+line[1]+" "+line[2]+" hs"+line[3]+" "+line[4]+" "+line[5]+"\t"+line[6].strip())
			elif line[4]>=start and line[5]<=end:
				if line[4]!=line[1] and int(line[4])+intt!=int(line[1]) and int(line[1])+intt!=int(line[4]):
					l.append(line[3]+" "+line[4]+" "+line[5]+" hs"+line[0]+" "+line[1]+" "+line[2]+"\t"+line[-1].strip())
	f.close()
	return l


def read_table_genes(infile):
	"""Reads the gene files, and returns a dictionary as
	dic[start coordinate]=end cordinate"""
	f=open(infile)
	cordinates={}
	for i in f:
		cords=i.split(" ")
		cordinates[int(cords[1])]=int(cords[2])
	f.close()
	return cordinates


def cross_intervals_genes (interval1, interval2, cords, outfile):
	"""crosses the two interval interactions, and the gene dictionary, 
	and writes the requisited interactions already formated for circos, to a file,.
	ATENTION: the values to add in the circos output, can be easely changed in the 
	colors dictionary"""
	out=open(outfile, "w")
	for el in interval1:
		cross1=el.split(" ")
		for le in interval2:
			cross2=le.split(" ")
			if cross1[4]==cross2[4]:
				for key, value in cords.items():
					if key>int(cross1[4]) and value<int(cross1[5].split("\t")[0]):
						if int(el.split("\t")[-1]) > 5:
							out.write("hs"+el.split("\t")[0]+" thickness=5p,color=red\n")
						if int(el.split("\t")[-1]) <= 5:
							out.write(("hs"+el.split("\t")[0]+" "+colors[el.split("\t")[-1]]))
						if int(le.split("\t")[-1]) > 5:
							out.write("hs"+le.split("\t")[0]+" thickness=5p,color=red\n")
						if int(le.split("\t")[-1]) <= 5:
							out.write(("hs"+le.split("\t")[0]+" "+colors[le.split("\t")[-1]]))							
					elif key<int(cross1[4]) and value>int(cross1[5].split("\t")[0]):
						if int(el.split("\t")[-1]) > 5:
							out.write("hs"+el.split("\t")[0]+" thickness=5p,color=red\n")
						if int(el.split("\t")[-1]) <= 5:
							out.write(("hs"+el.split("\t")[0]+" "+colors[el.split("\t")[-1]]))
						if int(le.split("\t")[-1]) > 5:
							out.write("hs"+le.split("\t")[0]+" thickness=5p,color=red\n")
						if int(le.split("\t")[-1]) <= 5:
							out.write(("hs"+le.split("\t")[0]+" "+colors[le.split("\t")[-1]]))
	out.close()

def cross_intervals(interval1, interval2, outfile):
	"""cross the interactions of the intervals and writes two a file, the 
	intercrossed interactions, already formated for circos"""
	out=open(outfile, "w")
	for el in interval1:
		cross1=el.split(" ")[4]
		for le in interval2:
			cross2=le.split(" ")[4]
			if cross1==cross2:
				if int(el.split("\t")[-1]) > 5:
					out.write("hs"+el.split("\t")[0]+" thickness=5p,color=red\n")
				if int(el.split("\t")[-1]) <= 5:
					out.write(("hs"+el.split("\t")[0]+" "+colors[el.split("\t")[-1]]))
				if int(le.split("\t")[-1]) > 5:
					out.write("hs"+le.split("\t")[0]+" thickness=5p,color=red\n")
				if int(le.split("\t")[-1]) <= 5:
					out.write(("hs"+le.split("\t")[0]+" "+colors[le.split("\t")[-1]]))
	out.close()


##############################
#Runs the program, in the desired form:
if argv[1]=="genes":
	cross_intervals_genes(make_table(argv[2], argv[3], argv[4], argv[5]), make_table(argv[2], argv[6], argv[7], argv[8]), read_table_genes(argv[9]), argv[10])
elif argv[1]=="all":
	cross_intervals(make_table(argv[2], argv[3], argv[4], argv[5]), make_table(argv[2], argv[6], argv[7], argv[8]), argv[9])
