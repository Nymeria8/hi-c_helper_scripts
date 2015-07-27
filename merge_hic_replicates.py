#merge the hi-c replicates. only the entries present in the two replicates
#are present in the final set
#i.e. only the lower value of interactions is present in the final output
#only the interactions for the given chromossome are present in the output
#usage: python merge_hic_replicates.py chromossome_number replicate1.tab replicate2.tab output

from sys import argv

def read_replicate (ch, infile):
	"""parse the replicate file to a dictionary"""
	f=open(infile)
	cores={}
	for i in f:
		line=i.split("\t")
		if ch==line[0] and ch==line[3]:
			cores[line[1]+"\t"+line[2]+","+line[4]+"\t"+line[5]]=int(line[6].strip())
	f.close()
	return cores


def merge_replicates(ch, replicate1, replicate2, mergedoutfile):
	"""merge the two replicates, and write a output similar to the replicates.
	If interactions from two fragments exists in both replicates, the lower number
	of replicates is present in the output. if it is similar in both replicates
	that is the value to use"""
	out=open(mergedoutfile, "w")
	for key, value in replicate1.items():
		if key in replicate2:
			if replicate1[key]>replicate2[key]:
				pair=key.split(",")
				out.write(ch+"\t"+pair[0]+"\t"+ch+"\t"+pair[1]+"\t"+str(replicate2[key])+"\n")
			else:
				pair=key.split(",")
				out.write(ch+"\t"+pair[0]+"\t"+ch+"\t"+pair[1]+"\t"+str(replicate1[key])+"\n")
	out.close()
	

merge_replicates(argv[1],read_replicate(argv[1],argv[2]), read_replicate(argv[1],argv[3]), argv[4])
			
