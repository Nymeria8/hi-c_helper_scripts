#This script uses hi-c matrices and transform them into the input
#required by circos, taking of the interactions out of the interval
#given by the user
#the color and thickness of the interactions is defined by the number of interactions
#between the intervals

#usage: python script_circos.py infile, chromossome, start, end, interval, outfile

from sys import argv

colors={"1":"thickness=1p,color=lgrey","2":"thickness=2p,color=green","3":"thickness=3p,color=blue","4":"thickness=4p,color=black","5":"thickness=5p,color=red"}


def make_table(infile, ch, start, end):
	l=[]
	f=open(infile)
	for i in f:
		line=i.split("\t")
		if line[0]==ch and line[3]==ch:
			if line[1]>=start and line[2]<=end:
				l.append(i)
			elif line[4]>=start and line[5]<=end:
				l.append(i)
	f.close()
	return l
	
def make_table2(infile, ch, start, end):
	l=[]
	f=open(infile)
	for i in f:
		line=i.split("\t")
		if line[0]==ch and line[3]==ch:
			if line[1]>=start and line[2]<=end and line[4]>=start and line[5]<=end:
				l.append(i)
	f.close()
	return l


def read_table(l):
	if len(l)==0:
		print("\nThe window that yo chosed do not have interactions! Please use a larger interval.")
		exit()
	else:
		intt=int(l[1].split("\t")[2])-int(l[1].split("\t")[1])
		"""le a tabela de input com as interacoes"""
		interactions={}
		for i in l:
			line=i.split("\t")
			if line[1]==line[4] or int(line[1])+intt==int(line[4]) or  int(line[4])+intt==int(line[1]):
				pass
			else:
				chf="hs"+line[0]
				font=chf+" "+line[1]+" "+line[2]
				chs="hs"+line[3]
				#stop=chs+" "+line[4]+" "+line[5]+" thickness="+str(round(float(line[-1].strip("\n"))))+"p"
				if line[6].strip() in colors:
					stop=chs+" "+line[4]+" "+line[5]+" "+colors[line[6].strip()]
					if font in interactions:
						interactions[font].append(stop)
					else:
						interactions[font]=[stop]	
				else:
					stop=chs+" "+line[4]+" "+line[5]+" thickness=5p,color=red"
					if font in interactions:
						interactions[font].append(stop)
					else:
						interactions[font]=[stop]	
		
		return interactions
	
def write_for_circos(interactions, outfile):
	out=open(outfile, "w")
	for key, value in interactions.items():
		for element in value:
			out.write(key+" "+element+"\n")
	out.close()


if argv[1]=="all":
	write_for_circos(read_table(make_table(argv[2], argv[3], argv[4],argv[5])), argv[6])
else:
	write_for_circos(read_table(make_table2(argv[2], argv[3], argv[4],argv[5])), argv[6])		
			
