#uses 2 hi-c matrixes and makes a table with the interactions in comum and not in common
#TO DO: use more than 2 matrices
#usage:
#python cross_matrix.py cromossome start end matrix1 matrix2 outfile

from sys import argv
from itertools import combinations

def read_matrix(ch, start, end, infile):
	"""selects from the interaction matrices the lines corresponding
	to the interval defined by start and end. Returns a dictionary as
	start1 start2=number of interactions"""
	f=open(infile)
	l={}
	for i in f:
		line=i.split()
		if int(line[0])==ch:
			if int(line[1])>=start and int(line[2])<=end:
					l[line[1]+" "+line[4]]=line[-1].strip()
		elif int(line[3])==ch:
			if int(line[4])>=start and int(line[5])<=end:
					l[line[4]+" "+line[1]]=line[-1].strip()	
	f.close()
	return l

def cross (matrix1, matrix2):
	"""crosses two matrices and returns a dictionary as
	onlymatrix1=interval1, interval2...; matrix1,matrix2=interval1, interval2..."""
	data=dict(mt1=set(matrix1.keys()), mt2=set(matrix2.keys()))
	variations={}
	results={}
	for i in range(len(data)):
		for v in combinations(data.keys(),i+1):
			vsets=[ data[x] for x in v]
			variations[tuple(sorted(v))]=reduce(lambda x,y: x.intersection(y), vsets)
	for k,v in variations.items():
		results[str(k).strip("(),'")]=v
	return results

def see_results(results, matrix1, matrix2):
	"""uses the dictionay from cross, and see if the number of interactions
	match in both matrices. if not, the rest are added to the entries of the
	non matching matrices. Returns a dictionary"""
	count_int=0
	count_no1=0
	count_no2=0
	final={}
	for key, value in results.items():
		if len(key.split(","))>1:
			for el in value:
				if matrix1[el]==matrix2[el]:	
					count_int+=int()
				else:
					if int(matrix1[el])>int(matrix2[el]):
						a=int(matrix1[el])-int(matrix2[el])
						count_no1+=a
						count_int+=int(matrix1[el])-a
					else:
						a=abs(int(matrix1[el])-int(matrix2[el]))
						count_no2+=a
						count_int+=int(matrix2[el])-a
			final[key]=str(count_int)+"\t"+str(count_no1)+"\t"+str(count_no2)
		else:
			final[key]=len(value)
	for key, value in final.items():
		print(key)
		if len(key.split(","))>1:
			s=key.split(",")
			print(s[1][1:])
			final[s[0][:-1]]+=int(value.split("\t")[1])
			final[s[1][2:]]+=int(value.split("\t")[2])
	return final
	
def writer (final, outfile):
	"""writes the cross dictionary as
	interval    number of interactions"""
	out=open(outfile, "w")
	for key, value in final.items():
		if len(key.split(","))>1:
			out.write(key+"\t"+value.split("\t")[0]+"\n")
		else:
			out.write(key+"\t"+str(value)+"\n")
	out.close()

#######################Runs the program#############################
					
mat1=read_matrix(argv[1], argv[2], argv[3], argv[4])
mat2=read_matrix(argv[1], argv[2], argv[3], argv[5])	
writer(see_results(cross(mat1,mat2),mat1,mat2),argv[6])
