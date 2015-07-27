#Testes in a file with interactions, if two sets of user given coordinates
#Have interactions where one interval is target of the other, and vice-versa
#Returns the number of reciprocate interactions



from sys import argv

f=open(argv[1])

interval1=[int(argv[2]), int(argv[3])]
interval2=[int(argv[4]), int(argv[5])]

for i in f:
	ff=i.split(" ")
	rr=0
	if int(ff[1])>= interval1[0] and  int(ff[2])<= interval1[1] and int(ff[4])>= interval2[0] and  int(ff[5])<= interval2[1]:
		rr+=1
	elif int(ff[4])>= interval1[0] and  int(ff[5])<= interval1[1] and int(ff[1])>= interval2[0] and  int(ff[2])<= interval2[1]:
		rr+=1

f.close()

print(rr)
