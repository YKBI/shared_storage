import sys,os

fi = sys.argv[1]
st1 = int(sys.argv[2])
ed1 = int(sys.argv[3])
#st2 = int(sys.argv[4])
#ed2 = int(sys.argv[5])

Lines = []
def ff(f,num1,num2):#,num3,num4):
	with open(f,'r') as F:
		for line in F.readlines():
			if line.startswith('ATOM'):
				if num1 <= int(line[22:26].strip()) <= num2:Lines.append(line)
				#elif num3 <= int(line[22:26].strip()) <= num4:Lines.append(line)
			if line.startswith('TER'):
				Lines.append(line)
				ter = int(line[7:12].strip())
				print ter
	with open(f,'r') as ff:
		for line in ff.readlines():
			if line.startswith('ATOM'):
				if ter <= int(line[7:12].strip()):Lines.append(line)
		Lines.append('END')
	with open(f +'_new','w') as W:
		for line in Lines:
			W.write(line)
ff(fi,st1,ed1)#,st2,ed2)	
