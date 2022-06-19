import random
import sys

lines = []
f = open("/bigtemp/rgq5aw/GeneratedInput/genRefrenceGenome.fa", "r")
while True:
	line = f.readline()
  
	if not line:
		break
	if (line.find('chr') == -1):
		lines.append(line)	

print(len(lines))
f.close()

fo = open("/bigtemp/rgq5aw/GeneratedInput/genQuerySampledOverlapping.fq", "a")
readCount = 0
for i in range (1,len(lines)):
	if i%1000==0:
		query =  ("A" * 100) + lines[i][100:150]
		sequence = ">chr" + str(readCount) + "\n" + query + "\n"
#		print("query and its raq format")
#		print(lines[i])
#		print(query)
        	readCount+=1
        	fo.write(sequence)
fo.close()

