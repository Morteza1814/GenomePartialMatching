import random
import sys

lines = []
f = open("/bigtemp/rgq5aw/GeneratedInput/genRefrenceGenome.txt", "r")
while True:
	line = f.readline()
  
	if not line:
		break
	if (line.find('chr') == -1):
		lines.append(line)	

print(len(lines))
f.close()

fo = open("/bigtemp/rgq5aw/GeneratedInput/genQuerySampled.fa", "a")
readCount = 0
for i in range (1,len(lines)):
	if i%1000==0:
		sequence = ">chr" + str(readCount) + "\n" + lines[i] + "\n"
        	readCount+=1
        	fo.write(sequence)
fo.close()

