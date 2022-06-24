import random
import sys
def generate_sequence (length) :
	return ''.join(random.choice('AGCT') for _ in range(length))


numberOfReads=0
readsList = []
while len(readsList) < 1000000:
	read = generate_sequence(150)
	if read in readsList :
		print("reapead read: " , readsList)
	else :
		readsList.append(read)
	if len(readsList)%10000 == 0 :
		print("10000")	
f = open("/bigtemp/rgq5aw/GeneratedInput/genRefrenceGenome.txt", "a")
readCount = 0 
for read in readsList :
	sequence = ">chr" + str(readCount) + "\n" + read + "\n"
	readCount+=1
	f.write(sequence)
f.close()
