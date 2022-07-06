import sys
import time

KMER_LENGTH = 10


readFile = open("/bigtemp/rgq5aw/batFiles/SRR1246727.fasta", "r")

start_time = time.time()
kmers = dict()
i=0
while True:
        line = readFile.readline()

        if not line:
                break
        if (line.find('SRR') != -1):
                continue
	i+=1

	if i%10==0:
		print(i)
	j=0
	while j<=(len(line)-KMER_LENGTH):
		kmer = line[j:j+KMER_LENGTH]
		if kmer in kmers.keys():
			kmerCount = kmers[kmer]
			kmerCount += 1
			kmers[kmer] = kmerCount
		if kmer not in kmers.keys():
			kmers[kmer] = 1
		j+=1

print("len of kmers map : " + str(len(kmers)))
print("--- %s seconds ---" % (time.time() - start_time))
