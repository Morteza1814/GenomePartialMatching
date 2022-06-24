import sys
import time

KMER_LENGTH = 0
ALLOWED_EDIT_DISTANCE = 0

def compare(a,b):
	numberOfMisMatches=0
	for i in range(len(a)):
		if a[i] !=  b[i]:
			numberOfMisMatches+=1
	return numberOfMisMatches

start_time = time.time()

args = sys.argv
KMER_LENGTH = int(args[1])
ALLOWED_EDIT_DISTANCE = int(args[2])

print("kmer_length is : " + str(KMER_LENGTH))
print("allowed edit distance is : " + str(ALLOWED_EDIT_DISTANCE))


readFile = open("/bigtemp/rgq5aw/naiveBitVectorData/sampledDataFromSRR.fa", "r")
queryFile = open("/bigtemp/rgq5aw/naiveBitVectorData/sampledQueryDataFromSRR.fa", "r")
#fo = open("/bigtemp/rgq5aw/naiveBitVectorData/kmerSize32FP.txt","a")

reads = dict()
key=''
value=''
while(True):
	line = readFile.readline()
	if 'chr' in line and key=='':
		 key = line[1:len(line)]

        if 'chr' not in line and len(line)!=0 and key!='':
                value = line

	if key!='' and value!='':	
		reads[key] = value
		key=''
		value=''
	
	if not line:
		break

i=0
for value in reads: 
	i+=1
print(i)

queries = dict()
key=''
value=''
while(True):
        line = queryFile.readline()
        if 'chr' in line and key=='':
                key = line[1:len(line)]

        if 'chr' not in line and len(line)!=0 and key!='':
                value = line

        if key!='' and value!='':
                queries[key] = value
		key=''
		value=''

        if not line:
                break

i=0
for key in queries:
	i+=1
print(i)

readFile.close()
queryFile.close()

numberOfBitVectorFullMatches=0
numberOfTrueFullMatches=0
queriesEvaluated=0
numberOfBitVectorPartialMatches=0
numberOfTruePartialMatches=0
for queryInd,query in queries.items():
	#print("query=\n" + query)
	for readInd, read in reads.items():
		i=0
		numberOfKmersFound=0
		while i<(len(query)-KMER_LENGTH):
			kmer = query[i:i+KMER_LENGTH]
			if kmer in read:
				numberOfKmersFound+=1
			i+=1
		if numberOfKmersFound==((len(query)-1)-KMER_LENGTH+1):
			numberOfBitVectorFullMatches+=1
                if query==read:	
			numberOfTrueFullMatches+=1
		if numberOfKmersFound==((len(query)-1)-KMER_LENGTH+1-(ALLOWED_EDIT_DISTANCE*KMER_LENGTH)):
			numberOfBitVectorPartialMatches+=1
		if compare(read, query)==ALLOWED_EDIT_DISTANCE:
			numberOfTruePartialMatches+=1
					
		#if numberOfKmersFound==((len(query)-1)-KMER_LENGTH+1) and query!=read:
			#fo.write("-------------------------")
			#fo.write("query" + query)
			#fo.write("read" + read)
	queriesEvaluated+=1
	print(queriesEvaluated)

print("--- %s seconds ---" % (time.time() - start_time))
print("numberOfBitVectorFullMatches = "+str(numberOfBitVectorFullMatches))
print("numberOfTrueFullMatches = "+str(numberOfTrueFullMatches))
print("numberOfBitVectorPartialMatches = "+str(numberOfBitVectorPartialMatches))
print("numberOfTruePartialMatches = "+str(numberOfTruePartialMatches))
