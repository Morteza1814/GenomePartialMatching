import sys

lines = []
f = open("/bigtemp/rgq5aw/uC_COLORADO_FILES/SRR.csv", "r")
fo = open("/bigtemp/rgq5aw/naiveBitVectorData/sampledDataFromSRR.fa", "a")

readCount = 0
queryCount=0
while readCount<100000:
        line = f.readline()
        if readCount%10==0:
                sequence = ">chr" + str(queryCount) + "\n" +  line + "\n"
                fo.write(sequence)
                queryCount+=1
        readCount += 1
f.close()
fo.close()
