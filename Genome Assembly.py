'''
Genome Assembly /Neta Sohlberg

In genome assembly many short sequences (reads) from a sequencing machine is assembled into long sequences –ultimately chromosomes.
 This is done by ordering overlapping reads so they together represent genomic sequence.
 For example given these three reads:
 AGGTCGTAG, CGTAGAGCTGGGAG, GGGAGGTTGAAA
 ordering them based on their overlap like this:
  AGGTCGTAG
      CGTAGAGCTGGGAG
               GGGAGGTTGAAA
produces the genomic sequence:
    AGGTCGTAGAGCTGGGAGGTTGAAA
'''
#The function read a file and create a dictionary includes every line and its number
def readDataFromFile(fileName):
    f=open(fileName)
    dic={}
    #read line after line. Line numbers are the keys in dic and the sequences are the values.
    for line in f:
        dic[line[:1]]=line[2:].strip()
    f.close()
    return dic #return the dictionary
#end function

#The function return the mean sequence length in the file
def meanLength(fileName):
    #create a dictionary,using 'readDataFromFile' function
    dic=readDataFromFile(fileName)
    sum=0.0
    #sum the lengthes of all sequences
    for i in dic:
        sum+=len(dic[i])
    return sum/len(dic) #return the averge of the lengthes
#end function

# The function get 2 sequences and returns the overlapping sequence
def getOverlap(left, right):
    #search the overlapping sequence (i runs from len(left) to 0. 
    for i in range(len(left))[::-1]:
        #if there is a match- return the overlapping sequence 
        if right[:i]== left[len(left)-i:]:
            return right[:i]
    #If there is no overlap the function returns an empty string
    return ""
#end function

# The function get a dictionary of sequences and create a dictionary
# of dictionaries with the length of all overlaps, using 'getOverlap' function 
def getAllOverlaps(reads):
    #d is a dictionary of dictionaries
    d={}
    for i in reads:
        d[i]={}
        for j in reads:
        #every value in d is the length of the overlapping sequence between sequence 'i' and sequence 'j' 
            if i!=j:
                d[i][j]=len(getOverlap(reads[i],reads[j]))
    return d
#end function

# The function get a dictionary of dictionaries and print it into a matrix
def prettyPrint(overlaps):
    #print first line
    print ("%1s" %"",end=' ')
    for num in range (len(overlaps)):
        print ("% 3d" %  (num+1),end=' ')
    print()
    # print the matrix, according to the dictionary 
    for i in range (len(overlaps)):
        print  (i+1,end=' ')
        for j in range (len(overlaps)):
            if i==j:
                print ("%3s"% "-" , end=' ')
            else:
                print ("% 3d" % overlaps[str(i+1)][str(j+1)],end=' ')
        print()
#end fuction

# The function get a dictionary of dictionaries and find the first read sequence 
def findFirstRead(overlaps):
    d={}
    #All values in d start with 'true'.
    for key in overlaps:
        d[key]=True
    #If one value (or more) at 'j' colum is higher than 2, its value in d convert to 'false'
    for i in overlaps:
        for j in overlaps:
            if i!=j and overlaps[i][j]>2:
                d[j]=False
    #The function return the key in d that his value is still 'true'
    for key in overlaps:
        if d[key]==True: return key
#end function

#The function get a dictionary and return the key with the highest value
def findKeyForLargestValue(d):
    highest=''
    num=2
    for i in d:
        if d[i]>num:
            num=d[i]
            highest=i
    return highest
#end function

#The function find the order of the sequences. It is a recursive function   
def findOrder(name, overlaps):
    #using 'findKeyForLargestValue' function to find the next sequence
    nextName=findKeyForLargestValue(overlaps[name])
    #if it is the last sequence- stop and return
    if nextName=='': return [name]
    else:
        #else- continue to find the next sequence
        return [name] + findOrder(nextName,overlaps)
#end function

#The function get list argument readOrder containing the order of reads,
#a dictionary argument reads returned from readDataFromFile
#and a dictionary argument overlaps returned fromgetAllOverlaps.
#The function return a string with the the full DNA sequence
def assembleGenome(readOrder, reads, overlaps):
    Genomic_sequence=""
    Genomic_sequence+=reads[readOrder[0]]#insert the first sequence into 'Genomic_sequence'
    for i in range(len(readOrder)):#insert the all sequence by their order
        if i+1==len(readOrder): break
        a=overlaps[readOrder[i]][readOrder[i+1]]#a is the place in a sequence that is not overlap with the previous sequence
        Genomic_sequence+=reads[readOrder[i+1]][a:]#insert the next sequence from the 'a' place
    return Genomic_sequence #return the full Genomic sequence
#end function
            

#main:
d=readDataFromFile("genome_assembly.txt")
ol=getAllOverlaps(d)
l= findOrder(findFirstRead(ol), ol)
print ("sequences:")
for i in d:
    print( i,":",d[i])
print ('\n',"matrix:")
prettyPrint( ol)
print ('\n',"order:")
for j in l: print (j,end=' ')
print ('\n\n',"final sequence:")
print (assembleGenome(l,d,ol))

#example:
#sequences:
#1 : GGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTCGTCCAGACCCCTAGC
#3 : GTCTTCAGTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGTCGTGAACACATCAGT
#2 : CTTTACCCGGAAGAGCGGGACGCTGCCCTGCGCGATTCCAGGCTCCCCACGGG
#5 : CGATTCCAGGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTC
#4 : TGCGAGGGAAGTGAAGTATTTGACCCTTTACCCGGAAGAGCG
#6 : TGACAGTAGATCTCGTCCAGACCCCTAGCTGGTACGTCTTCAGTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGT
#
#matrix:
#    1   2   3   4   5   6
#1   -   1   0   0   1  29
#2  13   -   1   0  21   0
#3   0   0   -   1   0   1
#4   1  17   1   -   2   0
#5  39   1   0   0   -  14
#6   0   0  43   1   0   -
#
#order:
#4 2 5 1 6 3 
#
#final sequence:
#TGCGAGGGAAGTGAAGTATTTGACCCTTTACCCGGAAGAGCGGGACGCTGCCCTGCGCGATTCCAGGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTCGTCCAGACCCCTAGCTGGTACGTCTTCAGTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGTCGTGAACACATCAGT
