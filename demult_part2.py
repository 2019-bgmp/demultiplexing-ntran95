#!/usr/bin/env python

print("The following table is in correspondance to Part 1, Question 1")
from prettytable import PrettyTable

x = PrettyTable()

x.field_names = ["File Name", "File Type"]
x.add_row(["1294_S1_L008_R1_001.fastq.gz", "Read 1"])
x.add_row(["1294_S1_L008_R4_001.fastq.gz", "Read 2"])
x.add_row(["1294_S1_L008_R2_001.fastq.gz", "Index 1"])
x.add_row(["1294_S1_L008_R3_001.fastq.gz", "Index 2"])

print(x)

#the following is answering Part 1, Question 2, c:
print("The number of indexes that have undetermined (N) base calls for R2 is: 3976613")
print("The number of indexes that have undetermined (N) base calls for R3 is: 3328051")


#--------------------------------------------------------------------------------------------------
'''This function is going to impliment the arg parse library so that I can set my quality score cut off, specify files to read in'''

#def get_arguements():
#parser = argparse.ArgumentParser(description="settting quality score cut off")
#parser.add_argument("-q", "--qual_score_cutoff", help="this argument specifies the quality score cut off", type =int, required=True)
#    parser.add_argument("-r1", "--read1", help="this argument specifies the filename", type =str, required=True)
##    parser.add_argument("-r2", "--read2", help="this argument specifies the filename", type =str, required=True)
##    parser.add_argument("-r3", "--read3", help="this argument specifies the filename", type =str, required=True)
##    parser.add_argument("-r4", "--read4", help="this argument specifies the filename", type =str, required=True)


#return parser.parse_args()

#args=get_arguements()
#qual_score_cutoff=args.qual_score_cutoff
#read1=args.read1
#read2=args.read2
##read3=args.read3
##read4=args.read4




'''Problem/Objective: given 4 files: 2 biological reads (R1 & R2), 2 index reads (R2, R3). Look through the biological files, use the index files as reference, and store the reads into their respective category:
(24 known indeces) x 2 (forward and reverse) index-paired = 48 fastq files
2 (forward and reverse) index hopped
2 (forward and reverse) unmatched/low qual
'''

#the two low qual files do not match the index and do not meet the qual score cut off
# my qual score cut off, based on part 1 is: 25, meaning all mean qual with the min 25 is kept (basically I want all my reads)

# I want to read in the two biological files and export the reads that are index-paired, index-hopped, and low quality (containing unknown N's or just low quality)'''

#since I made an arg parse function, I can read in all four files (R1-R4) simulataneously

#maybe make a function that executes the reverse compliment index/barcode and seq1 and seq2
#-------------------------------------------------------------------------------------------------
#create an empty list that holds all known index
#known_index_dict=[]
#when I read in the R1, R2, R3, R4 files, I'm going to compare the strings to the sequence in the list/dict of known indexes

#read in the indexes.txt file
#extract columns 4 and 5: index and index sequences
#store the ouputs into a dictionary
#make the index (i.g B1) the keys, make index sequences the value (i.g GTAGCGTA)

#-------------------------------------------------------------------------------------------------
''' Creating a function that converts ASCII symbols to Phred Scores '''
#will call this function when we want to see if the quality score is too low
#def convert_phred(letter):

    #return (ord(letter)) - 33


#test examples

#test: !
#return: 0
#-------------------------------------------------------------------------------------------------
'''this function is going to create the compliment string from the sequence lines,call this function when  '''
#def complement(seq):
    #create emply dictionary called complement={}
    #create a dictionary of compliment bases: {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    #make a variable called bases, set equal to seq
    #iterate through bases string and match compliment base inside dictionary
    #bases= [complement[base] for base in bases]
    #use the .join() function to turn 'bases' into a string again
    #return ''.join(bases)

#------ testing example-------
#test_sequence = "AACCTT"
#print(complement(test_sequence))

#output should be: TTGGAA

#-------------------------------------------------------------------------------------------------
#this function builds off of complement function, we will need to call this function when finding the reverse complement of R2 in R3
'''This function is used to create the reverse compliment of a sequence by calling the compliment function'''
#def reverse(sequence):
    #return compliment(seq[::-1])

#use print(reverse(variable_name)) to call the function, you can put any variable inside the reverse(__) to use this function
#------------------------------------------------------------------------------------------------
#reading in all four files simulataneously, with open(read1, read2, read3, read4) as file_handle:

    #I want to read in each line of the four fastq files and extract the sequence line to fulfill my conditions
    #While True: (I want to execute everything inside this while loop)
        ##I want to use .readline() to extract sequence line
        #L1 = file_handle.readline().strip()
        #if L1 == "":
        #    break
        #L2 = file_handle.readline().strip()
        #print(L2)
        #L3 = file_handle.readline().strip()
        #L4 = file_handle.readline().strip()

        #extract the sequence line
        #barcode = L2

        #this next section is going to use if conditions to match the indexes of R2 and R3 to the known indexes and extract the corresponding biological reads
        #if R1 has an index (R2) that matches a known barcode, and R4 has an index (R3) that is the reverse complement of the first index,
        #append the sequence of R2 and R3 to the R1's header, and append the same sequence of R2 and R3 to R4's header

        #if barcode of R2 in known_index_dict:
            #call reverse_complitment(barcode) function to see if R4 has a R3 index that is the reverse complement of first index
            #if R3 = reverse_complitment(barcode of R2):
                #append whole record of R1 into a new forward file corresponding to the barcode found, use a+ to continue appending over each iteration
                #append barcodes in R2 and R3 into each record's header
                #file_out = open(known_barcode+".fw", "a+"
                #print(L1, +"barcode in R2" + "barcode in R3" file = file_out)
                #print(L2, file = file_out)
                #print(L3, file = file_out)
                #print(L4, file = file_out)
                #close file


                #append whole record of R4 into a new reverse file corresponding to the known_barcode, use a+ to continue appending over each iteration
                #append barcodes in R2 and R3 into each record's header
                #file_out = open(known_barcode+".rv", "a+"
                #print(L1, +"barcode in R2" + "barcode in R3" file = file_out)
                #print(L2, file = file_out)
                #print(L3, file = file_out)
                #print(L4, file = file_out)
                #close file


        #this next section will set conditions for index-hopping
        #set scenerio where R2 matches known_index_dict but does not have a reverse complement in R3
        #elif barcode of R2 in known_index_dict:
            #if barcode in R3 != reverse_complitment(barcode in R2)
                #append whole record of R1 into a new forward, index-hopped file, use a+ to continue appending over each iteration
                #append barcodes in R2 and R3 into each record's header
                #file_out = open("index-hopped.fw", "a+")
                #print(L1, +"barcode in R2" + "barcode in R3" file = file_out)
                #print(L2, file = file_out)
                #print(L3, file = file_out)
                #print(L4, file = file_out)
                #close file

                ##append whole record of R4 into a new reverse, index-hopped file, use a+ to continue appending over each iteration
                #append barcodes in R2 and R3 into each record's header
                #file_out = open("index-hopped.rv", "a+"
                #print(L1, +"barcode in R2" + "barcode in R3" file = file_out)
                #print(L2, file = file_out)
                #print(L3, file = file_out)
                #print(L4, file = file_out)
                #close file

        #this next section will set conditions to find the unknown indexes
        #set condition where R2 indexes do not match the indexes in the known_index_dict
        #elif barcode not in known_index_dict:
            ##append whole record of R1 into a new forward, unknown file, use a+ to continue appending over each iteration
            #append barcodes in R2 and R3 into each record's header
            #file_out = open("unknown.fw", "a+")
            #print(L1, +"barcode in R2" + "barcode in R3" file = file_out)
            #print(L2, file = file_out)
            #print(L3, file = file_out)
            #print(L4, file = file_out)
            #close file

            ##append whole record of R4 into a new reverse, unknown file, use a+ to continue appending over each iteration
            #append barcodes in R2 and R3 into each record's header
            #file_out = open("unknown.rv", "a+")
            #print(L1, +"barcode in R2" + "barcode in R3" file = file_out)
            #print(L2, file = file_out)
            #print(L3, file = file_out)
            #print(L4, file = file_out)
            #close file


        #this next section will set conditions to extract low quality reads
        #first, set for loop to iterate (i) over each bases in line 4:
            #call convert_phred(i) to convert each base into phred score
            #calculate mean of each barcode phred score
            #store mean barcode into variable: qual_score
        #set condition where R2 indexes and R3 contains "N" or if qual_score < qual_score_cutoff, append to new files
        #elif "N" in barcode and qual_score < qual_score_cutoff:
            #append whole record of R1 into the forward, unknown file, use a+ to continue appending over each iteration
            #append barcodes in R2 and R3 into each record's header
            #file_out = open("unknown.fw", "a+")
            #print(L1, +"barcode in R2" + "barcode in R3" file = file_out)
            #print(L2, file = file_out)
            #print(L3, file = file_out)
            #print(L4, file = file_out)
            #close file


            #append whole record of R4 into the reverse, unknown file, use a+ to continue appending over each iteration
            #append barcodes in R2 and R3 into each record's header
            #file_out = open("unknown.rv", "a+")
            #print(L1, +"barcode in R2" + "barcode in R3" file = file_out)
            #print(L2, file = file_out)
            #print(L3, file = file_out)
            #print(L4, file = file_out)
            #close file
