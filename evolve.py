#BINF3020 Assignment 1 - Protein sequence evolution and alignment
#24/09/2022 -- Part 1 implement a protein sequence evolution simulator
#written by z5289835 Yuet Yat Chan
#
#-------------------------
# #Package used
#-------------------------
#Python --version 3.8.9
#Biopython --version 1.79
#-------------------------
 
import sys
from Bio import SeqIO
import re
import random


# Dictionary to convert AA to index for mutation matrix 
aaToIndex = {
    'A' : 0,  #alanine
    'R' : 1,  #arginine
    'N' : 2,  #asparagine
    'D' : 3,  #aspartate
    'C' : 4,  #cysteine
    'Q' : 5,  #glutamine
    'E' : 6,  #glutamate
    'G' : 7,  #glycine
    'H' : 8,  #histidine
    'I' : 9,  #isoleucine
    'L' : 10, #leucine
    'K' : 11, #lysine
    'M' : 12, #methionine
    'F' : 13, #phenylalanine
    'P' : 14, #proline
    'S' : 15, #serine
    'T' : 16, #threonine
    'W' : 17, #tryptophan
    'Y' : 18, #tyrosine
    'V' : 19, #valine
}

aaList = list(aaToIndex.keys())

# 2D array of the mutation matrix
pMatrix = [[9867,2,9,10,3,8,17,21,2,6,4,2,6,2,22,35,32,0,2,18],
           [1,9914,1,0,1,10,0,0,10,3,1,19,4,1,4,6,1,8,0,1],
           [4,1,9822,36,0,4,6,6,21,3,1,13,0,1,2,20,9,1,4,1],
           [6,0,42,9859,0,6,53,6,4,1,0,3,0,0,1,5,3,0,0,1],
           [1,1,0,0,9973,0,0,0,1,1,0,0,0,0,1,5,1,0,3,2],
           [3,9,4,5,0,9876,27,1,23,1,3,6,4,0,6,2,2,0,0,1],
           [10,0,7,56,0,35,9865,4,2,3,1,4,1,0,3,4,2,0,1,2],
           [21,1,12,11,1,3,7,9935,1,0,1,2,1,1,3,21,3,0,0,5],
           [1,8,18,3,1,20,1,0,9913,0,1,1,0,2,3,1,1,1,4,1],
           [2,2,3,1,2,1,2,0,0,9871,9,2,12,7,0,1,7,0,1,33],
           [3,1,3,0,0,6,1,1,4,22,9947,2,45,13,3,1,3,4,2,15],
           [2,37,25,6,0,12,7,2,2,4,1,9924,20,0,3,8,11,0,1,1],
           [1,1,0,0,0,2,0,0,0,5,8,4,9875,1,0,1,2,0,0,4],
           [1,1,1,0,0,0,0,1,2,8,6,0,4,9944,0,2,1,3,28,0],
           [13,5,2,1,1,8,3,2,5,1,2,2,1,1,9924,12,4,0,0,2],
           [28,11,34,7,11,4,6,16,2,2,1,7,4,3,17,9840,38,5,2,2],
           [22,2,13,4,1,3,2,2,1,11,2,8,6,1,5,32,9869,0,2,9],
           [0,2,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,9976,1,0],
           [1,0,3,0,3,0,1,0,4,1,1,0,0,21,0,1,1,2,9947,1],
           [13,2,1,1,3,2,2,3,3,57,11,1,17,1,3,2,10,0,2,9901],
]


#get input file as a list of fasta objects
#each fasta object as attribute to id, name, ... , sequence
input_fasta = list(SeqIO.parse(sys.stdin, "fasta"))


#error checking

#if the input does not contain exactly one sequence in FASTA format 
#or is not a fasta file
if len(input_fasta) != 1:
    sys.exit("Please input FASTA file with exactly one sequence!")

#get original sequence
original_seq = input_fasta[0].seq

#check if the seq is valid ()
for aminoAcids in original_seq:
    if aminoAcids not in aaToIndex:
        sys.exit("The sequence contains something other than the above 20 amino acid code letters.")

#first generation
print(">" + input_fasta[0].id + " | generation_0")
#keep line shorter than 80 characters
print(re.sub("(.{80})", "\\1\n", str(original_seq), 0, re.DOTALL))

#better variable naming
sequence = original_seq

#run for 500 generations
for generation in range (1,501):

    # first iteration uses input sequence, iterations after uses result sequence (mutated)
    if generation > 1:
        sequence = resultSeq
    #empty previous sequence
    resultSeq = ""

    #mutation loop for each generation of sequence
    for aminoAcids in sequence:
        #get the row number of the amino acid to access the probability
        row = aaToIndex.get(aminoAcids)
        
        #generate a random AA according to the respective porbability matrix row
        resultAA = random.choices(aaList, weights=pMatrix[row], k=1)
        #append the aa to the result
        resultSeq += resultAA[0]


    #print description of the sequence
    print(">" + input_fasta[0].id + " | generation_" + str(generation))
    #keep line shorter than 80 characters
    print(re.sub("(.{80})", "\\1\n", str(resultSeq), 0, re.DOTALL))