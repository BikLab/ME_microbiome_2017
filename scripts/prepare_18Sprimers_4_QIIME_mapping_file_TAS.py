#!/usr/bin/python3
# A program for creating primer combo labels and outputting corresponding sequences
# USAGE: python format_primer_files.py forward reverse intermediate_file final_file
# Author: Taruna A. Schuelke (Bik lab)
# Affiliation: University of California, Riverside
# Date: 06/13/2017

import sys

# create file variables
forward = open(sys.argv[1], "r")
reverse = open(sys.argv[2], "r")
output = open(sys.argv[3], "w")
final = open(sys.argv[4], "w")


# create a list of forward primer IDs
forward_list = []
for currentLine in forward:
    currentLine = currentLine.split("\t")[0]
    forward_list.append(currentLine)
#print(forward_list)

# create a list of reverse primer IDs
reverse_list = []
for currentLine in reverse:
    currentLine = currentLine.split("\t")[0]
    reverse_list.append(currentLine)
#print(reverse_list)

# create a multidimensional list containing possible forward and reverse primer ID combos
combos=[]
for i in forward_list:
    for j in reverse_list:
        combos.append([i,j])
#print(combos)
#print(len(combos))

# create separate columns with primer IDs and names for primer mixes
count = 0
for item in combos:
    f = item[0].split("-")[1]
    f = f.split("_")[0]
    f = f.replace("bc", "F")
    r = item[1].split("-")[1]
    r = r.split("_")[0]
    r = r.replace("bc", "R")
    output.write("18S.{0}\t{1}\t{2}\t{3}.{4}.{5}\n".format(count, item[0], item[1], "18S", f, r))
    count += 1

# close files
forward.close()
reverse.close()
output.close()

# create a dictionary for forward primer IDs and their corresponding barcodes
forward_barcodes = {}
forward_primer = {}
with open(sys.argv[1], "r") as infile:
    for eachLine in infile:
        f_primer_id = eachLine.split("\t")[0]
        f_barcode = eachLine.split("\t")[1]
        f_primer_seq = eachLine.split("\t")[2]
        forward_barcodes[f_primer_id] = f_barcode

# create a dictionary for reverse primer IDs and their corresponding barcodes
reverse_barcodes = {}
reverse_primer = {}
with open(sys.argv[2], "r") as infile:
    for eachLine in infile:
        r_primer_id = eachLine.split("\t")[0]
        r_barcode = eachLine.split("\t")[1]
        r_primer_seq = eachLine.split("\t")[2]
        reverse_barcodes[r_primer_id] = r_barcode

# open the intermediate file and compare primer IDs to primer dictionaries generated above
with open(sys.argv[3], "r") as file:
    myList = []
    for line in file:
        line = line.strip("\n")
        myList = line.split("\t")
        f_primer = myList[1]
        r_primer = myList[2]
        if f_primer in forward_barcodes.keys():
            F_barc = forward_barcodes[f_primer]
            myList.append(F_barc)
            myList.append(f_primer_seq)
            if r_primer in reverse_barcodes.keys():
                R_barc = reverse_barcodes[r_primer]
                myList.append(R_barc)
                myList.append(r_primer_seq)
                #sys.stdout.write("\t".join(myList))
                newLine = "\t".join(myList)
                final.write("{0}{1}".format(newLine, "\n"))
