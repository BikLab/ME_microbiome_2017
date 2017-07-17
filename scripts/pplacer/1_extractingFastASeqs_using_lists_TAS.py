#!/usr/bin/python3
# A program for pulling out seqs from a fasta file using a list of headers or parts of a header
# USAGE: python 1_extractingFastASeqs_using_lists_TAS.py list input output
# Author: Taruna Aggarwal
# Affiliation: University of California, Riverside
# Date: 04/13/2017
import sys

inFile = open(sys.argv[1], "r")
readFile = open(sys.argv[2], "r")
outFile = open(sys.argv[3], "w")

header_list = []
for currentLine in inFile:
    currentLine = currentLine.rstrip()
    header_list.append(currentLine)

flag = False
for eachLine in readFile:
    eachLine = eachLine.rstrip()
    if eachLine.startswith(">"):
        flag = False
        header = eachLine.split(" ")[1]
        if header[1:] in header_list:
            outFile.write(eachLine + "\n")
            flag = True
    else:
        if flag:
            outFile.write(eachLine + "\n")

inFile.close()
readFile.close()
outFile.close()

