#!/usr/bin/python3
# A program for removing spaces from sequences
# USAGE: python 4_removing_spaces_in_seqs_TAS.py input_file output_file
# Author: Taruna Aggarwal
# Affiliation: University of California, Riverside
# Date: 04/13/2017

import sys

inFile = open(sys.argv[1], "r")
outFile = open(sys.argv[2], "w")

flag = False
for currentLine in inFile:
    currentLine = currentLine.rstrip()
    if currentLine.startswith(">"):
        flag = False
        outFile.write("{0}{1}".format(currentLine, "\n"))
        flag = True
    else:
        if flag:
            currentLine = currentLine.replace(" ", "")
            outFile.write("{0}{1}".format(currentLine, "\n"))

inFile.close()
outFile.close()
