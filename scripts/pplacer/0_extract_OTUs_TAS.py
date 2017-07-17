#!/usr/bin/python3
# A program for extracting nematoda OTU IDs
# USAGE: python 0_extract_OTUs_TAS.py input_file output_file
# Author: Taruna Aggarwal
# Affiliation: University of California, Riverside
# Date: 07/14/2017


import sys

inFile = open(sys.argv[1], "r")
outFile = open(sys.argv[2], "w")


for currentLine in inFile:
    currentLine = currentLine.rstrip()
    if "__Nematoda;" in currentLine:
        OTU = currentLine.split("\t")[0]
        outFile.write("{0}{1}".format(OTU, "\n"))

inFile.close()
outFile.close()
