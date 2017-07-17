#!/usr/bin/python3
# A program for counting the len of each line in a file
# USAGE: python 2_count_lines.py input_file
# Author: Taruna Aggarwal
# Affiliation: University of California, Riverside
# Date: 04/13/2017

import sys

inFile = open(sys.argv[1], "r")

for currentLine in inFile:
    currentLine = currentLine.rstrip()
    currentLine = currentLine.replace(" ", "")
    length = len(currentLine)
    print(length)
