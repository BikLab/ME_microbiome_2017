#!/bin/bash

module load fastx_toolkit

nematoda_file=$1
nema_out_file=$2

fasta_formatter -i $1 -o $2
