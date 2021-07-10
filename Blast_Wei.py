#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 16:41:18 2020

@author: vivi
"""
## Import the modules
import Bio.Blast
from Bio.Blast import NCBIWWW
import os
# Change the directory to file location
os.getcwd()
os.chdir('/Users/weiwang') ## please change it to your directory
fasta_string = open("MultipleBlast.txt").read()
print('read in data')
save_file = open("out.xml", "w")
print('output file is created')
result_handle = NCBIWWW.qblast("blastn","nt",fasta_string)
print("result is generated")

save_file.write(result_handle.read())
save_file.close()
# save multiple sequence blast results in file out


# code for parsing
from Bio.Blast import NCBIXML
from Bio.SearchIO._model import hsp

result_handle = open("out.xml",'r')
blast_records = NCBIXML.parse(result_handle)
blast_record = next(blast_records)

max_E_value= 0.04

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < max_E_value:
            print("****Alignment****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...") ## print out first 75 results
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
