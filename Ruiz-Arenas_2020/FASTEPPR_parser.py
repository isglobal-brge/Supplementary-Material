# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 16:16:03 2016

@author: marta
"""

import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Chromosome parser")

parser.add_argument('-i', '--input_file',
                    dest = "input",
                    action = "store",
                    default = False, 
                    help = "Input file")
                    
parser.add_argument('-o', '--output_file',
                    dest = "output",
                    action = "store",
                    default = False, 
                    help = "Output file")
                
     
options = parser.parse_args()

input_file = options.input
output_file = options.output

#input_file = "chr_17"
#output_file = "chr_17_out"


infile = open(input_file, "r")

counter = 0
for line in infile:
    counter += 1
    
    
matrix = np.empty([int(counter)/2 + 1, 5], dtype="a10")
matrix[0,] = ["Start", "End", "Rho", "CIL", "CIR"]

infile = open(input_file, "r")

i = 0

for line in infile:
    
    line = line.strip()
    
    if line.startswith("Position"):
        i = i + 1
        matrix[i][0] = line.split()[1].split("-")[0]
        matrix[i][1] = line.split()[1].split("-")[1][0:-1]
        
        
    else:
        
        matrix[i][2] = line.split()[0].split(":")[1]
        try:
            matrix[i][3] = line.split()[1].split(":")[1]
            matrix[i][4] = line.split()[2].split(":")[1]
        except:
            matrix[i][3] = "NA"
            matrix[i][4] = "NA"

#print matrix      
np.savetxt(output_file, matrix, delimiter = "\t", fmt = "%s")
        

