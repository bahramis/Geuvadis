#!/usr/bin/env python

import os, sys

matrix=[]
for line in open(sys.argv[1]):
    line=line.rstrip('\n\r')
    a=line.split('\t')
    if (matrix==[]):
        matrix=a
    else:
        if (len(a)!=len(matrix)):
            print ['error', len(a), len(matrix)]
        for i in range(len(a)):
            matrix[i]+='\t'+a[i]

output=open(sys.argv[2], 'w')
for line in matrix:
    output.write(line+'\n')

