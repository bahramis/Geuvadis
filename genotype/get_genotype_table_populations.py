#!/usr/bin/env python

import sys, os

population={}
for line in open(sys.argv[1]):
    line=line.rstrip('\n\r')
    a=line.split('\t')
    sample=a[0]
    popu=a[2]
    if (not population.has_key(popu)):
        population[popu]=[sample]
    elif(not sample in population[popu]):
        population[popu].append(sample)
    else:
        pass


gender={}
for line in open(sys.argv[2]):
    line=line.rstrip('\n\r')
    a=line.split('\t')
    if (a[0]=='Family.ID'):
        continue
    gender[a[1]]=a[4]

chrom = int(sys.argv[3])
population_name = sys.argv[4]

indgeno={}
for sample in population[population_name]:
    indgeno[sample]=[]

samples=[]
for line in open('ALL.chr'+str(chrom)+'.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf'):
    line=line.rstrip('\n\r')
    a=line.split('\t')
    if (a[0]=='#CHROM'):
        samples=a[9:]
        break


mapout=open(str(population_name) + "/chr" + str(chrom) + "_" + str(population_name) + "_allsnp.map", 'w')
pedout=open(str(population_name) + "/chr" + str(chrom) + "_" + str(population_name) + "_allsnp.ped", 'w')

pedout.write('\t'.join( population[population_name] )+'\n')
pedout.write('\t'.join(['1']*len(population[population_name]))+'\n')
allsex=[]
for sample in population[population_name]:
    allsex.append(gender[sample])

pedout.write('\t'.join(allsex)+'\n')

print 'processing chr'+str(chrom)
for line in open('ALL.chr'+str(chrom)+'.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf'):
    line=line.rstrip('\n\r')
    if (line[:2]=='##'):
        continue
    a=line.split('\t')
    if (a[0]=='#CHROM'):
        continue
    pos=int(a[1])
    snpid=a[2]
    ref=a[3]
    alt=a[4]
    qual=a[5]
    fil=a[6]
    if (fil!= 'PASS' or len(ref)!=1 or len(alt)!=1):
        continue
    if (str(chrom) !=a[0]):
        print ['chr error', chrom, a[0]]
    alt_count=0
    ref_count=0
    rsgeno={}
    for i in range(9, len(a)):
        sample=samples[i-9]
        if(not sample in population[population_name]):
            continue
        gt=a[i].split(':')[0]
        if (gt=='0|0' or gt=='0/0'):
            gt='1 1'
            ref_count+=2
        elif(gt=='1|1' or gt =='1/1'):
            gt='2 2'
            alt_count+=2
        elif (gt=='0|1' or gt =='1|0' or gt=='0/1' or gt=='1/0'):
            gt='1 2'
            ref_count+=1
            alt_count+=1
        else:
            gt='0 0'
        rsgeno[sample]=gt
    allgt=[]
    for sample in population[population_name]:
        allgt.append(rsgeno[sample])

    pedout.write('\t'.join(allgt)+'\n')
    mapout.write('\t'.join([str(chrom), snpid, '0', str(pos)])+'\n')
 
