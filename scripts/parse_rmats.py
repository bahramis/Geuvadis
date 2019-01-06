#!/usr/bin/env python

import sys,os
import math
from collections import defaultdict

def averageCountSample1(record):
    return float(sum(record.IJC_SAMPLE_1) + sum(record.SJC_SAMPLE_1))/len(record.IJC_SAMPLE_1)

class PSL:
    def __init__(self, line):
        (self.ID,
         self.IJC_SAMPLE_1,
         self.SJC_SAMPLE_1,
         self.IJC_SAMPLE_2,
         self.SJC_SAMPLE_2,
         self.IncFormLen,
         self.SkipFormLen,
         self.IncLevel1,
         self.IncLevel2,
         self.IncLevelDifference) = line.strip().split()
        self.ID = int(self.ID)
        self.IJC_SAMPLE_1 = [int(x) for x in self.IJC_SAMPLE_1.split(",")]
        self.SJC_SAMPLE_1 = [int(x) for x in self.SJC_SAMPLE_1.split(",")]
        self.IJC_SAMPLE_2 = [int(x) for x in self.IJC_SAMPLE_2.split(",")]
        self.SJC_SAMPLE_2 = [int(x) for x in self.SJC_SAMPLE_2.split(",")]
        IncLevel1_temp = []
        for x in self.IncLevel1.split(",") :
            if (x == "NA") :
                IncLevel1_temp.append("NA")
            else :
                IncLevel1_temp.append(float(x))
        self.IncLevel1 = IncLevel1_temp
        IncLevel2_temp = []
        for x in self.IncLevel2.split(",") :
            if (x == "NA") :
                IncLevel2_temp.append("NA")
            else :
                IncLevel2_temp.append(float(x))
        self.IncLevel2 = IncLevel2_temp

sample_population = defaultdict()
all_exons = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
total_counts = {'CEU': 89, 'FIN': 92, 'GBR': 86, 'TSI': 91, 'YRI': 87}

for line in open(sys.argv[1]):
    line = line.strip()
    if line == "" : continue
    sample = line.split("\t")[0]
    population = line.split("\t")[2]
    sample_population[sample] = population

sample_order1 = []
for line in open(sys.argv[2]):
    for x in line.strip().split(","):
        sample_order1.append(os.path.splitext(os.path.basename(x))[0])

sample_order2 = []
for line in open(sys.argv[3]):
    for x in line.strip().split(","):
        sample_order2.append(os.path.splitext(os.path.basename(x))[0])


foundStart = False
fh = open(sys.argv[4])
for i in fh:
    if foundStart:
        x = PSL(i)
        all_exons['CEU'][x.ID]['zero_count'].append(0)
        all_exons['FIN'][x.ID]['zero_count'].append(0)
        all_exons['GBR'][x.ID]['zero_count'].append(0)
        all_exons['TSI'][x.ID]['zero_count'].append(0)
        all_exons['YRI'][x.ID]['zero_count'].append(0)
        for j in range(len(x.IJC_SAMPLE_1)):
            iindx = sample_population[sample_order1[j]]
            if (x.IncLevel1[j] != "NA") :
                all_exons[iindx][x.ID]['inclusion_count'].append(x.IJC_SAMPLE_1[j])
                all_exons[iindx][x.ID]['skipping_count'].append(x.SJC_SAMPLE_1[j])
                all_exons[iindx][x.ID]['psi'].append(x.IncLevel1[j])
            else :
                all_exons[iindx][x.ID]['zero_count'][0] += 1
 
        for j in range(len(x.IJC_SAMPLE_2)):
            iindx = sample_population[sample_order2[j]]
            if (x.IncLevel2[j] != "NA") :
                all_exons[iindx][x.ID]['inclusion_count'].append(x.IJC_SAMPLE_2[j])
                all_exons[iindx][x.ID]['skipping_count'].append(x.SJC_SAMPLE_2[j])
                all_exons[iindx][x.ID]['psi'].append(x.IncLevel2[j])
            else :
                all_exons[iindx][x.ID]['zero_count'][0] += 1

    foundStart = True

fh.close()

filtered_exons = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

for i in all_exons:
    for j in all_exons[i]:
        if len(all_exons[i][j]['psi']) > 0 :
            average_psi = float(sum(all_exons[i][j]['psi']))/float(len(all_exons[i][j]['psi']))
            average_total_read_count = float(sum(all_exons[i][j]['inclusion_count']) + sum(all_exons[i][j]['skipping_count'])) / float(len(all_exons[i][j]['inclusion_count']))
            range_of_exon_inclusion = float(max(all_exons[i][j]['psi'])) - float(min(all_exons[i][j]['psi']))
            proportion_zero_count_individuals = float(all_exons[i][j]['zero_count'][0]) / float(total_counts[i])
            if (average_psi >= 0.05 and average_psi <= 0.95 and average_total_read_count >= 10 and range_of_exon_inclusion >= 0.1 and proportion_zero_count_individuals <= 0.8) :
                print str(i) + "\t" + str(j)

