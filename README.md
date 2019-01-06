# Running GLiMMPS on Geuvadis RNA-seq data

There are several steps to prepare the input for GLiMMPS:

Preparing the genome and gene annotation:
---------------------------------------
The hg19 genome was downloaded from:

http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/ENSEMBL/homo_sapiens/ENSEMBL.homo_sapiens.release-75/

The ENSEMBL gene annotation was downloaded from here:

http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

These files have to be put in genome_dir/


Downloading, preprocessing and Mapping the RNA-seq data:
---------------------------------------

The data is part of RNA-sequencing of 465 lymphoblastoid 
cell lines from the 1000 Genomes downloaded from:

https://www.ebi.ac.uk/ena/data/view/PRJEB3366

There are 445 individuals that we have the genotype for,
these individuals, the population they belong to, and the 
links to download the fastq files are in the metadata.txt 
file. This file is used for all of script generation in 
running the pipleine.

The RNA-seq data for some individuals contains reads that 
are 76 nucleotides and some of them have reads that are 
75 nucleotides. The preprocessing step before mapping cuts 
all the 76 nucleotide read down to 75 for proper splicing 
analysis where the read length is required to be equal for 
all the samples.

To downaload, preprocess and map the data to the hg19 genome 
just run:

```bash
$ python scripts/run_preprocess.py metadata.txt
```

It generates all the scripts for all the RNA-seq data to be 
run on cluster.

running rMATS
---------------------------------------
The latest version of rMATS was ran on the data. To do this 
simply run:

```bash
$ python scripts/run_splicing.py metadata.txt
```

This generates a script to run rMATS on all the 445 RNA-seq 
data.

Filtering exons
---------------------------------------
We filtered out the alternative splicing events in each population 
separately, according to the following criteria: 
(1) the average exon inclusion level across all individuals is 
between 5% and 95% 
(2) the average total read count of all individuals is no less 
than 10 
(3) the range (maximum - minimum) of exon inclusion levels across all
individuals is greater than 20% 
(4) more than 20% of the individuals have non-zero read count. 

To do this simply run:

```bash
$ python scripts/parse_rmats.py metadata.txt splicing/b_1.txt splicing/b_2.txt splicing/output/SE.MATS.JC.txt > splicing/filtered_exons.txt
```

The output is a tab delimited file with two columns, first column 
is the population and the second column is the exon IDs that passed 
all the criteria.

then run the following commands:

```bash
$ for {CEU,FIN,GBR,TSI,YRI}; do grep ${f} splicing/filtered_exons.txt | cut -f2 > glimmps/${f}/input/filtered_exons.txt
```

```bash
$ scripts/exon_info_commands.sh
$ scripts/exon_commands.sh
```

Preparing genotypes
---------------------------------------
First download the phase3 genotype from:

http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

The files have to be in the genotype directory. Extract the 
gzipped .vcf files and run the following script:

```bash
for f in {CEU,FIN,GBR,TSI,YRI}; do for i in {1..22}; do python scripts/get_genotype_table_populations.py genotype/sample_information.txt genotype/sample_pedigree_from_1kgenome.txt ${i} ${f}; done; done
```

Then to creat the .raw files run plink. This can be done by running:

```bash
$ python scripts/prepare_genotype.py metadata_genotype.txt
```

Running GLiMMPS
---------------------------------------
At this stage all the input files are ready to run GLiMMPS:

Jobs for running GLiMMPS are in the glimmps directories and for each population 
in the batch_cis directories.

Basically, GLiMMPS is run by using the commands:

```bash
$ cd glimmps/CEU/batch_cis
$ ../../../scripts/sQTLregress.oneexon.cis.R ../input/Geuvadis_exon_info_SE.txt ../input/JC_Geuvadis_DP_SE.txt ../input/JC_Geuvadis_IC_SE.txt Geuvadis_allsnp chr1 1 100
```

The exons are run in batched of 100 exons, for the optimal number of jobs and performance.

The results will be in the Glimmps_each_exon_cis directory for each population.

Running permutation tests
---------------------------------------
The scripts for running the permutation tests are in the random scripts. In brief, 
we shuffle the genotype of the individuals to test the association under null hypothesis.

GLiMMPS runs on the random genotypes using the commands:
```bash
$ cd glimmps/CEU/batch_perm1
$ ../../../scripts/sQTLregress.oneexon.cis.perm.R ../input/Geuvadis_exon_info_SE.txt ../input/JC_Geuvadis_DP_SE.txt ../input/JC_Geuvadis_IC_SE.txt Geuvadis_allsnp chr1 1 100 1
```

The output will be in the Glimmps_each_exon_cis_perm1 directory for each population.

Contacts and bug reports
------------------------
Emad Bahrami-Samani
bahramise1@email.chop.edu

Copyright and License Information
---------------------------------
Copyright (C) 2019 Children's Hospital of Philadelphia
Emad Bahrami-Samani

Authors: Emad Bahrami-Samani

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.

