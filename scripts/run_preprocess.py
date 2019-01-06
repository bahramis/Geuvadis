import sys, os

for line in open(sys.argv[1]) :

  line = line.strip()
  if line == "" : continue 
  fn = line.split("\t")[0]
  fid = line.split("\t")[1]
  fpop = line.split("\t")[2]
  fweb_1 = line.split("\t")[3]
  fweb_2 = line.split("\t")[4] 

  home_dir = "/home/bahramise1/scratch/Geuvadis_2/"
  input_dir = os.path.join(home_dir, "input/")
  mapping_dir =  os.path.join(home_dir, "mapping/")

  qsub_fn = os.path.join(home_dir, "job_" + fn + "_pipeline.qsub")

  fh = open(qsub_fn, "w")

  fh.write("cd " + input_dir + "\n")
  fh.write("wget " + fweb_1 + "\n")
  fh.write("gzip -d " + fid + "_1.fastq.gz\n")
  fh.write("mv " + fid + "_1.fastq " + fn + "_1.fastq\n")
  fh.write("awk '{if (NR%4==1 || NR%4==3) print $0; else print substr($0,1,75);}' " + fn + "_1.fastq > " + fn + "_1_temp.fastq\n")
  fh.write("mv " + fn + "_1_temp.fastq " + fn + "_1.fastq\n")
  fh.write("gzip " + fn + "_1.fastq\n")

  fh.write("wget " + fweb_2 + "\n")
  fh.write("gzip -d " + fid + "_2.fastq.gz\n")
  fh.write("mv " + fid + "_2.fastq " + fn + "_2.fastq\n")
  fh.write("awk '{if (NR%4==1 || NR%4==3) print $0; else print substr($0,1,75);}' " + fn + "_2.fastq > " + fn + "_2_temp.fastq\n")
  fh.write("mv " + fn + "_2_temp.fastq " + fn + "_2.fastq\n")
  fh.write("gzip " + fn + "_2.fastq\n")


  fh.write("cd " + mapping_dir + "\n")
  fh.write("mkdir " + fn + "\n")
  fh.write("cd " + fn + "\n")
  fh.write("STAR --readFilesCommand zcat --chimSegmentMin 2 --outFilterMismatchNmax 3 --alignEndsType EndToEnd --runThreadN 4 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 6 --alignIntronMax 299999 --genomeDir ../../genome_dir/ENSEMBL.homo_sapiens.release-75 --sjdbGTFfile ../../genome_dir/Homo_sapiens.GRCh37.75.gtf --readFilesIn ../../input/" + fn + "_1.fastq.gz ../../input/" + fn + "_2.fastq.gz\n")
  fh.write("cd " + mapping_dir + "\n")
  fh.write("mv " + fn + "/Aligned.sortedByCoord.out.bam " + fn + ".bam\n")
  fh.write("samtools index " + fn + ".bam\n")

  print fn
  

