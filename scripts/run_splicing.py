import sys, os

home_dir = "/home/bahramise1/scratch/Geuvadis_2/"
splicing_dir =  os.path.join(home_dir, "splicing/")
bam_list_file =  open(os.path.join(home_dir, "splicing/bam_list.txt"), "w")
b1_file =  open(os.path.join(home_dir, "splicing/b_1.txt"), "w")
b2_file =  open(os.path.join(home_dir, "splicing/b_2.txt"), "w")

mapping_dir =  os.path.join(home_dir, "mapping/")
qsub_fn = os.path.join(home_dir, "run_splicing.qsub")

i = 0
for line in open(sys.argv[1]) :
  i += 1
  line = line.strip()
  if line == "" : continue 
  fn = line.split("\t")[0]

  if (i == 1) :
    b1_file.write(mapping_dir + fn + ".bam")
  elif (i <= 223) :
    b1_file.write("," + mapping_dir + fn + ".bam")
  elif (i == 224) :
    b2_file.write(mapping_dir + fn + ".bam")
  else :
    b2_file.write("," + mapping_dir + fn + ".bam")

  bam_list_file.write(str(i) + "\t" + fn + "\n")
  print fn

b1_file.write("\n")
b2_file.write("\n")

bam_list_file.close()
b1_file.close()
b2_file.close()

fh = open(qsub_fn, "w")
fh.write("cd " + splicing_dir + "\n")
fh.write("python /home/bahramise1/programs/rmats_pipeline/rmats.py --b1 b_1.txt --b2 b_2.txt --gtf ../genome_dir/Homo_sapiens.GRCh37.75.gtf --task prep --od output -t paired --nthread 60 --readLength 75 --anchorLength 1 --tmp temp --statoff\n")
fh.write("python /home/bahramise1/programs/rmats_pipeline/rmats.py --b1 b_1.txt --b2 b_2.txt --gtf ../genome_dir/Homo_sapiens.GRCh37.75.gtf --task post --od output -t paired --nthread 60 --readLength 75 --anchorLength 1 --tmp temp --statoff\n")

fh.write("python ../scripts/parse_rmats.py ../metadata.txt b_1.txt b_2.txt output/SE.MATS.JC.txt > filtered_exons.txt\n")

fh.close()

