import sys, os

for line in open(sys.argv[1]) :

  line = line.strip()
  if line == "" : continue 
  f = line.split("\t")[0]
  popu = line.split("\t")[1]

  base_fn = os.path.splitext(os.path.basename(f))[0]
  qsub_dir = "/home/bahramise1/scratch/Geuvadis/genotype_parallel"
  qsub_fn = os.path.join(qsub_dir, "job_" + base_fn + "_genotype.qsub")

  fh = open(qsub_fn, "w")
  fh.write("cd "+qsub_dir+"\n\n")

  fh.write("python transpose_ped.py " + popu + "/" + base_fn + ".ped " + popu + "/" + base_fn + "_tposed.ped\n")
  fh.write("mv " + popu + "/" + base_fn + "_tposed.ped " + popu + "/" + base_fn + ".ped\n")
  fh.write("cd " + popu + "/\n")
  fh.write("plink --file " + base_fn + " --noweb --recodeA --out " + base_fn + "_recode --no-parents --no-pheno\n")
  fh.write("cd ..\n")
  fh.write("python transpose_plink_raw.py " + popu + "/" + base_fn + "_recode.raw " + popu + "/" + base_fn + "_recode_tposed.raw\n")
  fh.write("mv " + popu + "/" + base_fn + "_recode_tposed.raw " + popu + "/" + base_fn + "_tposed.raw\n")

  print f
  

