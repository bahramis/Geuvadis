import sys, os

for line in open(sys.argv[1]) :

  line = line.strip()
  if line == "" : continue
  fn = line.split("\t")[0]
  fpop = line.split("\t")[1]

  home_dir = "/home/bahramise1/scratch/Geuvadis_2/"
  genotype_dir =  os.path.join(home_dir, "glimmps", fpop, "genotype/")

  qsub_fn = os.path.join(home_dir, "job_" + fpop+ "_" + fn + "_genotype.qsub")

  fh = open(qsub_fn, "w")

  fh.write("python scripts/transpose_ped.py " + genotype_dir + fn + ".ped " + genotype_dir + fn + "_tposed.ped\n")
  fh.write("mv " + genotype_dir + fn + "_tposed.ped " + genotype_dir + fn + ".ped\n")
  fh.write("cd " + genotype_dir + "\n")
  fh.write("plink --file " + fn + " --noweb --recodeA --out " + fn + "_recode --no-parents --no-pheno\n")
  fh.write("cd ..\n")
  fh.write("python scripts/transpose_plink_raw.py " + genotype_dir + fn + "_recode.raw " + genotype_dir + fn + "_recode_tposed.raw\n")
  fh.write("mv " + genotype_dir + fn + "_recode_tposed.raw " + genotype_dir + fn + "_tposed.raw\n")

  print fn
  

