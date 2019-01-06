for f in {CEU,FIN,GBR,TSI,YRI}; do awk 'FILENAME==ARGV[1] {if (NR==1) split($0,individuals," "); next} FILENAME==ARGV[2] {indiv_order[$2]=$1; next} FILENAME==ARGV[3] {filtered_exons[$1]=1; next} FILENAME==ARGV[4] {kindx+=1; if (kindx==1) s="FID\tIID"; else if ($1 in filtered_exons) {iindx+=1; s=s"\tSE_"$1; split($0,c,"\t"); split(c[2],a,","); split(c[3],b,","); for (i=1;i<=223;i+=1) IC[i,iindx]=a[i]; for (i=1;i<=223;i+=1) SC[i,iindx]=b[i]; split(c[4],a2,","); split(c[5],b2,","); for (i=1;i<=222;i+=1) IC[i+223,iindx]=a2[i]; for (i=1;i<=222;i+=1) SC[i+223,iindx]=b2[i];}} END{print s; for (i=2;i<=length(individuals);i+=1) {s=individuals[i]"\t1"; for (j=1; j<=iindx; j+=1) s=s"\t"IC[indiv_order[individuals[i]],j]+SC[indiv_order[individuals[i]],j]; print s;}}' ${f}/genotype/header ../splicing/bam_list.txt ${f}/input/filtered_exons.txt ../splicing/output/JC.raw.input.SE.txt > ${f}/input/JC_Geuvadis_DP_SE.txt; done

for f in {CEU,FIN,GBR,TSI,YRI}; do awk 'FILENAME==ARGV[1] {if (NR==1) split($0,individuals," "); next} FILENAME==ARGV[2] {indiv_order[$2]=$1; next} FILENAME==ARGV[3] {filtered_exons[$1]=1; next} FILENAME==ARGV[4] {kindx+=1; if (kindx==1) s="FID\tIID"; else if ($1 in filtered_exons) {iindx+=1; s=s"\tSE_"$1; split($0,c,"\t"); split(c[2],a,","); split(c[3],b,","); for (i=1;i<=223;i+=1) IC[i,iindx]=a[i]; for (i=1;i<=223;i+=1) SC[i,iindx]=b[i]; split(c[4],a2,","); split(c[5],b2,","); for (i=1;i<=222;i+=1) IC[i+223,iindx]=a2[i]; for (i=1;i<=222;i+=1) SC[i+223,iindx]=b2[i];}} END{print s; for (i=2;i<=length(individuals);i+=1) {s=individuals[i]"\t1"; for (j=1; j<=iindx; j+=1) s=s"\t"IC[indiv_order[individuals[i]],j]; print s;}}' ${f}/genotype/header ../splicing/bam_list.txt ${f}/input/filtered_exons.txt ../splicing/output/JC.raw.input.SE.txt > ${f}/input/JC_Geuvadis_IC_SE.txt; done
