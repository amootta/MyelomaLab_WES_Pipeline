REF=/home/regma03/Scratch/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
P=COS_001_063

cat ${P}_PRG.tsv ${P}_BSL.tsv | cut -f1-4 | sort -k1,1 -k2,2n -u > ${P}.union_sites.tsv
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' ${P}.union_sites.tsv | sort -k1,1 -k2,2n -u > ${P}.union.bed

for T in PRG BSL; do
  bcftools mpileup -f "$REF" -R ${P}.union.bed -a AD -I \
      -q 20 -Q 20 -d 10000 /myriadfs/home/regma03/Scratch/WES_tumor_batch1/${P}_${T}.recal.cram \
  | bcftools norm -m -any -f "$REF" \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD{0}\t%AD{1}]\n' \
  > ${T}.pileup.tsv
done

cut -f1,2 ${P}.union_sites.tsv > ${P}.loci.txt

for T in PRG BSL; do
  alleleCounter -l ${P}.loci.txt -r "$REF" -m 20 -q 20 \
    -b /myriadfs/home/regma03/Scratch/WES_tumor_batch1/${P}_${T}.recal.cram \
    -o ${P}_${T}.allelecount.txt
done

for T in PRG BSL; do
  awk 'BEGIN{OFS="\t"; b["A"]=3; b["C"]=4; b["G"]=5; b["T"]=6}
       NR==FNR {ref[$1":"$2]=$3; alt[$1":"$2]=$4; next}
       /^#/ {next}
       {k=$1":"$2; if (k in ref) print $1,$2,ref[k],alt[k],$(b[ref[k]]),$(b[alt[k]])}' \
      ${P}.union_sites.tsv ${P}_${T}.allelecount.txt > ./counts/${P}_${T}.forced.tsv
done